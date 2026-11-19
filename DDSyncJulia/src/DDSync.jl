module DDSync
# DDSync.jl — Julia implementation of DDSync standalone (MATLAB) package
#
# Design goals
#   - Minimal dependencies: only Julia standard libraries.
#   - Memory-safe for huge dt.cc via two-pass streaming + per-(station,phase) spooling.
#   - Output file formats intended to match the MATLAB DDSync standalone package:
#       * dt_sync.cc (HypoDD style)
#       * theta/theta_STA_PH.txt
#       * thetastd/std_theta_STA_PH.txt
#       * sync_metrics.txt
#
# Notes
#   - This is a direct, readable implementation. It is not tuned for peak speed.
#   - For very large catalogs, consider cfg[:runtime][:store_dense_group_arrays] = false
#     to avoid nEvents × nGroups dense storage in memory.

using SparseArrays
using LinearAlgebra
using Statistics
using Random
using Printf

include("EventIDs.jl")
include("Uncertainty.jl")
include("Metrics.jl")

# ---------------------------- Public API ----------------------------

"""
    config_default() -> Dict{Symbol,Any}

Return a nested configuration dictionary whose fields mirror the MATLAB standalone
package `+ddsync/config_default.m`. Julia uses separate `irls` and `runtime` sections.
"""
function config_default()
    cfg = Dict{Symbol,Any}()

    cfg[:io] = Dict{Symbol,Any}(
        :infile_dt     => "dt.cc",
        :catalog_file  => "catalog.txt",
        :tmpdir        => "ddsync_tmp",
        :thetadir      => "theta",
        :thetastd_dir  => "thetastd",
        :metrics_file  => "sync_metrics.txt",
        :out_dt_sync   => "dt_sync.cc",
    )

    cfg[:weights] = Dict{Symbol,Any}(
        :base_mode => "cc",   # "cc" | "cc2" | "ones"
    )

    cfg[:robust] = Dict{Symbol,Any}(
        :K_SIGMA => 20.0,
        :MIN_EDGES => 30,
        :min_scale => 5e-4,
    )

    cfg[:irls] = Dict{Symbol,Any}(
        :iters   => 10,
        :C_HUBER => 1.345,
        :rel_tol => 1e-3,
    )

    cfg[:numeric] = Dict{Symbol,Any}(
        :RIDGE_EPS => 1e-10,
    )

    cfg[:output] = Dict{Symbol,Any}(
        :write_dt_sync        => true,
        :weight_mode          => "thetaStd", # "base" | "robust" | "combined" | "thetaStd"
        :thetastd_scale_mode  => "fixed",    # "fixed" | "median"
        :thetastd_scale_fixed => 500.0,
        :thetastd_weight_cap  => 1.0,
        :dt_decimals          => 5,
        :dt_weight_decimals   => 4,
        :theta_decimals       => 9,
        :thetastd_decimals    => 9,
        :station_field_width  => 8,
        :dt_field_width       => 10,
        :weight_field_width   => 8,
        :write_pruned_edges   => false,
    )

    cfg[:std] = Dict{Symbol,Any}(
        :export                           => true,
        :mode                             => "auto",          # exact for small components
        :auto_probe                       => "hutch",         # sampler for larger components in auto mode
        :fallback                         => "pseudo_degree", # only entries still unusable after selected exact repair
        :probes                           => 4096,           # maximum; stop early when diagnostics pass
        :min_probes                       => 128,
        :rel_tol                          => 0.10,           # maximum per-node MC standard error / diagonal
        :exact_max_nred                   => 2048,
        :repair_max_nodes                 => 128,            # selected full-system solves; reuses Hutch factorization
        :unconverged                      => "keep",         # retain positive estimates with relative MC SE < 1; or "fallback"
        :probe_dist                       => "rademacher",    # "rademacher" | "gaussian"
        :batch                            => 128,
        :report_every_batch               => 10,
        :max_nred                         => 40000,
        # sigma_hat measures scatter of (theta_i - theta_j) - measured_dt_ij after fitting.
        # It is estimated separately for each connected station-phase component.
        :min_sigma                        => 5e-4,           # minimum assumed mismatch scatter (s), not minimum std_theta
        :apply_min_sigma                  => true,           # use max(sigma_hat, min_sigma); false uses sigma_hat
        :min_diag_rel                     => 1e-12,          # deprecated; ignored (no variance clipping)
        :pseudo_weight_source             => "combined",      # "base" | "robust" | "combined"
        :pseudo_weight_eps                => 1e-6,
        :thetastd_write_weightcol         => true,
        :thetastd_write_alt_node_weightcol => false,
    )

    cfg[:runtime] = Dict{Symbol,Any}(
        :prog_every_lines        => 2_000_000,
        :print_new_groups        => false,
        # If true, store per-group theta/std as dense Float32 vectors length nEvents (fast, memory heavy).
        # If false, store per-group theta/std as Dict{Int,Float32} (memory light, slower in pass 2).
        :store_dense_group_arrays => false,
        :decision_chunk_records   => 200_000,
        :median_scale_cap         => 2_000_000,  # max kept obs sampled when computing median thetaStd scale
    )

    return cfg
end

"""
    run(cfg::Dict=config_default())

Run DDSync using the provided configuration dictionary. All outputs are written to disk.
"""
function run(cfg::Dict{Symbol,Any}=config_default())
    # ---- unpack config ----
    io_cfg     = cfg[:io]
    w_cfg      = cfg[:weights]
    robust_cfg = cfg[:robust]
    irls_cfg   = cfg[:irls]
    num_cfg    = cfg[:numeric]
    out_cfg    = cfg[:output]
    std_cfg    = cfg[:std]
    rt_cfg     = cfg[:runtime]
    validateStdConfig(std_cfg)
    isfinite(num_cfg[:RIDGE_EPS]) && num_cfg[:RIDGE_EPS] >= 0 || error("numeric.RIDGE_EPS must be finite and nonnegative")

    infile_dt    = String(io_cfg[:infile_dt])
    catalog_file = String(io_cfg[:catalog_file])
    tmpdir       = String(io_cfg[:tmpdir])
    thetadir     = String(io_cfg[:thetadir])
    thetastd_dir = String(io_cfg[:thetastd_dir])
    metrics_file = String(io_cfg[:metrics_file])
    out_dt_sync  = String(io_cfg[:out_dt_sync])

    mkpath(tmpdir)
    mkpath(thetadir)
    if Bool(std_cfg[:export])
        mkpath(thetastd_dir)
    end

    @printf("Reading event IDs from %s...\n", catalog_file)
    event_index = readEventIndex(catalog_file)
    maxEventID = length(event_index.original)
    index_paths = indexOutputPaths(out_dt_sync)
    write_dt = Bool(get(out_cfg, :write_dt_sync, true))
    mkpath(dirname(abspath(out_dt_sync)))
    mkpath(dirname(abspath(metrics_file)))
    destinations = String[metrics_file]
    write_dt && push!(destinations, out_dt_sync)
    if event_index.remapped
        append!(destinations, [index_paths.catalog, index_paths.mapping])
        write_dt && push!(destinations, index_paths.dt)
        println("Event IDs are not the complete range 1:N. Mapping $maxEventID events to 1:N in catalog row order.")
        println("  theta/thetastd will use the new IDs; dt_sync retains the original IDs.")
        println("  GraphSplit catalog: $(index_paths.catalog)")
        println("  ID correspondence: $(index_paths.mapping)")
    else
        println("Event IDs already form 1:$maxEventID; preserving the existing IDs.")
    end
    checkIndexOutputPaths(destinations, [infile_dt, catalog_file])

    @printf("Pass 1: streaming %s and spooling groups...\n", infile_dt)
    groupKeys, groupFiles, groupCounts = spoolGroups(infile_dt, tmpdir;
        event_index = event_index,
        prog_every = Int(rt_cfg[:prog_every_lines]),
        print_new_groups = Bool(rt_cfg[:print_new_groups])
    )
    nGroups = length(groupKeys)
    event_index.remapped && writeIndexFiles(catalog_file, event_index, index_paths)
    @printf("Found %d station+phase groups.\n", nGroups)

    key2g = Dict{String,Int}()
    for (g,key) in enumerate(groupKeys)
        key2g[key] = g
    end

    # Base weight function
    base_mode = lowercase(String(w_cfg[:base_mode]))
    weight_fun_base = makeBaseWeightFun(base_mode)

    # Robust params
    K_SIGMA   = Float64(robust_cfg[:K_SIGMA])
    MIN_EDGES = Int(robust_cfg[:MIN_EDGES])
    ROBUST_MIN_SCALE = Float64(robust_cfg[:min_scale])

    # IRLS
    IRLS_ITERS   = Int(irls_cfg[:iters])
    C_HUBER      = Float64(irls_cfg[:C_HUBER])
    IRLS_REL_TOL = Float64(irls_cfg[:rel_tol])

    RIDGE_EPS = Float64(num_cfg[:RIDGE_EPS])

    # Output
    WRITE_DT_SYNC        = Bool(get(out_cfg, :write_dt_sync, true))
    OUTPUT_WEIGHT_MODE   = lowercase(String(out_cfg[:weight_mode]))
    THETASTD_SCALE_MODE  = lowercase(String(out_cfg[:thetastd_scale_mode]))
    THETASTD_SCALE_FIXED = Float64(out_cfg[:thetastd_scale_fixed])
    THETASTD_WEIGHT_CAP  = Float64(out_cfg[:thetastd_weight_cap])
    THETA_DECIMALS       = max(0, Int(round(Float64(get(out_cfg, :theta_decimals, 9)))))
    THETASTD_DECIMALS    = max(0, Int(round(Float64(get(out_cfg, :thetastd_decimals, 9)))))
    WRITE_PRUNED_EDGES   = Bool(out_cfg[:write_pruned_edges])

    # Std
    EXPORT_THETA_STD = Bool(std_cfg[:export])
    THETASTD_MODE     = lowercase(String(std_cfg[:mode]))
    THETASTD_FALLBACK = lowercase(String(std_cfg[:fallback]))
    STD_PROBES        = Int(std_cfg[:probes])
    STD_PROBE_DIST    = lowercase(String(std_cfg[:probe_dist]))
    STD_BATCH         = Int(std_cfg[:batch])
    STD_REPORT_EVERY  = Int(std_cfg[:report_every_batch])
    STD_MAX_NRED      = Int(std_cfg[:max_nred])
    STD_MIN_SIGMA     = Float64(std_cfg[:min_sigma])
    STD_APPLY_MIN_SIGMA = Bool(std_cfg[:apply_min_sigma]) && (STD_MIN_SIGMA > 0.0)
    STD_MIN_DIAGREL   = Float64(std_cfg[:min_diag_rel])
    PSEUDO_WEIGHT_SOURCE = lowercase(String(std_cfg[:pseudo_weight_source]))
    PSEUDO_WEIGHT_EPS    = Float64(std_cfg[:pseudo_weight_eps])
    WRITE_THETASTD_WEIGHTCOL = Bool(std_cfg[:thetastd_write_weightcol])
    WRITE_THETASTD_ALT_NODEW_COL = Bool(std_cfg[:thetastd_write_alt_node_weightcol])
    # Runtime
    store_dense = Bool(rt_cfg[:store_dense_group_arrays])
    decision_chunk = Int(rt_cfg[:decision_chunk_records])
    median_cap = Int(rt_cfg[:median_scale_cap])

    # Metrics file
    fidm = open(metrics_file, "w")
    @printf(fidm, "# DDSync metrics
")
    writeIndexMetrics(fidm, event_index, index_paths, WRITE_DT_SYNC)
    mode_print = String(out_cfg[:weight_mode])
    dist_print = String(std_cfg[:probe_dist])
    thetastd_mode_print = String(std_cfg[:mode])
    @printf(fidm, "# K_SIGMA=%g MIN_EDGES=%d IRLS_ITERS=%d C_HUBER=%g MODE=%s EXPORT_THETA_STD=%d STD_PROBES=%d DIST=%s THETASTD_MODE=%s
",
        K_SIGMA, MIN_EDGES, IRLS_ITERS, C_HUBER, mode_print,
        EXPORT_THETA_STD ? 1 : 0, STD_PROBES, dist_print, thetastd_mode_print)
    @printf(fidm, "# STD_MIN_PROBES=%d STD_REL_TOL=%g STD_BATCH=%d STD_FALLBACK=%s STD_EXACT_MAX_NRED=%d\n",
        std_cfg[:min_probes], std_cfg[:rel_tol], STD_BATCH, THETASTD_FALLBACK, std_cfg[:exact_max_nred])
    println(fidm, "# STD_AUTO_PROBE=",get(std_cfg,:auto_probe,"hutch"))
    @printf(fidm, "# STD_REPAIR_MAX_NODES=%d STD_UNCONVERGED=%s\n", std_cfg[:repair_max_nodes], std_cfg[:unconverged])
    @printf(fidm, "# STD_MIN_SIGMA=%g STD_APPLY_MIN_SIGMA=%s\n",STD_MIN_SIGMA,string(STD_APPLY_MIN_SIGMA))
    metrics_rows = NamedTuple[]
    metric_components = NamedTuple[]
    flush(fidm)

    overall_edges = 0
    overall_kept  = 0
    nPseudoGroups = 0
    nPseudoComps  = 0

    # Per-group storage for pass 2
    decFiles = Vector{String}(undef, nGroups)
    theta_store = Vector{Any}(undef, nGroups)  # Dict{Int,Float32} or Vector{Float32}
    std_store   = Vector{Any}(undef, nGroups)  # Dict{Int,Float32} or Vector{Float32} or nothing

    for g in 1:nGroups
        key = groupKeys[g]
        n_edges = groupCounts[g]
        @printf("[%d/%d] Processing group %s (edges=%d)\n", g, nGroups, key, n_edges)

        X = readGroupFile(groupFiles[g])  # 4 x nEdges Float32
        gi = Int.(round.(X[1, :]))
        gj = Int.(round.(X[2, :]))
        gd = Float64.(X[3, :])
        gc = Float64.(X[4, :])

        overall_edges += length(gd)

        w_base = similar(gd)
        @inbounds for k in eachindex(gc)
            cc = gc[k]
            wb = (isfinite(cc) ? weight_fun_base(cc) : 0.0)
            w_base[k] = (isfinite(wb) && wb > 0.0) ? wb : 0.0
        end

        (ev_all, theta_local, ref_local, std_local, deg_local, nodew_local,
         keep_mask, wrob, sigma_hat, pseudo_info) = processStationPhaseGroup(
            gi, gj, gd, w_base, maxEventID;
            K_SIGMA=K_SIGMA, MIN_EDGES=MIN_EDGES, ROBUST_MIN_SCALE=ROBUST_MIN_SCALE,
            IRLS_ITERS=IRLS_ITERS, C_HUBER=C_HUBER, IRLS_REL_TOL=IRLS_REL_TOL,
            RIDGE_EPS=RIDGE_EPS,
            EXPORT_THETA_STD=EXPORT_THETA_STD,
            THETASTD_MODE=THETASTD_MODE, THETASTD_FALLBACK=THETASTD_FALLBACK,
            STD_PROBES=STD_PROBES, STD_PROBE_DIST=STD_PROBE_DIST,
            STD_BATCH=STD_BATCH, STD_REPORT_EVERY=STD_REPORT_EVERY,
            STD_MAX_NRED=STD_MAX_NRED, STD_MIN_SIGMA=STD_MIN_SIGMA, STD_APPLY_MIN_SIGMA=STD_APPLY_MIN_SIGMA, STD_MIN_DIAGREL=STD_MIN_DIAGREL,
            PSEUDO_WEIGHT_SOURCE=PSEUDO_WEIGHT_SOURCE, PSEUDO_WEIGHT_EPS=PSEUDO_WEIGHT_EPS,
            THETASTD_SCALE_FIXED=THETASTD_SCALE_FIXED,
            WRITE_ALT_NODEW_COL=WRITE_THETASTD_ALT_NODEW_COL,
            STD_AUTO_PROBE=lowercase(String(get(std_cfg,:auto_probe,"hutch"))),
            STD_MIN_PROBES=Int(std_cfg[:min_probes]), STD_REL_TOL=Float64(std_cfg[:rel_tol]),
            STD_EXACT_MAX_NRED=Int(std_cfg[:exact_max_nred]),
            STD_REPAIR_MAX_NODES=Int(std_cfg[:repair_max_nodes]), STD_UNCONVERGED=lowercase(String(std_cfg[:unconverged])),
        )

        for info in pseudo_info.components
            # Keep only scalar diagnostics; do not retain per-node arrays for reporting.
            compact = (; (k=>getproperty(info,k) for k in propertynames(info)
                if !(k in (:event_ids,:node_methods,:node_rel_se)))...)
            push!(metric_components, merge(compact,(group=key,)))
        end
        n_kept = count(identity, keep_mask)
        overall_kept += n_kept

        if pseudo_info.used_pseudo
            nPseudoGroups += 1
            nPseudoComps  += pseudo_info.n_pseudo_comps
        end

        # Store theta/std for pass 2
        if store_dense
            theta_dense = fill(Float32(NaN), maxEventID)
            @inbounds for (eid, val) in zip(ev_all, theta_local)
                if 1 <= eid <= maxEventID && isfinite(val)
                    theta_dense[eid] = Float32(val)
                end
            end
            theta_store[g] = theta_dense

            if EXPORT_THETA_STD
                std_dense = fill(Float32(NaN), maxEventID)
                @inbounds for (eid, val) in zip(ev_all, std_local)
                    if 1 <= eid <= maxEventID && isfinite(val)
                        std_dense[eid] = Float32(val)
                    end
                end
                std_store[g] = std_dense
            else
                std_store[g] = nothing
            end
        else
            theta_dict = Dict{Int,Float32}()
            @inbounds for (eid, val) in zip(ev_all, theta_local)
                if isfinite(val)
                    theta_dict[eid] = Float32(val)
                end
            end
            theta_store[g] = theta_dict

            if EXPORT_THETA_STD
                std_dict = Dict{Int,Float32}()
                @inbounds for (eid, val) in zip(ev_all, std_local)
                    if isfinite(val)
                        std_dict[eid] = Float32(val)
                    end
                end
                std_store[g] = std_dict
            else
                std_store[g] = nothing
            end
        end

        # Write theta + std files for internal IDs 1:N (one row per catalog event).
        sta, ph = splitKey(key)

        theta_fn = joinpath(thetadir, @sprintf("theta_%s_%s.txt", sta, ph))
        writeThetaTriples(theta_fn, maxEventID, ev_all, theta_local, ref_local; theta_decimals=THETA_DECIMALS)

        if EXPORT_THETA_STD
            std_fn = joinpath(thetastd_dir, @sprintf("std_theta_%s_%s.txt", sta, ph))
            writeStdQuality(joinpath(thetastd_dir, @sprintf("std_quality_%s_%s.txt",sta,ph)), pseudo_info.components)
            writeThetaStdWithDegree(std_fn, maxEventID,
                ev_all, std_local, ref_local, deg_local, nodew_local;
                write_weightcol=WRITE_THETASTD_WEIGHTCOL,
                write_alt_nodew_col=WRITE_THETASTD_ALT_NODEW_COL,
                thetastd_scale_fixed=THETASTD_SCALE_FIXED,
                thetastd_weight_cap=THETASTD_WEIGHT_CAP,
                thetastd_decimals=THETASTD_DECIMALS
            )
        end

        # Decisions file (edge order == spool order)
        dec_fn = joinpath(tmpdir, @sprintf("dec_%s.bin", sanitizeKey(key)))
        writeDecisionFile(dec_fn, keep_mask, wrob)
        decFiles[g] = dec_fn

        # Metrics row
        pruned_pct = 100.0 * (length(gd) - n_kept) / max(length(gd), 1)
        if n_kept > 0
            rob_kept = wrob[keep_mask]
            r50 = median(rob_kept)
            r05 = prctile(rob_kept, 5.0)
            r95 = prctile(rob_kept, 95.0)
            frac_dw = mean(rob_kept .< 0.99)
        else
            r50 = NaN; r05 = NaN; r95 = NaN; frac_dw = NaN
        end

        push!(metrics_rows,(group=key,n_edges=length(gd),n_kept=n_kept,pruned_pct=pruned_pct,
            r50=r50,r05=r05,r95=r95,frac_dw=frac_dw,sigma_hat=sigma_hat))

        @printf("  kept=%d (%.2f%% pruned), sigma_hat=%.4g s, robust_w p05/p50/p95=%.3g/%.3g/%.3g\n",
            n_kept, pruned_pct, sigma_hat, r05, r50, r95)
    end

    overall_pruned_pct = 100.0 * (overall_edges - overall_kept) / max(overall_edges, 1)

    writeMetricsReport(fidm,metrics_rows,metric_components,
        (edges=overall_edges,kept=overall_kept,pruned_pct=overall_pruned_pct,
         pseudo_groups=nPseudoGroups,pseudo_components=nPseudoComps))
    close(fidm)
    println()
    writeUncertaintySummary(stdout,metric_components;prefix="")
    by_group=uncertaintyByGroup(metric_components)
    for g in metrics_rows
        c=get(by_group,g.group,uncertaintyCounts(()))
        if needsStdReview(c)
            @printf("[STD GROUP ADVICE] %s: %s require review. See the group recommendations in %s.\n",
                g.group,stdPercent(c.review,c.total),metrics_file)
        end
    end

    # Optional median scaling for thetaStd weights
    thetaStdScale = THETASTD_SCALE_FIXED
    if OUTPUT_WEIGHT_MODE == "thetastd" && THETASTD_SCALE_MODE == "median"
        @printf("Computing global median scale for thetaStd weights (extra streaming pass)...\n")
        thetaStdScale = computeThetaStdMedianScale(infile_dt, key2g, decFiles, std_store, store_dense, maxEventID;
            cap=median_cap, chunk=decision_chunk, event_index=event_index)
        if !(isfinite(thetaStdScale) && thetaStdScale > 0)
            thetaStdScale = THETASTD_SCALE_FIXED
        end
        @printf("thetaStdScale (median of 1/std_dt) = %.6g\n", thetaStdScale)
    end

    # Pass 2: write dt_sync.cc
    if WRITE_DT_SYNC
        @printf("Pass 2: writing %s ...\n", out_dt_sync)
    readers = openDecisionReaders(decFiles; chunk=decision_chunk)

    fin  = open(infile_dt, "r")
    fout = open(out_dt_sync, "w")
    fseq = event_index.remapped ? open(index_paths.dt, "w") : nothing
    outputs = fseq === nothing ? (fout,) : (fout, fseq)

    # Formatting (match MATLAB standalone)
    STA_W = Int(out_cfg[:station_field_width])
    DT_W  = Int(out_cfg[:dt_field_width])
    W_W   = Int(out_cfg[:weight_field_width])
    DT_DEC = Int(out_cfg[:dt_decimals])
    W_DEC  = Int(out_cfg[:dt_weight_decimals])

    pendingHeader = ""
    pendingPrinted = false
    cur_i = 0
    cur_j = 0

    try
    for (line_number, L) in enumerate(eachline(fin))
        Lt = strip(L)
        if isempty(Lt)
            continue
        end

        if startswith(Lt, "#")
            pendingHeader = L
            pendingPrinted = false

            cur_i, cur_j = mappedPair(Lt, event_index, "$infile_dt:$line_number")
            continue
        end

        toks = split(Lt)
        if length(toks) < 4
            continue
        end
        sta = toks[1]
        dt_in = tryparse(Float64, toks[2])
        cc_in = tryparse(Float64, toks[3])
        ph  = toks[4]
        if dt_in === nothing || cc_in === nothing
            continue
        end

        key = string(sta, "_", ph)

        # Unknown group: pass through using base weight and original dt
        if !haskey(key2g, key)
            if !pendingPrinted && !isempty(pendingHeader)
                println(fout, pendingHeader)
                fseq === nothing || println(fseq, sequentialHeader(pendingHeader, cur_i, cur_j))
                pendingPrinted = true
            end
            dt_corr = dt_in
            w_out = max(0.0, isfinite(cc_in) ? weight_fun_base(cc_in) : 0.0)

            for output in outputs
                @printf(output, "    %-*s %*.*f %*.*f %s\n",
                    STA_W, sta, DT_W, DT_DEC, dt_corr, W_W, W_DEC, w_out, ph)
            end
            continue
        end

        g = key2g[key]

        keep, wrob = readNextDecision!(readers[g])
        if !keep && !WRITE_PRUNED_EDGES
            continue
        end

        # dt correction
        dt_corr = dt_in
        if keep
            tstore = theta_store[g]
            ti = getStoreValue(tstore, cur_i, store_dense)
            tj = getStoreValue(tstore, cur_j, store_dense)
            if isfinite(ti) && isfinite(tj)
                dt_corr = ti - tj
            end
        end

        # output weight
        w_out = 0.0
        if keep
            if OUTPUT_WEIGHT_MODE == "base"
                w_out = max(0.0, isfinite(cc_in) ? weight_fun_base(cc_in) : 0.0)
            elseif OUTPUT_WEIGHT_MODE == "robust"
                w_out = wrob
            elseif OUTPUT_WEIGHT_MODE == "combined"
                w_out = max(0.0, isfinite(cc_in) ? weight_fun_base(cc_in) : 0.0) * wrob
            elseif OUTPUT_WEIGHT_MODE == "thetastd"
                sstore = std_store[g]
                si = getStoreValue(sstore, cur_i, store_dense)
                sj = getStoreValue(sstore, cur_j, store_dense)
                std_dt = hypot(si, sj)
                if std_dt <= 0.0 || !isfinite(std_dt)
                    w_out = 0.0
                else
                    w_out = (1.0 / std_dt) / thetaStdScale
                    if isfinite(THETASTD_WEIGHT_CAP) && THETASTD_WEIGHT_CAP > 0.0
                        w_out = min(w_out, THETASTD_WEIGHT_CAP)
                    end
                end
            else
                # default
                w_out = max(0.0, isfinite(cc_in) ? weight_fun_base(cc_in) : 0.0) * wrob
            end
        else
            # pruned edge written: weight 0
            w_out = 0.0
        end


        if !pendingPrinted && !isempty(pendingHeader)
            println(fout, pendingHeader)
            fseq === nothing || println(fseq, sequentialHeader(pendingHeader, cur_i, cur_j))
            pendingPrinted = true
        end

        for output in outputs
            @printf(output, "    %-*s %*.*f %*.*f %s\n",
                STA_W, sta, DT_W, DT_DEC, dt_corr, W_W, W_DEC, w_out, ph)
        end
    end

    finally
        close(fin)
        foreach(close, outputs)
        closeDecisionReaders(readers)
    end
    end

    @printf("Done.\n")
    if WRITE_DT_SYNC
        @printf("  Synced dt: %s\n", out_dt_sync)
        event_index.remapped && println("  Sequential synced dt: $(index_paths.dt)")
    end
    @printf("  Metrics:  %s\n", metrics_file)
    if EXPORT_THETA_STD
        @printf("  Theta std dir: %s/\n", thetastd_dir)
    end
end

# ---------------------------- Helper functions ----------------------------

# Base weight function
function makeBaseWeightFun(mode::String)
    if mode == "cc"
        return (cc::Float64)->cc
    elseif mode == "cc2"
        return (cc::Float64)->cc*cc
    elseif mode == "ones"
        return (cc::Float64)->1.0
    else
        error("weights.base_mode must be one of: cc, cc2, ones")
    end
end

sanitizeKey(key::String) = replace(key, r"[^A-Za-z0-9_]" => "_")

function splitKey(key::String)
    parts = split(key, "_"; limit=2)
    if length(parts) < 2
        return key, ""
    end
    return parts[1], parts[2]
end

# -------- Spooling dt.cc --------

"""
    spoolGroups(infile_dt, tmpdir; prog_every=2e6, print_new_groups=false)

Stream dt.cc and spool records per (station,phase) into binary group files.
Each record written as 4×Float32: [i, j, dt, cc].
"""
function spoolGroups(infile_dt::String, tmpdir::String; prog_every::Int=2_000_000,
    print_new_groups::Bool=false, event_index::Union{Nothing,EventIndex}=nothing)
    keys = String[]
    files = String[]
    counts = Int[]

    key2idx = Dict{String,Int}()
    fids = IOStream[]

    cur_i = 0
    cur_j = 0
    nLine = 0

    try
    open(infile_dt, "r") do io
        for L in eachline(io)
            nLine += 1
            if prog_every > 0 && (nLine % prog_every == 0)
                @printf("  ... streamed %d lines\n", nLine)
            end
            Lt = strip(L)
            isempty(Lt) && continue

            if startswith(Lt, "#")
                cur_i, cur_j = mappedPair(Lt, event_index, "$infile_dt:$nLine")
                continue
            end

            toks = split(Lt)
            length(toks) < 4 && continue

            sta = toks[1]
            dt  = tryparse(Float64, toks[2])
            cc  = tryparse(Float64, toks[3])
            ph  = toks[4]
            (dt === nothing || cc === nothing) && continue
            cur_i > 0 && cur_j > 0 || error("$infile_dt:$nLine: observation precedes an event-pair header.")

            key = string(sta, "_", ph)

            g = get(key2idx, key, 0)
            if g == 0
                g = length(keys) + 1
                key2idx[key] = g
                push!(keys, key)
                fn = joinpath(tmpdir, @sprintf("grp_%s.bin", sanitizeKey(key)))
                push!(files, fn)
                push!(counts, 0)
                push!(fids, open(fn, "w"))
                if print_new_groups
                    @printf("  new group: %s\n", key)
                end
            end

            # record: [i, j, dt, cc] as Float32 (MATLAB-compatible spool format)
            rec = Float32[cur_i, cur_j, dt, cc]
            write(fids[g], rec)
            counts[g] += 1
        end
    end

    finally
        foreach(close, fids)
    end

    return keys, files, counts
end

function readGroupFile(fn::String)
    nbytes = filesize(fn)
    nfloat = Int(nbytes ÷ 4)
    if nfloat % 4 != 0
        error("Group file $fn size not divisible by 4 floats.")
    end
    buf = Vector{Float32}(undef, nfloat)
    open(fn, "r") do io
        read!(io, buf)
    end
    nrec = nfloat ÷ 4
    return reshape(buf, 4, nrec)
end

# -------- Decisions file I/O --------

"""\
    writeDecisionFile(fn, keep_mask, wrob)

Write the per-edge decision stream used in pass 2.

Each record is 8 bytes: two Float32 values:
  (1) keep flag encoded as 1.0 or 0.0
  (2) robust weight factor wrob

We accept any AbstractVector{Bool} for keep_mask (e.g., BitVector)
and any real-valued vector for wrob.
"""
function writeDecisionFile(fn::String, keep_mask::AbstractVector{Bool}, wrob::AbstractVector{<:Real})
    @assert length(keep_mask) == length(wrob) "keep_mask and wrob must have the same length"
    open(fn, "w") do io
        @inbounds for k in eachindex(keep_mask)
            write(io, Float32(keep_mask[k] ? 1.0f0 : 0.0f0))
            write(io, Float32(wrob[k]))
        end
    end
end

mutable struct DecisionReader
    io::IOStream
    buf::Vector{Float32}     # length 2*chunk
    idx::Int                 # next record index within buffer (1-based)
    nrec::Int                # number of records currently in buffer
    chunk::Int               # records per fill
    total::Int               # total records in file
    read_so_far::Int         # records read so far
end

function openDecisionReaders(decFiles::Vector{String}; chunk::Int=200_000)
    readers = Vector{DecisionReader}(undef, length(decFiles))
    for (g,fn) in enumerate(decFiles)
        io = open(fn, "r")
        total = Int(filesize(fn) ÷ 8)  # 2*Float32 per record
        buf = Vector{Float32}(undef, 2*chunk)
        readers[g] = DecisionReader(io, buf, 1, 0, chunk, total, 0)
    end
    return readers
end

function closeDecisionReaders(readers::Vector{DecisionReader})
    for r in readers
        close(r.io)
    end
end

function fillDecisionBuffer!(r::DecisionReader)
    remaining = r.total - r.read_so_far
    if remaining <= 0
        r.nrec = 0
        r.idx = 1
        return
    end
    n_to = min(r.chunk, remaining)
    # read 2*n_to floats
    viewbuf = view(r.buf, 1:(2*n_to))
    read!(r.io, viewbuf)
    r.nrec = n_to
    r.idx = 1
    r.read_so_far += n_to
end

function readNextDecision!(r::DecisionReader)
    if r.idx > r.nrec
        fillDecisionBuffer!(r)
        if r.nrec == 0
            return false, 0.0
        end
    end
    k = r.idx
    keep = r.buf[2*(k-1)+1] > 0.5f0
    wrob = Float64(r.buf[2*(k-1)+2])
    r.idx += 1
    return keep, wrob
end

# -------- Store access helpers --------

@inline function getStoreValue(store, eid::Int, store_dense::Bool)
    if store === nothing
        return NaN
    end
    if store_dense
        v = (1 <= eid <= length(store)) ? store[eid] : Float32(NaN)
        return Float64(v)
    else
        return Float64(get(store::Dict{Int,Float32}, eid, Float32(NaN)))
    end
end

# -------- Median scale computation --------

function computeThetaStdMedianScale(infile_dt::String,
    key2g::Dict{String,Int},
    decFiles::Vector{String},
    std_store::Vector{Any},
    store_dense::Bool,
    maxEventID::Int;
    cap::Int=2_000_000,
    chunk::Int=200_000,
    event_index::Union{Nothing,EventIndex}=nothing)

    readers = openDecisionReaders(decFiles; chunk=chunk)
    vals = Float64[]
    sizehint!(vals, min(cap, 200_000))

    cur_i = 0
    cur_j = 0

    try
    open(infile_dt, "r") do io
        for (line_number, L) in enumerate(eachline(io))
            Lt = strip(L)
            isempty(Lt) && continue

            if startswith(Lt, "#")
                cur_i, cur_j = mappedPair(Lt, event_index, "$infile_dt:$line_number")
                continue
            end

            toks = split(Lt)
            length(toks) < 4 && continue
            (tryparse(Float64, toks[2]) === nothing || tryparse(Float64, toks[3]) === nothing) && continue
            sta = toks[1]
            ph  = toks[4]
            key = string(sta, "_", ph)
            if !haskey(key2g, key)
                continue
            end
            g = key2g[key]
            keep, _ = readNextDecision!(readers[g])
            if !keep
                continue
            end
            sstore = std_store[g]
            si = getStoreValue(sstore, cur_i, store_dense)
            sj = getStoreValue(sstore, cur_j, store_dense)

            (isfinite(si) && isfinite(sj)) || continue
            std_dt = hypot(si, sj)
            (isfinite(std_dt) && std_dt > 0.0) || continue
            v = 1.0 / std_dt
            if isfinite(v) && v > 0.0
                push!(vals, v)
            end

            if length(vals) >= cap
                break
            end
        end
    end

    finally
        closeDecisionReaders(readers)
    end

    if isempty(vals)
        return 1.0
    end
    return median(vals)
end

# -------- Output writers --------

function writeThetaTriples(fn::String, maxEventID::Int, ev_all::Vector{Int}, theta_local::Vector{Float64}, ref_local::Vector{Float64};
    theta_decimals::Int=9)

    theta_decimals = max(theta_decimals, 0)
    theta_fmt = Printf.Format("%d %." * string(theta_decimals) * "f %.0f\n")

    open(fn, "w") do io
        idx = 1
        n = length(ev_all)
        for ev in 1:maxEventID
            if idx <= n && ev_all[idx] == ev
                th = theta_local[idx]
                rf = ref_local[idx]
                idx += 1
            else
                th = NaN
                rf = NaN
            end
            print(io, Printf.format(theta_fmt, ev, th, rf))
        end
    end
end

function writeThetaStdWithDegree(fn::String, maxEventID::Int,
    ev_all::Vector{Int}, std_local::Vector{Float64}, ref_local::Vector{Float64},
    deg_local::Vector{Int}, nodew_local::Vector{Float64};
    write_weightcol::Bool=true,
    write_alt_nodew_col::Bool=false,
    thetastd_scale_fixed::Float64=500.0,
    thetastd_weight_cap::Float64=1.0,
    thetastd_decimals::Int=9)

    thetastd_decimals = max(thetastd_decimals, 0)
    fmt4 = Printf.Format("%d\t%." * string(thetastd_decimals) * "f\t%.0f\t%d\n")
    fmt5 = Printf.Format("%d\t%." * string(thetastd_decimals) * "f\t%.0f\t%d\t%.6f\n")
    fmt6 = Printf.Format("%d\t%." * string(thetastd_decimals) * "f\t%.0f\t%d\t%.6f\t%.6f\n")

    open(fn, "w") do io
        idx = 1
        n = length(ev_all)
        for ev in 1:maxEventID
            if idx <= n && ev_all[idx] == ev
                s  = std_local[idx]
                rf = ref_local[idx]
                dg = deg_local[idx]
                nw = nodew_local[idx]
                idx += 1
            else
                s  = NaN
                rf = NaN
                dg = 0
                nw = NaN
            end

            if write_weightcol
                if !(isfinite(s) && s > 0.0)
                    wtheta = 0.0
                else
                    wtheta = (1.0 / s) / thetastd_scale_fixed
                    if isfinite(thetastd_weight_cap) && thetastd_weight_cap > 0.0
                        wtheta = min(wtheta, thetastd_weight_cap)
                    end
                end
            end

            if write_weightcol && write_alt_nodew_col
                print(io, Printf.format(fmt6, ev, s, rf, dg, wtheta, nw))
            elseif write_weightcol
                print(io, Printf.format(fmt5, ev, s, rf, dg, wtheta))
            elseif write_alt_nodew_col
                print(io, Printf.format(fmt5, ev, s, rf, dg, nw))
            else
                print(io, Printf.format(fmt4, ev, s, rf, dg))
            end
        end
    end
end

# -------- Robust utilities --------

mad1(x::AbstractVector{<:Real}) = median(abs.(x .- median(x)))

function prctile(x::AbstractVector{<:Real}, p::Real)
    y = collect(Float64.(x))
    y = y[isfinite.(y)]
    if isempty(y)
        return NaN
    end
    sort!(y)
    q = clamp(p / 100.0, 0.0, 1.0)
    n = length(y)
    if n == 1
        return y[1]
    end
    h = (n - 1) * q + 1
    j = floor(Int, h)
    g = h - j
    if j >= n
        return y[n]
    end
    return (1 - g) * y[j] + g * y[j + 1]
end

huberWeight(rs::Vector{Float64}, c::Float64) = begin
    ar = abs.(rs)
    w = ones(Float64, length(ar))
    idx = ar .> c
    @inbounds for k in eachindex(ar)
        if idx[k]
            w[k] = c / max(ar[k], eps())
        end
    end
    w
end

# -------- Union-find components --------

function uf_components(n::Int, a::Vector{Int}, b::Vector{Int})
    parent = collect(1:n)
    rankv  = zeros(Int, n)

    for k in eachindex(a)
        union!(parent, rankv, a[k], b[k])
    end

    # Path compress
    for i in 1:n
        parent[i] = find_root(parent, i)
    end

    roots = parent
    u = sort(unique(roots))
    root2c = Dict{Int,Int}()
    for (i,r) in enumerate(u)
        root2c[r] = i
    end
    comps = similar(roots)
    for i in 1:n
        comps[i] = root2c[roots[i]]
    end
    return comps
end

function find_root(parent::Vector{Int}, x::Int)
    r = x
    while parent[r] != r
        r = parent[r]
    end
    while parent[x] != x
        px = parent[x]
        parent[x] = r
        x = px
    end
    return r
end

function union!(parent::Vector{Int}, rankv::Vector{Int}, x::Int, y::Int)
    rx = find_root(parent, x)
    ry = find_root(parent, y)
    rx == ry && return
    if rankv[rx] < rankv[ry]
        parent[rx] = ry
    elseif rankv[rx] > rankv[ry]
        parent[ry] = rx
    else
        parent[ry] = rx
        rankv[rx] += 1
    end
end

# -------- Core per-group processing --------

function processStationPhaseGroup(gi::Vector{Int}, gj::Vector{Int}, gd::Vector{Float64}, w_base::Vector{Float64}, maxEventID::Int;
    K_SIGMA::Float64,
    MIN_EDGES::Int,
    ROBUST_MIN_SCALE::Float64,
    IRLS_ITERS::Int,
    C_HUBER::Float64,
    IRLS_REL_TOL::Float64,
    RIDGE_EPS::Float64,
    EXPORT_THETA_STD::Bool,
    THETASTD_MODE::String,
    THETASTD_FALLBACK::String,
    STD_PROBES::Int,
    STD_PROBE_DIST::String,
    STD_BATCH::Int,
    STD_REPORT_EVERY::Int,
    STD_MAX_NRED::Int,
    STD_MIN_SIGMA::Float64,
    STD_APPLY_MIN_SIGMA::Bool,
    STD_MIN_DIAGREL::Float64,
    PSEUDO_WEIGHT_SOURCE::String,
    PSEUDO_WEIGHT_EPS::Float64,
    THETASTD_SCALE_FIXED::Float64,
    WRITE_ALT_NODEW_COL::Bool=false,
    STD_MIN_PROBES::Int=128, STD_REL_TOL::Float64=0.10, STD_EXACT_MAX_NRED::Int=2048,
    STD_REPAIR_MAX_NODES::Int=128, STD_UNCONVERGED::String="keep", STD_AUTO_PROBE::String="hutch")

    ev_all = sort(unique(vcat(gi, gj)))
    nNodes = length(ev_all)
    map = Dict(eid => i for (i,eid) in enumerate(ev_all))
    li = [map[i] for i in gi]; lj = [map[j] for j in gj]
    theta_local = fill(NaN, nNodes); ref_local = fill(NaN, nNodes)
    std_local = fill(NaN, nNodes); nodew_local = fill(NaN, nNodes)
    deg_local = zeros(Int, nNodes)
    keep_mask = isfinite.(w_base) .& (w_base .> 0.0) .& isfinite.(gd)
    wrob_full = zeros(Float64, length(gd))
    diagnostics = NamedTuple[]
    component_sigmas = Float64[]
    n_pseudo_comps = 0

    # First pass: initial component gate, then residual-based pruning.
    comps = uf_components(nNodes, li[keep_mask], lj[keep_mask])
    for c in unique(comps)
        nodes = findall(==(c), comps)
        edges = findall(keep_mask .& (comps[li] .== c) .& (comps[lj] .== c))
        if length(edges) < MIN_EDGES || length(edges) < length(nodes)
            keep_mask[edges] .= false
            continue
        end
        localmap = Dict(node => i for (i,node) in enumerate(nodes))
        a = [localmap[li[e]] for e in edges]; b = [localmap[lj[e]] for e in edges]
        degree = nodeDegrees(length(nodes), a, b)
        pin = argmax(degree)
        free = [i for i in eachindex(nodes) if i != pin]
        theta = solvePinnedLS(length(nodes), a, b, gd[edges], w_base[edges], pin, free, RIDGE_EPS)
        r = theta[a] .- theta[b] .- gd[edges]
        scale = 1.4826 * mad1(r)
        if scale == 0.0 || !isfinite(scale)
            scale = max(ROBUST_MIN_SCALE, std(r))
        end
        keep_mask[edges] .= abs.(r) .<= K_SIGMA * max(scale, ROBUST_MIN_SCALE)
    end

    # References, fits and error estimates must describe the surviving graph.
    comps = uf_components(nNodes, li[keep_mask], lj[keep_mask])
    for c in unique(comps)
        nodes_c = findall(==(c), comps)
        idx_e = findall(keep_mask .& (comps[li] .== c) .& (comps[lj] .== c))
        n_c = length(nodes_c)
        # Repeat the microfamily gate after pruning; isolated rows stay NaN.
        if length(idx_e) < MIN_EDGES || length(idx_e) < n_c
            keep_mask[idx_e] .= false
            continue
        end
        mapc = Dict(node => i for (i,node) in enumerate(nodes_c))
        a = [mapc[li[e]] for e in idx_e]; b = [mapc[lj[e]] for e in idx_e]
        dt_k = gd[idx_e]; wb_k = w_base[idx_e]
        deg_kept = nodeDegrees(n_c, a, b)
        pin = argmax(deg_kept) # highest surviving degree; ties use smallest event ID
        keep_idx = [i for i in 1:n_c if i != pin]
        w_eff = copy(wb_k); wrob_c = ones(length(idx_e))
        theta_c = solvePinnedLS(n_c, a, b, dt_k, w_eff, pin, keep_idx, RIDGE_EPS)
        Jprev = Inf
        for it in 1:IRLS_ITERS
            r = theta_c[a] .- theta_c[b] .- dt_k
            scale = max(residualScale(r, w_eff), ROBUST_MIN_SCALE)
            wrob_c = huberWeight(r ./ scale, C_HUBER)
            w_eff = wb_k .* wrob_c
            # Always solve with the weights that will be exported/used for covariance.
            theta_c = solvePinnedLS(n_c, a, b, dt_k, w_eff, pin, keep_idx, RIDGE_EPS)
            r = theta_c[a] .- theta_c[b] .- dt_k
            J = sum(w_eff .* r.^2)
            if isfinite(Jprev) && abs(Jprev - J) <= IRLS_REL_TOL * max(Jprev, eps())
                break
            end
            Jprev = J
        end
        r = theta_c[a] .- theta_c[b] .- dt_k
        sigma_hat = residualScale(r, w_eff)
        push!(component_sigmas, sigma_hat)
        deg_local[nodes_c] .= deg_kept

        nodew_c = fill(NaN, n_c)
        if THETASTD_MODE == "pseudo_weight" || THETASTD_FALLBACK == "pseudo_weight" || WRITE_ALT_NODEW_COL
            src = PSEUDO_WEIGHT_SOURCE == "base" ? wb_k : PSEUDO_WEIGHT_SOURCE == "robust" ? wrob_c : w_eff
            sums = zeros(n_c)
            for k in eachindex(a)
                sums[a[k]] += src[k]; sums[b[k]] += src[k]
            end
            nodew_c = sums ./ deg_kept
            nodew_local[nodes_c] .= nodew_c
        end

        opts=(mode=EXPORT_THETA_STD ? THETASTD_MODE : "none",auto_probe=STD_AUTO_PROBE,fallback=THETASTD_FALLBACK,
            probes=STD_PROBES,probe_dist=STD_PROBE_DIST,batch=STD_BATCH,report_every_batch=STD_REPORT_EVERY,
            max_nred=STD_MAX_NRED,min_probes=STD_MIN_PROBES,rel_tol=STD_REL_TOL,
            exact_max_nred=STD_EXACT_MAX_NRED,repair_max_nodes=STD_REPAIR_MAX_NODES,unconverged=STD_UNCONVERGED)
        weight_proxy=1.0 ./ (max.(nodew_c,PSEUDO_WEIGHT_EPS).*THETASTD_SCALE_FIXED)
        std_c,info=estimateThetaStd(n_c,a,b,w_eff,pin,keep_idx,RIDGE_EPS,
            sigmaForThetaStd(sigma_hat,STD_MIN_SIGMA,STD_APPLY_MIN_SIGMA),opts;pseudo_weight=weight_proxy)
        EXPORT_THETA_STD || (info=merge(info,(method="disabled",)))
        any(m->m in ("pseudo_degree","pseudo_weight"),info.node_methods) && (n_pseudo_comps+=1)
        ref=ev_all[nodes_c[pin]]
        push!(diagnostics,merge(info,(ref=ref,nodes=n_c,event_ids=ev_all[nodes_c],sigma_hat=sigma_hat,
            sigma_std=sigmaForThetaStd(sigma_hat,STD_MIN_SIGMA,STD_APPLY_MIN_SIGMA))))
        theta_local[nodes_c] .= theta_c
        ref_local[nodes_c] .= Float64(ref)
        std_local[nodes_c] .= std_c
        wrob_full[idx_e] .= wrob_c
    end
    finite_sigmas = filter(isfinite, component_sigmas)
    sigma_hat_global = isempty(finite_sigmas) ? NaN : median(finite_sigmas)
    pseudo_info = (used_pseudo=n_pseudo_comps > 0, n_pseudo_comps=n_pseudo_comps, components=diagnostics)
    return ev_all, theta_local, ref_local, std_local, deg_local, nodew_local,
           keep_mask, wrob_full, sigma_hat_global, pseudo_info
end

function nodeDegrees(n, a, b)
    degree = zeros(Int, n)
    for k in eachindex(a)
        degree[a[k]] += 1; degree[b[k]] += 1
    end
    return degree
end

function residualScale(r, w)
    pos = isfinite.(w) .& (w .> 0)
    any(pos) || return NaN
    use = pos .& (w .>= 0.5 * median(w[pos]))
    s = 1.4826 * mad1(r[use])
    if !isfinite(s) || s == 0.0
        s = std(r[pos])
    end
    return (isfinite(s) && s >= 0.0) ? s : NaN
end

function solvePinnedLS(n_c::Int, a::Vector{Int}, b::Vector{Int},
    dt_k::Vector{Float64}, w::Vector{Float64},
    pin::Int, keep_idx::Vector{Int}, RIDGE_EPS::Float64)

    theta = zeros(Float64, n_c)

    valid = w .> 0.0
    if !any(valid)
        theta[pin] = 0.0
        return theta
    end

    aa = a[valid]; bb = b[valid]
    ww = w[valid]; dd = dt_k[valid]

    # Build weighted Laplacian
    Wadj = sparse(vcat(aa, bb), vcat(bb, aa), vcat(ww, ww), n_c, n_c)
    dsum = vec(sum(Wadj, dims=2))
    L = spdiagm(0 => dsum) - Wadj

    # bvec = A' W d (incidence form)
    bvec = zeros(Float64, n_c)
    @inbounds for k in eachindex(aa)
        i = aa[k]; j = bb[k]
        wd = ww[k] * dd[k]
        bvec[i] += wd
        bvec[j] -= wd
    end

    if isempty(keep_idx)
        theta[pin] = 0.0
        return theta
    end

    Lred = L[keep_idx, keep_idx] + (RIDGE_EPS * speye(length(keep_idx)))

    # Prefer Cholesky (SPD) if possible
    x = nothing
    try
        F = cholesky(Symmetric(Lred))
        x = F \ bvec[keep_idx]
    catch
        x = Lred \ bvec[keep_idx]
    end

    theta[pin] = 0.0
    theta[keep_idx] .= x
    return theta
end

speye(n::Int) = spdiagm(0 => ones(Float64, n))

function sigmaForThetaStd(sigma_hat::Float64, minSigma::Float64, applyMinSigma::Bool)
    sigma_std = (isfinite(sigma_hat) && sigma_hat >= 0.0) ? sigma_hat : NaN
    if applyMinSigma && isfinite(minSigma) && minSigma > 0.0
        sigma_std = max(sigma_std, minSigma)
    end
    return sigma_std
end

end # module
