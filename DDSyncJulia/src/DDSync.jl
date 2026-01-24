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
#     to avoid maxEventID × nGroups dense storage in memory.

using SparseArrays
using LinearAlgebra
using Statistics
using Random
using Printf

# ---------------------------- Public API ----------------------------

"""
    config_default() -> Dict{Symbol,Any}

Return a nested configuration dictionary whose fields mirror the MATLAB standalone
package `+ddsync/config_default.m` (as shipped in DDSync_standalone_fixed.zip).
"""
function config_default()
    cfg = Dict{Symbol,Any}()

    cfg[:io] = Dict{Symbol,Any}(
        :infile_dt     => "dt.cc",
        :catalog_file  => "catalog.txt",
        :tmpdir        => "synchro_tmp",
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
        :weight_mode          => "thetaStd", # "base" | "robust" | "combined" | "thetaStd"
        :thetastd_scale_mode  => "fixed",    # "fixed" | "median"
        :thetastd_scale_fixed => 500.0,
        :thetastd_weight_cap  => 1.0,
        :dt_decimals          => 5,
        :dt_weight_decimals   => 4,
        :station_field_width  => 8,
        :dt_field_width       => 10,
        :weight_field_width   => 8,
        :write_pruned_edges   => false,
    )

    cfg[:std] = Dict{Symbol,Any}(
        :export                           => true,
        :mode                             => "hutch",         # "hutch" | "pseudo_degree" | "pseudo_weight" | "none"
        :fallback                         => "pseudo_degree", # fallback if "hutch" is skipped
        :probes                           => 100000,
        :probe_dist                       => "rademacher",    # "rademacher" | "gaussian"
        :batch                            => 250,
        :report_every_batch               => 10,
        :max_nred                         => 40000,
        :min_sigma                        => 5e-4,
        :min_diag_rel                     => 1e-12,
        :pseudo_weight_source             => "combined",      # "base" | "robust" | "combined"
        :pseudo_weight_eps                => 0.05,
        :thetastd_write_weightcol         => true,
        :thetastd_write_alt_node_weightcol => false,
        :thetastd_fallback_weight_scale   => 1000.0,          # shared with output.thetastd_scale_fixed by default
    )

    cfg[:runtime] = Dict{Symbol,Any}(
        :prog_every_lines        => 2_000_000,
        :print_new_groups        => false,
        # If true, store per-group theta/std as dense Float32 vectors length maxEventID (fast, memory heavy).
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

    @printf("Reading %s to get max EventID...\n", catalog_file)
    maxEventID = readMaxEventID(catalog_file)
    @printf("Max EventID = %d\n", maxEventID)

    @printf("Pass 1: streaming %s and spooling groups...\n", infile_dt)
    groupKeys, groupFiles, groupCounts = spoolGroups(infile_dt, tmpdir;
        prog_every = Int(rt_cfg[:prog_every_lines]),
        print_new_groups = Bool(rt_cfg[:print_new_groups])
    )
    nGroups = length(groupKeys)
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

    # IRLS
    IRLS_ITERS   = Int(irls_cfg[:iters])
    C_HUBER      = Float64(irls_cfg[:C_HUBER])
    IRLS_REL_TOL = Float64(irls_cfg[:rel_tol])

    RIDGE_EPS = Float64(num_cfg[:RIDGE_EPS])

    # Output
    OUTPUT_WEIGHT_MODE   = lowercase(String(out_cfg[:weight_mode]))
    THETASTD_SCALE_MODE  = lowercase(String(out_cfg[:thetastd_scale_mode]))
    THETASTD_SCALE_FIXED = Float64(out_cfg[:thetastd_scale_fixed])
    THETASTD_WEIGHT_CAP  = Float64(out_cfg[:thetastd_weight_cap])
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
    STD_MIN_DIAGREL   = Float64(std_cfg[:min_diag_rel])
    PSEUDO_WEIGHT_SOURCE = lowercase(String(std_cfg[:pseudo_weight_source]))
    PSEUDO_WEIGHT_EPS    = Float64(std_cfg[:pseudo_weight_eps])
    WRITE_THETASTD_WEIGHTCOL = Bool(std_cfg[:thetastd_write_weightcol])
    WRITE_THETASTD_ALT_NODEW_COL = Bool(std_cfg[:thetastd_write_alt_node_weightcol])
    THETASTD_FALLBACK_WEIGHT_SCALE = Float64(std_cfg[:thetastd_fallback_weight_scale])

    # Runtime
    store_dense = Bool(rt_cfg[:store_dense_group_arrays])
    decision_chunk = Int(rt_cfg[:decision_chunk_records])
    median_cap = Int(rt_cfg[:median_scale_cap])

    if EXPORT_THETA_STD && (THETASTD_MODE == "hutch") && STD_PROBES > 0
        @printf("\n[STD WARNING]\n")
        @printf("  Hutchinson diag(inv(L)) convergence is slow: estimator error ~ O(1/sqrt(K)).\n")
        @printf("  If you want ms-level stability in std(theta), K may need to be VERY large.\n\n")
    end

    # Metrics file
    fidm = open(metrics_file, "w")
    @printf(fidm, "# DDSync metrics
")
    mode_print = String(out_cfg[:weight_mode])
    dist_print = String(std_cfg[:probe_dist])
    thetastd_mode_print = String(std_cfg[:mode])
    @printf(fidm, "# K_SIGMA=%g MIN_EDGES=%d IRLS_ITERS=%d C_HUBER=%g MODE=%s EXPORT_THETA_STD=%d STD_PROBES=%d DIST=%s THETASTD_MODE=%s
",
        K_SIGMA, MIN_EDGES, IRLS_ITERS, C_HUBER, mode_print,
        EXPORT_THETA_STD ? 1 : 0, STD_PROBES, dist_print, thetastd_mode_print)
    @printf(fidm, "group\tn_edges\tn_kept\tpruned_pct\trobust_p50\trobust_p05\trobust_p95\tfrac_robust_lt_0p99\tsigma_hat\n")
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
            K_SIGMA=K_SIGMA, MIN_EDGES=MIN_EDGES,
            IRLS_ITERS=IRLS_ITERS, C_HUBER=C_HUBER, IRLS_REL_TOL=IRLS_REL_TOL,
            RIDGE_EPS=RIDGE_EPS,
            EXPORT_THETA_STD=EXPORT_THETA_STD,
            THETASTD_MODE=THETASTD_MODE, THETASTD_FALLBACK=THETASTD_FALLBACK,
            STD_PROBES=STD_PROBES, STD_PROBE_DIST=STD_PROBE_DIST,
            STD_BATCH=STD_BATCH, STD_REPORT_EVERY=STD_REPORT_EVERY,
            STD_MAX_NRED=STD_MAX_NRED, STD_MIN_SIGMA=STD_MIN_SIGMA, STD_MIN_DIAGREL=STD_MIN_DIAGREL,
            PSEUDO_WEIGHT_SOURCE=PSEUDO_WEIGHT_SOURCE, PSEUDO_WEIGHT_EPS=PSEUDO_WEIGHT_EPS,
            THETASTD_SCALE_FIXED=THETASTD_SCALE_FIXED,
            WRITE_ALT_NODEW_COL=WRITE_THETASTD_ALT_NODEW_COL,
        )

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

        # Write theta + std files (full 1..maxEventID rows, MATLAB-compatible)
        sta, ph = splitKey(key)

        theta_fn = joinpath(thetadir, @sprintf("theta_%s_%s.txt", sta, ph))
        writeThetaTriples(theta_fn, maxEventID, ev_all, theta_local, ref_local)

        if EXPORT_THETA_STD
            std_fn = joinpath(thetastd_dir, @sprintf("std_theta_%s_%s.txt", sta, ph))
            writeThetaStdWithDegree(std_fn, maxEventID,
                ev_all, std_local, ref_local, deg_local, nodew_local;
                write_weightcol=WRITE_THETASTD_WEIGHTCOL,
                write_alt_nodew_col=WRITE_THETASTD_ALT_NODEW_COL,
                thetastd_scale_fixed=THETASTD_SCALE_FIXED,
                thetastd_weight_cap=THETASTD_WEIGHT_CAP
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

        @printf(fidm, "%s\t%d\t%d\t%.3f\t%.4f\t%.4f\t%.4f\t%.3f\t%.6g\n",
            key, length(gd), n_kept, pruned_pct, r50, r05, r95, frac_dw, sigma_hat)
        flush(fidm)

        @printf("  kept=%d (%.2f%% pruned), sigma_hat=%.4g s, robust_w p05/p50/p95=%.3g/%.3g/%.3g\n",
            n_kept, pruned_pct, sigma_hat, r05, r50, r95)
    end

    overall_pruned_pct = 100.0 * (overall_edges - overall_kept) / max(overall_edges, 1)

    @printf(fidm, "\n# OVERALL\n")
    @printf(fidm, "overall_edges=%d\n", overall_edges)
    @printf(fidm, "overall_kept=%d\n", overall_kept)
    @printf(fidm, "overall_pruned_pct=%.3f\n", overall_pruned_pct)
    @printf(fidm, "pseudo_groups_used=%d\n", nPseudoGroups)
    @printf(fidm, "pseudo_components_used=%d\n", nPseudoComps)
    close(fidm)

    # Optional median scaling for thetaStd weights
    thetaStdScale = THETASTD_SCALE_FIXED
    if OUTPUT_WEIGHT_MODE == "thetastd" && THETASTD_SCALE_MODE == "median"
        @printf("Computing global median scale for thetaStd weights (extra streaming pass)...\n")
        thetaStdScale = computeThetaStdMedianScale(infile_dt, key2g, decFiles, std_store, store_dense, maxEventID; cap=median_cap, chunk=decision_chunk)
        if !(isfinite(thetaStdScale) && thetaStdScale > 0)
            thetaStdScale = THETASTD_SCALE_FIXED
        end
        @printf("thetaStdScale (median of 1/std_dt) = %.6g\n", thetaStdScale)
    end

    # Pass 2: write dt_sync.cc
    @printf("Pass 2: writing %s ...\n", out_dt_sync)
    readers = openDecisionReaders(decFiles; chunk=decision_chunk)

    fin  = open(infile_dt, "r")
    fout = open(out_dt_sync, "w")

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

    for L in eachline(fin)
        Lt = strip(L)
        if isempty(Lt)
            continue
        end

        if startswith(Lt, "#")
            pendingHeader = L
            pendingPrinted = false

            nums = split(strip(Lt[2:end]))
            if length(nums) >= 2
                try
                    cur_i = parse(Int, nums[1])
                    cur_j = parse(Int, nums[2])
                catch
                    cur_i = 0; cur_j = 0
                end
            else
                cur_i = 0; cur_j = 0
            end
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
                pendingPrinted = true
            end
            dt_corr = dt_in
            w_out = max(0.0, isfinite(cc_in) ? weight_fun_base(cc_in) : 0.0)

            @printf(fout, "    %-*s %*.*f %*.*f %s\n",
                STA_W, sta, DT_W, DT_DEC, dt_corr, W_W, W_DEC, w_out, ph)
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
                std_dt = sqrt(max(0.0, si^2 + sj^2))
                std_dt = max(std_dt, 1e-12)
                w_raw = 1.0 / std_dt
                w_out = w_raw / max(thetaStdScale, eps())
            else
                # default
                w_out = max(0.0, isfinite(cc_in) ? weight_fun_base(cc_in) : 0.0) * wrob
            end
        else
            # pruned edge written: weight 0
            w_out = 0.0
        end

        if isfinite(THETASTD_WEIGHT_CAP)
            w_out = min(w_out, THETASTD_WEIGHT_CAP)
        end

        if !pendingPrinted && !isempty(pendingHeader)
            println(fout, pendingHeader)
            pendingPrinted = true
        end

        @printf(fout, "    %-*s %*.*f %*.*f %s\n",
            STA_W, sta, DT_W, DT_DEC, dt_corr, W_W, W_DEC, w_out, ph)
    end

    close(fin); close(fout)
    closeDecisionReaders(readers)

    @printf("Done.\n  Synced dt: %s\n  Metrics:  %s\n", out_dt_sync, metrics_file)
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

# Read max EventID (last token per nonempty line)
function readMaxEventID(catalog_file::String)
    maxID = 0
    open(catalog_file, "r") do io
        for L in eachline(io)
            Lt = strip(L)
            isempty(Lt) && continue
            toks = split(Lt)
            id = tryparse(Int, toks[end])
            if id !== nothing
                maxID = max(maxID, id)
            end
        end
    end
    return maxID
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
function spoolGroups(infile_dt::String, tmpdir::String; prog_every::Int=2_000_000, print_new_groups::Bool=false)
    keys = String[]
    files = String[]
    counts = Int[]

    key2idx = Dict{String,Int}()
    fids = IOStream[]

    cur_i = 0
    cur_j = 0
    nLine = 0

    open(infile_dt, "r") do io
        for L in eachline(io)
            nLine += 1
            if prog_every > 0 && (nLine % prog_every == 0)
                @printf("  ... streamed %d lines\n", nLine)
            end
            Lt = strip(L)
            isempty(Lt) && continue

            if startswith(Lt, "#")
                nums = split(strip(Lt[2:end]))
                if length(nums) >= 2
                    i = tryparse(Int, nums[1])
                    j = tryparse(Int, nums[2])
                    if i !== nothing && j !== nothing
                        cur_i = i
                        cur_j = j
                    end
                end
                continue
            end

            toks = split(Lt)
            length(toks) < 4 && continue

            sta = toks[1]
            dt  = tryparse(Float64, toks[2])
            cc  = tryparse(Float64, toks[3])
            ph  = toks[4]
            (dt === nothing || cc === nothing) && continue

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

    for fid in fids
        close(fid)
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
    chunk::Int=200_000)

    readers = openDecisionReaders(decFiles; chunk=chunk)
    vals = Float64[]
    sizehint!(vals, min(cap, 200_000))

    cur_i = 0
    cur_j = 0

    open(infile_dt, "r") do io
        for L in eachline(io)
            Lt = strip(L)
            isempty(Lt) && continue

            if startswith(Lt, "#")
                nums = split(strip(Lt[2:end]))
                if length(nums) >= 2
                    i = tryparse(Int, nums[1])
                    j = tryparse(Int, nums[2])
                    if i !== nothing && j !== nothing
                        cur_i = i
                        cur_j = j
                    end
                end
                continue
            end

            toks = split(Lt)
            length(toks) < 4 && continue
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

            std_dt = sqrt(max(0.0, si^2 + sj^2))
            std_dt = max(std_dt, 1e-12)
            v = 1.0 / std_dt
            if isfinite(v) && v > 0.0
                push!(vals, v)
            end

            if length(vals) >= cap
                break
            end
        end
    end

    closeDecisionReaders(readers)

    if isempty(vals)
        return 1.0
    end
    return median(vals)
end

# -------- Output writers --------

function writeThetaTriples(fn::String, maxEventID::Int, ev_all::Vector{Int}, theta_local::Vector{Float64}, ref_local::Vector{Float64})
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
            @printf(io, "%d %.6f %.0f\n", ev, th, rf)
        end
    end
end

function writeThetaStdWithDegree(fn::String, maxEventID::Int,
    ev_all::Vector{Int}, std_local::Vector{Float64}, ref_local::Vector{Float64},
    deg_local::Vector{Int}, nodew_local::Vector{Float64};
    write_weightcol::Bool=true,
    write_alt_nodew_col::Bool=false,
    thetastd_scale_fixed::Float64=1000.0,
    thetastd_weight_cap::Float64=1.0)

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

            # Optional weight column derived from std(theta)
            if write_weightcol
                if !(isfinite(s) && s > 0.0)
                    wtheta = 0.0
                else
                    wtheta = (1.0 / s) / max(thetastd_scale_fixed, eps())
                    if isfinite(thetastd_weight_cap)
                        wtheta = min(wtheta, thetastd_weight_cap)
                    end
                end
            end

            if write_weightcol && write_alt_nodew_col
                @printf(io, "%d\t%.6f\t%.0f\t%d\t%.6f\t%.6f\n", ev, s, rf, dg, wtheta, nw)
            elseif write_weightcol
                @printf(io, "%d\t%.6f\t%.0f\t%d\t%.6f\n", ev, s, rf, dg, wtheta)
            elseif write_alt_nodew_col
                @printf(io, "%d\t%.6f\t%.0f\t%d\t%.6f\n", ev, s, rf, dg, nw)
            else
                @printf(io, "%d\t%.6f\t%.0f\t%d\n", ev, s, rf, dg)
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
    STD_MIN_DIAGREL::Float64,
    PSEUDO_WEIGHT_SOURCE::String,
    PSEUDO_WEIGHT_EPS::Float64,
    THETASTD_SCALE_FIXED::Float64,
    WRITE_ALT_NODEW_COL::Bool=false)

    # Unique events in this group (MATLAB unique() => sorted ascending)
    ev_all = sort(unique(vcat(gi, gj)))
    nNodes = length(ev_all)

    # Map global EventID -> local 1..nNodes
    map = Dict{Int,Int}()
    for (i,eid) in enumerate(ev_all)
        map[eid] = i
    end
    li = [map[i] for i in gi]
    lj = [map[j] for j in gj]

    theta_local = fill(NaN, nNodes)
    ref_local   = fill(NaN, nNodes)   # store GLOBAL ref EventID (float, printed with %.0f)
    std_local   = fill(NaN, nNodes)
    nodew_local = fill(NaN, nNodes)
    deg_local   = zeros(Int, nNodes)

    keep_mask = (w_base .> 0.0)
    wrob_full = ones(Float64, length(gd))

    valid_edges = findall(>(0.0), w_base)
    if isempty(valid_edges)
        # No usable edges: return all-NaN theta, keep_mask already false
        pseudo_info = (used_pseudo=false, n_pseudo_comps=0)
        return ev_all, theta_local, ref_local, std_local, deg_local, nodew_local, keep_mask, wrob_full, NaN, pseudo_info
    end

    comps = uf_components(nNodes, li[valid_edges], lj[valid_edges])
    nComp = maximum(comps)

    comp_sigma = fill(NaN, nComp)
    comp_w     = zeros(Int, nComp)

    used_pseudo_any = false
    n_pseudo_comps = 0

    for c in 1:nComp
        nodes_c = findall(==(c), comps)
        n_c = length(nodes_c)

        # Edges within this component AND positive base weight
        mE = BitVector(undef, length(gd))
        @inbounds for k in eachindex(gd)
            mE[k] = (w_base[k] > 0.0) && (comps[li[k]] == c) && (comps[lj[k]] == c)
        end

        if !any(mE)
            # isolated nodes / no valid edges
            @inbounds for idx in nodes_c
                theta_local[idx] = 0.0
                ref_local[idx] = Float64(ev_all[nodes_c[1]])
                if EXPORT_THETA_STD
                    std_local[idx] = 0.0
                end
            end
            continue
        end

        idx_e = findall(mE)
        nE = length(idx_e)
        redundancy = nE - (n_c - 1)

        # Microfamily gate: require enough edges and at least one cycle
        if nE < MIN_EDGES || redundancy < 1
            keep_mask[idx_e] .= false
            wrob_full[idx_e] .= 0.0
            continue
        end

        # Extract component edges
        a_all = li[idx_e]
        b_all = lj[idx_e]
        dt_k  = gd[idx_e]
        wb_k  = w_base[idx_e]

        # Map component node indices -> 1..n_c
        mapc = Dict{Int,Int}()
        for (ii,node) in enumerate(nodes_c)
            mapc[node] = ii
        end
        a = [mapc[x] for x in a_all]
        b = [mapc[x] for x in b_all]

        # Degree for pin selection (unweighted)
        deg0 = zeros(Int, n_c)
        @inbounds for k in eachindex(a)
            deg0[a[k]] += 1
            deg0[b[k]] += 1
        end
        pin = argmax(deg0)  # first max => stable tie-break
        keep_idx = [i for i in 1:n_c if i != pin]

        # LS (all edges)
        theta_c = solvePinnedLS(n_c, a, b, dt_k, wb_k, pin, keep_idx, RIDGE_EPS)

        # residuals
        r = (theta_c[a] .- theta_c[b]) .- dt_k

        # Prune by MAD
        if length(r) >= MIN_EDGES
            s0 = 1.4826 * mad1(r)
            if s0 == 0.0 || !isfinite(s0)
                s0 = max(STD_MIN_SIGMA, std(r))
            end
            bad = abs.(r) .> (K_SIGMA * s0)
        else
            bad = falses(length(r))
        end
        keep_c = .!bad

        # LS (pruned)
        w_eff = wb_k .* Float64.(keep_c)
        theta_c = solvePinnedLS(n_c, a, b, dt_k, w_eff, pin, keep_idx, RIDGE_EPS)

        # IRLS (Huber)
        wrob_c = ones(Float64, length(dt_k))
        if IRLS_ITERS > 0
            Jprev = Inf
            for it in 1:IRLS_ITERS
                r = (theta_c[a] .- theta_c[b]) .- dt_k
                ww = w_eff
                pos = ww .> 0.0
                if any(pos)
                    medw = median(ww[pos])
                    use = pos .& (ww .>= 0.5 * medw)
                    if any(use)
                        s = 1.4826 * mad1(r[use])
                    else
                        s = 1.4826 * mad1(r[pos])
                    end
                    if s == 0.0 || !isfinite(s)
                        s = max(STD_MIN_SIGMA, std(r[pos]))
                    end
                else
                    s = STD_MIN_SIGMA
                end
                rs = r ./ s

                wrob_c = huberWeight(rs, C_HUBER)
                w_eff  = wb_k .* Float64.(keep_c) .* wrob_c

                J = sum(w_eff .* (r .^ 2))
                if isfinite(Jprev) && (Jprev - J) <= IRLS_REL_TOL * max(Jprev, eps())
                    break
                end
                Jprev = J

                theta_c = solvePinnedLS(n_c, a, b, dt_k, w_eff, pin, keep_idx, RIDGE_EPS)
            end
        end

        # Final sigma_hat (robust scale, high-weight subset)
        r = (theta_c[a] .- theta_c[b]) .- dt_k
        ww = w_eff
        pos = ww .> 0.0
        if any(pos)
            medw = median(ww[pos])
            use = pos .& (ww .>= 0.5 * medw)
            if any(use)
                sigma_hat = 1.4826 * mad1(r[use])
            else
                sigma_hat = 1.4826 * mad1(r[pos])
            end
            if sigma_hat == 0.0 || !isfinite(sigma_hat)
                sigma_hat = max(STD_MIN_SIGMA, std(r[pos]))
            end
        else
            sigma_hat = STD_MIN_SIGMA
        end
        sigma_hat = max(sigma_hat, STD_MIN_SIGMA)
        comp_sigma[c] = sigma_hat

        kept = w_eff .> 0.0
        comp_w[c] = count(identity, kept)

        # Degree on kept graph (unweighted)
        deg_kept = zeros(Int, n_c)
        @inbounds for k in eachindex(a)
            if kept[k]
                deg_kept[a[k]] += 1
                deg_kept[b[k]] += 1
            end
        end
        # Write degree back to local arrays (global node indices)
        @inbounds for (ii,node) in enumerate(nodes_c)
            deg_local[node] = deg_kept[ii]
        end

        # Node weights for pseudo_weight and/or alt output column
        need_nodew = (THETASTD_MODE == "pseudo_weight") ||
                     ((THETASTD_MODE == "hutch") && (THETASTD_FALLBACK == "pseudo_weight")) ||
                     WRITE_ALT_NODEW_COL
        nodew_c = fill(NaN, n_c)
        if need_nodew
            kept_for_nodew = kept
            if any(kept_for_nodew)
                ww_src = similar(w_eff)
                if PSEUDO_WEIGHT_SOURCE == "base"
                    ww_src .= wb_k
                elseif PSEUDO_WEIGHT_SOURCE == "robust"
                    ww_src .= wrob_c
                else
                    ww_src .= w_eff
                end
                sumw = zeros(Float64, n_c)
                cntw = zeros(Int, n_c)
                @inbounds for k in eachindex(a)
                    if kept_for_nodew[k]
                        wv = ww_src[k]
                        ai = a[k]; bi = b[k]
                        sumw[ai] += wv; sumw[bi] += wv
                        cntw[ai] += 1;  cntw[bi] += 1
                    end
                end
                @inbounds for i in 1:n_c
                    nodew_c[i] = sumw[i] / max(cntw[i], 1)
                end
                @inbounds for (ii,node) in enumerate(nodes_c)
                    nodew_local[node] = nodew_c[ii]
                end
            end
        end

        # Theta std
        std_c = fill(NaN, n_c)
        std_c[pin] = 0.0

        run_hutch = EXPORT_THETA_STD &&
                    (THETASTD_MODE == "hutch") &&
                    (STD_PROBES > 0) &&
                    (length(keep_idx) <= STD_MAX_NRED)

        if EXPORT_THETA_STD
            if run_hutch
                std_c = estimateThetaStdHutch(n_c, a, b, w_eff, pin, keep_idx, RIDGE_EPS,
                    sigma_hat, STD_PROBES, STD_PROBE_DIST, STD_BATCH, STD_REPORT_EVERY, STD_MIN_DIAGREL)
            else
                use_mode = THETASTD_MODE
                if use_mode == "hutch"
                    use_mode = THETASTD_FALLBACK
                end
                if use_mode == "pseudo_degree"
                    degv = max.(Float64.(deg_kept), 1.0)
                    std_c = (sigma_hat ./ sqrt.(degv))
                    std_c[pin] = 0.0
                    used_pseudo_any = true
                    n_pseudo_comps += 1
                elseif use_mode == "pseudo_weight"
                    wnode = max.(nodew_c, PSEUDO_WEIGHT_EPS)
                    std_c = 1.0 ./ (wnode .* THETASTD_SCALE_FIXED)
                    std_c[pin] = 0.0
                    used_pseudo_any = true
                    n_pseudo_comps += 1
                elseif use_mode == "none"
                    std_c .= NaN
                    std_c[pin] = 0.0
                else
                    std_c .= NaN
                    std_c[pin] = 0.0
                end
            end
        end

        # Assign theta/ref/std to global local arrays for nodes in this component
        ref_global = Float64(ev_all[nodes_c[pin]])
        @inbounds for (ii,node) in enumerate(nodes_c)
            theta_local[node] = theta_c[ii]
            ref_local[node]   = ref_global
            if EXPORT_THETA_STD
                std_local[node] = std_c[ii]
            end
        end

        # Update per-edge keep + wrob in original edge order
        keep_mask[idx_e] .= keep_c
        @inbounds for (kk, eidx) in enumerate(idx_e)
            if keep_c[kk]
                wrob_full[eidx] = wrob_c[kk]
            else
                wrob_full[eidx] = 0.0
            end
        end
    end

    # Representative sigma_hat per group (reporting only): median over comps with kept edges
    good = [(comp_w[c] > 0) && isfinite(comp_sigma[c]) for c in 1:nComp]
    sigma_hat_global = NaN
    if any(good)
        sigma_hat_global = median(comp_sigma[good])
    end

    pseudo_info = (used_pseudo=used_pseudo_any, n_pseudo_comps=n_pseudo_comps)

    return ev_all, theta_local, ref_local, std_local, deg_local, nodew_local,
           keep_mask, wrob_full, sigma_hat_global, pseudo_info
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

function estimateThetaStdHutch(n_c::Int, a::Vector{Int}, b::Vector{Int}, w_eff::Vector{Float64},
    pin::Int, keep_idx::Vector{Int}, RIDGE_EPS::Float64,
    sigma_hat::Float64, K::Int, dist::String, batch::Int, reportEvery::Int, minDiagRel::Float64)

    valid = w_eff .> 0.0
    aa = a[valid]; bb = b[valid]; ww = w_eff[valid]

    Wadj = sparse(vcat(aa, bb), vcat(bb, aa), vcat(ww, ww), n_c, n_c)
    dsum = vec(sum(Wadj, dims=2))
    L = spdiagm(0 => dsum) - Wadj

    nred = length(keep_idx)
    Lred = L[keep_idx, keep_idx] + (RIDGE_EPS * speye(nred))

    # Factorization
    F = nothing
    use_chol = true
    try
        F = cholesky(Symmetric(Lred))
    catch
        use_chol = false
        F = lu(Lred)
    end

    diag_est = zeros(Float64, nred)
    nDone = 0

    nBatches = Int(ceil(K / batch))
    for bi in 1:nBatches
        thisK = min(batch, K - nDone)
        thisK <= 0 && break

        Z = zeros(Float64, nred, thisK)
        if dist == "rademacher"
            Z .= sign.(randn(nred, thisK))
            Z[Z .== 0.0] .= 1.0
        elseif dist == "gaussian"
            Z .= randn(nred, thisK)
        else
            error("STD_PROBE_DIST must be rademacher or gaussian")
        end

        X = F \ Z
        diag_est .+= vec(sum(X .* Z, dims=2))

        nDone += thisK
        if reportEvery > 0 && (bi % reportEvery == 0)
            tmp = diag_est ./ max(nDone, 1)
            @printf("    Hutch diag: %d/%d probes (median diag ~ %.3e)\n", nDone, K, median(abs.(tmp)))
        end
    end

    diag_est ./= max(nDone, 1)

    medDiag = median(abs.(diag_est))
    floorDiag = max(medDiag * minDiagRel, 0.0)
    diag_est = max.(diag_est, floorDiag)

    std_red = sigma_hat .* sqrt.(max.(diag_est, 0.0))

    std_c = fill(NaN, n_c)
    std_c[pin] = 0.0
    std_c[keep_idx] .= std_red
    return std_c
end

end # module
