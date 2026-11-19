# Reporting counts exported non-reference station-phase/event entries, not
# components or unique catalog events. Exact includes selected exact repairs.
const STD_REVIEW_THRESHOLD_PCT = 10.0

function uncertaintyCounts(components)
    exact = hutch = cholesky = unconverged = cholesky_kept = fallback = disabled = 0
    for c in components
        exact += c.exact_nodes
        hutch += c.hutch_converged
        cholesky += get(c, :cholesky_converged, 0)
        cholesky_kept += get(c, :cholesky_unconverged, 0)
        unconverged += c.unconverged_nodes
        fallback += c.fallback_nodes
        c.method in ("none", "disabled") && (disabled += c.nodes - 1)
    end
    total = exact + hutch + cholesky + unconverged + fallback
    review = unconverged + fallback
    return (exact=exact, hutch=hutch, cholesky=cholesky, unconverged=unconverged, cholesky_kept=cholesky_kept, fallback=fallback,
        total=total, review=review, disabled=disabled, references=length(components))
end

# Index once: scanning all components for every station is costly on DAS data.
function uncertaintyByGroup(components)
    grouped = Dict{String,Vector{NamedTuple}}()
    for c in components
        push!(get!(() -> NamedTuple[], grouped, c.group), c)
    end
    return Dict(group=>uncertaintyCounts(entries) for (group,entries) in grouped)
end

stdPercent(n, total) = total == 0 ? "N/A" : @sprintf("%.2f%%",100.0*n/total)
stdCountPercent(n,total) = @sprintf("%d (%s)",n,stdPercent(n,total))
needsStdReview(c) = c.total > 0 && 100.0*c.review/c.total >= STD_REVIEW_THRESHOLD_PCT

function writeUncertaintySummary(io, components; prefix="# ")
    c = uncertaintyCounts(components)
    println(io,prefix,"UNCERTAINTY SUMMARY")
    println(io,prefix,"Denominator: ",c.total," supported non-reference std entries (station-phase/event rows).")
    println(io,prefix,"Excluded: ",c.references," gauge references; ",c.disabled," disabled entries; rejected/absent events.")
    @printf(io,"%s%-30s %12s %12s\n",prefix,"Category","Count","Percent")
    for (label,n) in (("Exact (including repairs)",c.exact),("Converged Hutch",c.hutch),
                      ("Kept unconverged Hutch",c.unconverged-c.cholesky_kept),
                      ("Converged Cholesky",c.cholesky),("Kept unconverged Cholesky",c.cholesky_kept),("Fallback / proxy",c.fallback))
        @printf(io,"%s%-30s %12d %12s\n",prefix,label,n,stdPercent(n,c.total))
    end
    println(io,prefix,"A draw is one random sample. Relative MC SE is the estimated sampling error divided by the event timing variance estimate.")
    println(io,prefix,"For example, 0.10 means 10% sampling error on the variance. The uncertainty exported in seconds is its square root after scaling.")
    println(io,prefix,"Requires review (kept unconverged + fallback/proxy): ",stdCountPercent(c.review,c.total))
    println(io,prefix,"Numerically checked (exact + converged estimates): ",stdCountPercent(c.exact+c.hutch+c.cholesky,c.total))
    println(io,prefix,"Fallback/proxy includes explicitly requested proxies/bounds. Checked does not mean physically calibrated.")
    if c.total == 0
        println(io,prefix,"No assessed std entries; percentages are N/A. Enable uncertainty if required by the application.")
    elseif needsStdReview(c)
        println(io,prefix,"ADVICE: At least 10% require review. Check the categories above and std_quality_*.txt.")
        println(io,prefix,"        If many Hutch estimates miss the target, try Cholesky before increasing the number of probes.")
        println(io,prefix,"        Under [std], set auto_probe = \"cholesky\" in auto mode, or change explicit hutch mode to mode = \"cholesky\".")
        println(io,prefix,"        For a few unresolved entries, consider raising std.repair_max_nodes. For fallbacks, check the component reason.")
    else
        println(io,prefix,"Below the 10% review advisory threshold. Inspect flagged entries and group summaries before weighting.")
    end
    println(io,prefix,"Sampled entries with relative MC SE >= 1, invalid estimates or too few probes cannot be kept;")
    println(io,prefix,"after selected repairs these entries use fallback. They are included in the fallback/proxy counts above.")
    return c
end

function writeMetricsReport(io, groups, components, overall)
    println(io)
    writeUncertaintySummary(io, components)
    println(io,"# Timing mismatch = fitted (theta_i - theta_j) minus measured differential travel time, for kept pairs.")
    println(io,"# sigma_hat: robust scatter of these mismatches per component (seconds); group entries report the median across components.")
    println(io,"# sigma_std = max(sigma_hat, min_sigma) when the minimum is enabled, otherwise sigma_hat; this is sigma in std formulas.")
    println(io,"# Neither scale is a theta std or a Monte Carlo error. pseudo_weight does not use sigma_std.")
    width = max(16, isempty(groups) ? 0 : maximum(length(g.group) for g in groups))
    println(io,"\n# STATION-PHASE OBSERVATIONS")
    @printf(io,"%-*s %10s %10s %11s %12s %12s %12s %20s %13s\n",width,"group",
        "n_edges","n_kept","pruned_pct","robust_p50","robust_p05","robust_p95","frac_robust_lt_0p99","sigma_hat")
    for g in groups
        @printf(io,"%-*s %10d %10d %11.3f %12.4f %12.4f %12.4f %20.3f %13.6g\n",width,g.group,
            g.n_edges,g.n_kept,g.pruned_pct,g.r50,g.r05,g.r95,g.frac_dw,g.sigma_hat)
    end
    println(io,"\n# OVERALL")
    @printf(io,"overall_edges=%d\noverall_kept=%d\noverall_pruned_pct=%.3f\npseudo_groups_used=%d\npseudo_components_used=%d\n",
        overall.edges,overall.kept,overall.pruned_pct,overall.pseudo_groups,overall.pseudo_components)
    by_group = uncertaintyByGroup(components)
    empty_counts = uncertaintyCounts(())
    println(io,"\n# UNCERTAINTY BY GROUP (counts and percentages use each group's assessed std entries)")
    @printf(io,"# %-*s %8s %18s %18s %18s %20s %20s %18s %12s %10s\n",width,"group","n_std",
        "exact","hutch_converged","hutch_unconverged","cholesky_converged","cholesky_unconverged","fallback/proxy","review_pct","disabled")
    for g in groups
        c = get(by_group,g.group,empty_counts)
        @printf(io,"# %-*s %8d %18s %18s %18s %20s %20s %18s %12s %10d\n",width,g.group,c.total,
            stdCountPercent(c.exact,c.total),stdCountPercent(c.hutch,c.total),stdCountPercent(c.unconverged-c.cholesky_kept,c.total),
            stdCountPercent(c.cholesky,c.total),stdCountPercent(c.cholesky_kept,c.total),
            stdCountPercent(c.fallback,c.total),stdPercent(c.review,c.total),c.disabled)
    end
    flagged = [g.group for g in groups if needsStdReview(get(by_group,g.group,empty_counts))]
    if !isempty(flagged)
        println(io,"# GROUP ADVICE: At least 10% require review in each group below, even if the overall share is small.")
        println(io,"# If many Hutch estimates miss the target, try Cholesky before increasing the number of probes.")
        println(io,"# Under [std], set auto_probe = \"cholesky\" in auto mode, or change explicit hutch mode to mode = \"cholesky\".")
        println(io,"# Check the component reason for fallback entries.")
        for group in flagged;println(io,"#   ",group);end
    end
    println(io,"\n# COMPONENT DETAILS (all uncertainty methods)")
    println(io,"# max_rel_se is the pre-repair diagnostic, not the error of exact repairs; scales are seconds.")
    for c in components
        @printf(io,"# std_component group=%s ref=%d nodes=%d method=%s\n",c.group,c.ref,c.nodes,c.method)
        @printf(io,"#   reason=%s probes=%d max_rel_se=%.6g\n",c.reason,c.probes,c.max_rel_se)
        @printf(io,"#   sigma_hat=%.6g sigma_std=%.6g repair_attempted=%d repaired_nodes=%d\n",c.sigma_hat,c.sigma_std,c.repair_attempted,c.repaired_nodes)
        @printf(io,"#   exact_nodes=%d hutch_converged=%d cholesky_converged=%d cholesky_unconverged=%d unconverged_nodes=%d fallback_nodes=%d\n\n",
            c.exact_nodes,c.hutch_converged,get(c,:cholesky_converged,0),get(c,:cholesky_unconverged,0),c.unconverged_nodes,c.fallback_nodes)
    end
end
