@testset "Readable uncertainty accounting" begin
    function component(group,ref; exact=0,hutch=0,kept=0,fallback=0,disabled=0,method="hybrid",reason="not_converged",repaired=0)
        (group=group,ref=ref,nodes=exact+hutch+kept+fallback+disabled+1,method=method,reason=reason,
         probes=4096,max_rel_se=1.2,exact_nodes=exact,hutch_converged=hutch,
         unconverged_nodes=kept,fallback_nodes=fallback,repair_attempted=repaired,repaired_nodes=repaired,
         sigma_hat=1e-5,sigma_std=5e-4)
    end
    comps=[component("GOOD_P",1;exact=800,hutch=190,repaired=40),
           component("BAD_S",2;kept=3,fallback=7),component("OFF_P",3;disabled=5,method="disabled")]
    c=DDSync.uncertaintyCounts(comps)
    @test c.total==1000 && c.exact==800 && c.review==10 && c.disabled==5 && c.references==3
    @test !DDSync.needsStdReview(c) # row-weighted 1%, not mean of 0% and 100%
    @test DDSync.needsStdReview(DDSync.uncertaintyCounts(comps[2:2]))
    @test DDSync.needsStdReview(DDSync.uncertaintyCounts([component("T",1;exact=90,kept=10)]))
    @test !DDSync.needsStdReview(DDSync.uncertaintyCounts([component("T",1;exact=91,fallback=9)]))
    summary=sprint(io->DDSync.writeUncertaintySummary(io,comps))
    for (label,n,pct) in (("Exact (including repairs)",800,"80.00%"),("Converged Hutch",190,"19.00%"),
                         ("Kept unconverged Hutch",3,"0.30%"),("Fallback / proxy",7,"0.70%"))
        line=only(filter(l->occursin(label,l),split(summary,'\n')))
        @test join(split(line)[end-1:end]," ")=="$n $pct"
    end
    groups=[(group=g,n_edges=100,n_kept=90,pruned_pct=10.0,r50=1.0,r05=0.8,r95=1.0,frac_dw=0.1,sigma_hat=1e-5)
        for g in ("GOOD_P","BAD_S","OFF_P","A_VERY_LONG_STATION_NAME_P")]
    overall=(edges=400,kept=360,pruned_pct=10.0,pseudo_groups=1,pseudo_components=1)
    report=sprint(io->DDSync.writeMetricsReport(io,groups,comps,overall))
    @test findfirst("# UNCERTAINTY SUMMARY",report)<findfirst("# STATION-PHASE OBSERVATIONS",report)<findfirst("# COMPONENT DETAILS",report)
    @test !occursin("# std_component",first(split(report,"# COMPONENT DETAILS")))
    @test occursin("# GROUP ADVICE",report) && occursin("#   BAD_S",report)
    @test !occursin("# ADVICE:",report) # overall does not hide per-group warning
    @test occursin("relative mc se >= 1",lowercase(report))
    @test occursin("sigma_hat=1e-05 sigma_std=0.0005",report)
    lines=split(report,'\n');header=only(filter(l->startswith(l,"group "),lines))
    # Right edges of all numeric columns coincide with their headers, including sigma_hat.
    spans=[m.offset+length(m.match)-1 for m in eachmatch(r"\S+",header)]
    for line in filter(l->startswith(l,"GOOD_P ")||startswith(l,"A_VERY_LONG_STATION_NAME_P "),lines)
        rowspans=[m.offset+length(m.match)-1 for m in eachmatch(r"\S+",line)]
        @test rowspans[2:end]==spans[2:end]
    end
    for entries in (NamedTuple[],[component("OFF",1;disabled=20,method="none")])
        out=sprint(io->DDSync.writeUncertaintySummary(io,entries))
        @test occursin("N/A",out) && !occursin("NaN%",out) && !occursin("# ADVICE:",out)
    end
    all_good=[component("OK",1;hutch=20,method="hutch",reason="converged")]
    out=sprint(io->DDSync.writeUncertaintySummary(io,all_good;prefix=""))
    @test occursin("Converged Hutch",out) && occursin("100.00%",out)
    @test !occursin("ADVICE:",out)
    @test DDSync.config_default()[:std][:probes]==4096
    # End-to-end accounting against actual quality-file rows and scale floor.
    mktempdir() do root
        cfg=split_config(root,mode="hutch")
        cfg[:std][:probes]=1;cfg[:std][:repair_max_nodes]=2
        execute(cfg)
        metrics=read(cfg[:io][:metrics_file],String)
        methods=String[]
        for path in readdir(cfg[:io][:thetastd_dir];join=true)
            startswith(basename(path),"std_quality_") || continue
            append!(methods,[split(line)[3] for line in readlines(path) if !startswith(line,"#")])
        end
        @test count(==("exact"),methods)==8
        @test count(==("pseudo_degree"),methods)==28
        @test occursin("36 supported non-reference",metrics)
        @test occursin("77.78%",metrics) && occursin("22.22%",metrics)
        @test occursin("# ADVICE:",metrics) && occursin("# std_component",metrics)
        @test !occursin("# std_component",first(split(metrics,"# COMPONENT DETAILS")))
    end
end
