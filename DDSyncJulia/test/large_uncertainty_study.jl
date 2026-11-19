# Large-component study: Julia, three seeds, 2048 probes, unit residual scale.
# Timings include JIT on first use and are illustrative, not a benchmark.
using LinearAlgebra,Random,Statistics,Printf
include(joinpath(@__DIR__,"..","src","DDSync.jl"))
BLAS.set_num_threads(2)
function graph_data(name)
    n=name=="dense2400" ? 2400 : 5000
    a=Int[];b=Int[];w=Float64[]
    edge(i,j,v=1.0)=(push!(a,i);push!(b,j);push!(w,v))
    if startswith(name,"dense")
        for i=1:n-1,j=i+1:n;edge(i,j);end
    elseif name=="mostly_good5000"
        split=4960
        for j=2:split;edge(1,j);j>2 && edge(j-1,j);end
        for i=split+1:n-1,j=i+1:n;edge(i,j);end
        edge(2,split+1,0.001)
    elseif name=="weak_bridge5000"
        for start in (1,2501),j=start+1:start+2499
            edge(start,j);j>start+1 && edge(j-1,j)
        end
        edge(1,2501,1.0)
    end
    n,a,b,w
end
function study()
 println("graph,seed,policy,seconds,converged,exact,kept,fallback,median_std_ratio,p95_abs_std_error,max_abs_std_error,graph_bound_median_ratio,graph_bound_p95_ratio,reason,probes")
 for name in ("mostly_good5000","weak_bridge5000","dense2400")
    n,a,b,w=graph_data(name);free=collect(2:n);ridge=1e-10
    oracle,ok=DDSync.estimateThetaStdExact(n,a,b,w,1,free,ridge,1.0,0.0,false,128);@assert ok
    name=="dense2400" && @assert maximum(abs.(oracle[free]./sqrt(2/n).-1))<1e-6
    bound=DDSync.estimateThetaStdGraphBound(n,a,b,w,1,1.0)
    br=bound[free]./oracle[free]
    base=DDSync.config_default()[:std];base[:report_every_batch]=0;base[:probes]=2048 # historical November 12 experiment
    for seed in 1:3, policy in ("keep","fallback")
        opts=(; (k=>v for (k,v) in base)...)
        opts=merge(opts,(mode="hutch",unconverged=policy))
        Random.seed!(seed)
        t=@elapsed sd,info=DDSync.estimateThetaStd(n,a,b,w,1,free,ridge,1.0,opts)
        ratio=sd[free]./oracle[free];err=abs.(ratio.-1)
        @printf("%s,%d,%s,%.3f,%d,%d,%d,%d,%.6g,%.6g,%.6g,%.6g,%.6g,%s,%d\n",name,seed,policy,t,info.hutch_converged,info.exact_nodes,info.unconverged_nodes,info.fallback_nodes,median(ratio),quantile(err,.95),maximum(err),median(br),quantile(br,.95),info.reason,info.probes);flush(stdout)
    end
 end
end
study()
