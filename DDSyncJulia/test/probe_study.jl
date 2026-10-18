# Reproducible sampling study, separate from fast regression tests.
# Run from the repository root: julia DDSyncJulia/test/probe_study.jl
using Random, LinearAlgebra, Printf
include(joinpath(@__DIR__,"..","src","DDSync.jl"))
println("graph,budget,accepted,seeds,median_used,max_relative_std_error_accepted")
for graph in ("dense20","dense100","ring100","weak_bridge100","hub100")
    n=graph=="dense20" ? 20 : 100
    a=Int[];b=Int[];w=Float64[]
    for i in 1:n-1,j in i+1:n
        keep = startswith(graph,"dense") || (graph=="ring100" && (j==i+1 || (i==1 && j==n))) ||
            (graph=="hub100" && (i==1 || j==i+1)) ||
            (graph=="weak_bridge100" && ((i<=50 && j<=50) || (i>50 && j>50) || (i==50 && j==51)))
        keep || continue
        push!(a,i);push!(b,j);push!(w,graph=="weak_bridge100" && i==50 && j==51 ? 1e-3 : 1.0)
    end
    pin=argmax(DDSync.nodeDegrees(n,a,b));free=setdiff(1:n,[pin])
    oracle,ok=DDSync.estimateThetaStdExact(n,a,b,w,pin,free,1e-10,1.0,0.0,false,128)
    @assert ok
    for budget in (128,512,2048,8192)
        accepted=0;used=Int[];errors=Float64[]
        for seed in 1:10
            Random.seed!(seed)
            sd,info=DDSync.estimateThetaStdHutch(n,a,b,w,pin,free,1e-10,1.0,0.0,false,budget,"rademacher",128,0,1e-12)
            push!(used,info.probes)
            if info.reason=="converged"
                accepted+=1;push!(errors,maximum(abs.(sd[free]./oracle[free].-1)))
            end
        end
        @printf("%s,%d,%d,10,%.0f,%.6f\n",graph,budget,accepted,DDSync.median(used),isempty(errors) ? NaN : maximum(errors))
    end
end
