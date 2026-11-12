@testset "Per-event uncertainty and selected exact repair" begin
    # A well-constrained hub with a weakly connected dense tail. Most entries
    # converge quickly; the tail's marginal covariance includes the weak bridge.
    n=300; a=Int[]; b=Int[]; w=Float64[]
    edge(i,j,v=1.0)=(push!(a,i);push!(b,j);push!(w,v))
    for j=2:260;edge(1,j);j>2 && edge(j-1,j);end
    for i=261:299,j=i+1:300;edge(i,j);end
    edge(2,261,0.001)
    free=collect(2:n);ridge=1e-10
    A=zeros(length(a),n)
    for k in eachindex(a);A[k,a[k]]=1;A[k,b[k]]=-1;end
    L=A[:,free]'*Diagonal(w)*A[:,free]+ridge*I
    expected=sqrt.(diag(inv(L)))
    cfg=DDSync.config_default()[:std]
    opts=(;(k=>v for (k,v) in cfg)...)
    opts=merge(opts,(mode="hutch",probes=512,report_every_batch=0,repair_max_nodes=64,exact_max_nred=10))
    Random.seed!(42)
    raw,sample=DDSync.estimateThetaStdHutch(n,a,b,w,1,free,ridge,1.0,0.0,false,512,"rademacher",128,0,0.0)
    @test sum(sample.converged)>250
    @test 1<=sum(.!sample.converged)<=64
    Random.seed!(42)
    sd,info=DDSync.estimateThetaStd(n,a,b,w,1,free,ridge,1.0,opts)
    good=free[sample.converged];bad=free[.!sample.converged]
    @test sd[good]==raw[good] # another event failing must not alter these
    @test all(info.node_methods[good].=="hutch")
    @test all(info.node_methods[bad].=="exact")
    @test info.repaired_nodes==length(bad) && info.fallback_nodes==0
    @test sd[bad]≈expected[bad.-1] rtol=1e-8
    @test info.method=="hybrid"
    # Full-system selected solves, not a local neighborhood inverse.
    selected=DDSync.selectedInverseDiagonal(sparse(L),cholesky(Symmetric(sparse(L))),[299,2],2)
    @test sqrt.(selected)≈expected[[299,2]] rtol=1e-8
    @test selected[1]>100 # local degree proxy would be approximately 1/39
    # Checked columns survive even if a different selected solve is unusable.
    partial=DDSync.selectedInverseDiagonal(Matrix{Float64}(I,3,3),Diagonal([1.0,2.0,1.0]),[1,2,3],3)
    @test partial[[1,3]]==[1.0,1.0] && isnan(partial[2])
    # Numerically negligible remote entries can underflow without invalidating
    # an accurate selected diagonal. The infinite-chain endpoint has this limit.
    chain=spdiagm(-1=>fill(-1.0,899),0=>fill(3.0,900),1=>fill(-1.0,899))
    endpoint=DDSync.selectedInverseDiagonal(chain,cholesky(Symmetric(chain)),[1],1)
    @test endpoint[1]≈(3-sqrt(5))/2 rtol=1e-12
    @test DDSync.hutchNodeDiagnostics([-1.0,2.0],[1.0,0.5],128,128,0.1)[1]==[false,true]
    # Zero repair budget: retain positive unresolved values with clear labels.
    Random.seed!(42)
    kept,ki=DDSync.estimateThetaStd(n,a,b,w,1,free,ridge,1.0,merge(opts,(repair_max_nodes=0,)))
    @test kept[good]==raw[good]
    @test ki.unconverged_nodes>0 && ki.repair_attempted==0
    for i in free
        if ki.node_methods[i]=="hutch_unconverged"
            @test kept[i]==raw[i] && 0.1<ki.node_rel_se[i]<1
        end
    end
    Random.seed!(42)
    proxy,pi=DDSync.estimateThetaStd(n,a,b,w,1,free,ridge,1.0,merge(opts,(repair_max_nodes=0,unconverged="fallback")))
    @test proxy[good]==raw[good]
    @test pi.fallback_nodes==length(bad) && pi.unconverged_nodes==0
    @test proxy[bad]≈1 ./sqrt.(DDSync.nodeDegrees(n,a,b)[bad])
    # Very few probes cannot be called converged or retained as usable sampling.
    Random.seed!(1)
    low,li=DDSync.estimateThetaStd(n,a,b,w,1,free,ridge,1.0,merge(opts,(probes=1,repair_max_nodes=2)))
    @test li.exact_nodes==2 && li.fallback_nodes==n-3 && li.unconverged_nodes==0
    @test all(isfinite,low) && all(low[free].>0)
    # Auto does no random work below the cap; preserve explicit hutch behavior.
    Random.seed!(10);next_random=rand();Random.seed!(10)
    exact,ei=DDSync.estimateThetaStd(n,a,b,w,1,free,ridge,1.0,merge(opts,(mode="auto",exact_max_nred=400)))
    @test rand()==next_random
    @test ei.method=="exact" && ei.probes==0 && ei.exact_nodes==n-1
    @test exact[free]≈expected rtol=1e-8
    Random.seed!(42)
    large,ai=DDSync.estimateThetaStd(n,a,b,w,1,free,ridge,1.0,merge(opts,(mode="auto",)))
    @test large==sd && ai.method=="hybrid"
    skip,si=DDSync.estimateThetaStd(n,a,b,w,1,free,ridge,1.0,merge(opts,(max_nred=0,)))
    @test si.probes==0 && si.method=="pseudo_degree" && all(isfinite,skip)
    skipped,sk=DDSync.estimateThetaStd(n,a,b,w,1,free,ridge,1.0,merge(opts,(probes=0,fallback="exact",exact_max_nred=400)))
    @test sk.method=="exact" && sk.reason=="disabled_probes" && sk.probes==0
    @test skipped[free]≈expected rtol=1e-8
    # Invalid scale must never be silently replaced with a finite precision floor.
    @test_throws ErrorException DDSync.estimateThetaStd(n,a,b,w,1,free,ridge,NaN,opts)
    for (key,value) in ((:repair_max_nodes,-1),(:unconverged,"typo"))
        c=copy(cfg);c[key]=value
        @test_throws ErrorException DDSync.validateStdConfig(c)
    end
    mktempdir() do root
        c=split_config(root,mode="auto")
        execute(c)
        rows=filter(l->!startswith(l,"#"),readlines(joinpath(c[:io][:thetastd_dir],"std_quality_TEST_P.txt")))
        @test length(rows)==20
        @test count(l->split(l)[3]=="exact",rows)==18
        @test count(l->split(l)[3]=="reference",rows)==2
        @test occursin("exact_nodes=9",read(c[:io][:metrics_file],String))
    end
end
