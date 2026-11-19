@testset "Shared timing variance and adaptive Cholesky draws" begin
    # The weak link lets events on its far side move together. Force a nontrivial
    # factor order so that comparing to the dense inverse checks event identity.
    d = [0.4, 2.0, 0.01, 3.0, 1.5, 0.8]
    incidence = Matrix{Float64}(I,6,6) - diagm(-1=>ones(5))
    A = sparse(incidence' * Diagonal(d) * incidence)
    F = cholesky(Symmetric(A);perm=[5,3,4,2,6,1])
    truth = diag(inv(Matrix(A)))
    for dist in ("rademacher","gaussian")
        Random.seed!(309)
        result = DDSync.choleskySamples(A,F,1024,dist,37,0,128,0.25)
        @test result.probes == 256 && all(result.converged)
        @test all(abs.(result.mu ./ truth .- 1) .< 0.25)
    end
    Random.seed!(17)
    limited = DDSync.choleskySamples(A,F,613,"gaussian",37,0,128,1e-12)
    @test limited.probes == 613 && limited.reason == "not_converged"
    single = sparse(reshape([4.],1,1))
    shared = DDSync.choleskySamples(single,cholesky(single),600,"gaussian",19,0,128,0.1)
    @test shared.probes == 256 && shared.mu ≈ [0.25] && shared.rel_se[1] < 1e-14
    cfg = DDSync.config_default()[:std]
    @test cfg[:auto_probe] == "hutch"
    @test_throws ErrorException DDSync.validateStdConfig(merge(cfg,Dict(:auto_probe=>"invalid")))
    opts = merge((;(k=>v for (k,v) in cfg)...),
        (mode="auto",auto_probe="cholesky",exact_max_nred=1,rel_tol=1e-12,
         probes=300,min_probes=128,batch=41,repair_max_nodes=0,report_every_batch=0))
    edges = (7,collect(1:6),collect(2:7),d,1,collect(2:7),1e-10,0.01)
    Random.seed!(19)
    sd,info = DDSync.estimateThetaStd(edges...,opts)
    @test all(isfinite,sd) && info.probes == 300 && info.node_methods[1]=="reference"
    @test info.cholesky_unconverged>0 && info.cholesky_unconverged==info.unconverged_nodes
    @test info.fallback_nodes==0 && info.hutch_converged==0
    c = DDSync.uncertaintyCounts([merge(info,(nodes=7,))])
    @test c.total==6 && c.review==info.cholesky_unconverged
    _,skipped = DDSync.estimateThetaStd(edges...,merge(opts,(max_nred=5,)))
    @test skipped.probes==0 && skipped.reason=="size_limit" && skipped.fallback_nodes==6
    mktempdir() do root
        cfg=split_config(root,mode="auto")
        cfg[:std][:auto_probe]="cholesky"
        cfg[:std][:exact_max_nred]=1
        execute(cfg)
        report=read(cfg[:io][:metrics_file],String)
        @test occursin("Converged Cholesky",report) && occursin("probes=256",report)
        @test occursin(" cholesky ",read(joinpath(cfg[:io][:thetastd_dir],"std_quality_TEST_P.txt"),String))
    end
end
