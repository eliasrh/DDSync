using LinearAlgebra, SparseArrays

include("fixture_config.jl")

@testset "Post-pruning references and safe exports" begin
    mktempdir() do root
        for irls in (0,10), dense in (false,true)
            cfg=split_config(joinpath(root,"split_$(irls)_$dense"))
            cfg[:irls][:iters]=irls
            cfg[:runtime][:store_dense_group_arrays]=dense
            execute(cfg)
            for phase in ("P","S")
                theta=readdlm(joinpath(cfg[:io][:thetadir],"theta_TEST_$phase.txt"),Float64)
                sd=readdlm(joinpath(cfg[:io][:thetastd_dir],"std_theta_TEST_$phase.txt"),Float64)
                @test theta[1:10,3]==fill(3,10)
                @test theta[11:20,3]==fill(13,10)
                @test theta[[3,13],2]==[0,0]
                @test all(isnan,theta[21:23,2:3])
                @test all(isnan,sd[21:23,2:3])
                @test sd[21:23,4]==[0,0,0]
                @test all(sd[setdiff(1:20,[3,13]),2].>1e-5)
                @test sd[[1,2,11,12],4]==[8,8,8,8]
                @test sd[[3,13],4]==[9,9]
                @test maximum(abs.(theta[1:10,2] .- 0.013.*((1:10).-3).*(phase=="P" ? 1 : 1.7)))<0.0002
            end
            lines=readlines(cfg[:io][:out_dt_sync])
            @test count(l->startswith(l,"#"),lines)==88
            for line in filter(l->startswith(l,"#"),lines)
                fields=split(line); i=parse(Int,fields[2]); j=parse(Int,fields[3])
                @test (i<=10 && j<=10) || (11<=i<=20 && 11<=j<=20)
            end
            @test occursin("ref=3 nodes=10 method=exact",read(cfg[:io][:metrics_file],String))
        end
        # Surviving fragments must pass the gate themselves, not borrow old redundancy.
        cfg=split_config(joinpath(root,"gate"));cfg[:robust][:MIN_EDGES]=45
        execute(cfg)
        @test isempty(read(cfg[:io][:out_dt_sync],String))
        @test all(isnan,readdlm(joinpath(cfg[:io][:thetadir],"theta_TEST_P.txt"),Float64)[:,2:3])
        # Nonconvergence triggers exact fallback, not tiny standard deviations.
        cfg=split_config(joinpath(root,"fallback"),mode="hutch");cfg[:std][:probes]=1
        execute(cfg)
        metrics=read(cfg[:io][:metrics_file],String)
        @test occursin("method=exact",metrics)
        @test occursin("probes=1",metrics)
        sd=readdlm(joinpath(cfg[:io][:thetastd_dir],"std_theta_TEST_P.txt"),Float64)
        @test minimum(sd[setdiff(1:20,[3,13]),2])>1e-5
        # Exceeding the exact size cap preserves finite conservative uncertainties.
        cfg=split_config(joinpath(root,"bound"),mode="hutch")
        cfg[:std][:probes]=1;cfg[:std][:exact_max_nred]=0
        execute(cfg)
        sd=readdlm(joinpath(cfg[:io][:thetastd_dir],"std_theta_TEST_P.txt"),Float64)
        @test all(isfinite,sd[1:20,2])
        @test all(sd[setdiff(1:20,[3,13]),2].>0)
        @test all(sd[setdiff(1:20,[3,13]),5].>0)
        @test all(isnan,sd[21:23,2])
        @test occursin("exact_size_limit",read(cfg[:io][:metrics_file],String))
        @test all(l->parse(Float64,split(l)[3])>0,filter(l->!startswith(l,"#"),readlines(cfg[:io][:out_dt_sync])))
        @test occursin("method=graph_bound",read(cfg[:io][:metrics_file],String))
        # A tree separated from cyclic components has no residual redundancy.
        cfg=split_config(joinpath(root,"split_tree"))
        cfg[:io][:infile_dt]=joinpath(dirname(cfg[:io][:infile_dt]),"split_tree.cc")
        cfg[:robust][:MIN_EDGES]=3
        execute(cfg)
        theta=readdlm(joinpath(cfg[:io][:thetadir],"theta_TEST_P.txt"),Float64)
        @test theta[1:3,3]==fill(1,3)
        @test theta[4:6,3]==fill(4,3)
        @test all(isnan,theta[7:23,2:3])
        # Mixed exact/bounded uncertainties preserve finite median weights.
        cfg=split_config(joinpath(root,"median"))
        cfg[:io][:infile_dt]=joinpath(dirname(cfg[:io][:infile_dt]),"median.cc")
        cfg[:robust][:MIN_EDGES]=3
        cfg[:std][:exact_max_nred]=2
        cfg[:output][:thetastd_scale_mode]="median"
        execute(cfg)
        for line in filter(l->!startswith(l,"#"),readlines(cfg[:io][:out_dt_sync]))
            fields=split(line);weight=parse(Float64,fields[3])
            @test weight>0.5
        end
        cfg=split_config(joinpath(root,"columns"),mode="pseudo_weight")
        cfg[:std][:thetastd_write_weightcol]=false
        cfg[:std][:thetastd_write_alt_node_weightcol]=true
        execute(cfg)
        sd=readdlm(joinpath(cfg[:io][:thetastd_dir],"std_theta_TEST_P.txt"),Float64)
        @test size(sd,2)==5
        @test all(abs.(sd[1:20,5].-0.9).<1e-6)
        cfg=split_config(joinpath(root,"no_export"))
        cfg[:std][:export]=false;cfg[:output][:weight_mode]="combined"
        execute(cfg)
        @test !isdir(cfg[:io][:thetastd_dir])
        @test !isempty(read(cfg[:io][:out_dt_sync],String))
        # Explicit proxies, disabled/size-gated Hutch, no-std modes, retained pruned rows.
        for (mode,fallback,probes,max_nred) in (("hutch","pseudo_degree",0,40000),("hutch","pseudo_weight",128,0),("none","exact",128,40000))
            cfg=split_config(joinpath(root,"$mode-$fallback"),mode=mode)
            cfg[:std][:fallback]=fallback;cfg[:std][:probes]=probes;cfg[:std][:max_nred]=max_nred
            cfg[:output][:write_pruned_edges]=true
            execute(cfg)
            @test count(l->startswith(l,"#"),readlines(cfg[:io][:out_dt_sync]))==95
            sd=readdlm(joinpath(cfg[:io][:thetastd_dir],"std_theta_TEST_P.txt"),Float64)
            @test mode=="none" ? all(isnan,sd[setdiff(1:20,[3,13]),2]) : all(isfinite,sd[1:20,2])
        end
    end
end

@testset "Uncertainty against exact inverse and failure diagnostics" begin
    a=Int[]; b=Int[]
    for i in 1:9,j in i+1:10
        push!(a,i);push!(b,j)
    end
    w=0.5 .+ collect(1:length(a))./length(a)
    free=collect(2:10);ridge=1e-10
    # Independent incidence-matrix oracle.
    A=zeros(length(a),10)
    for k in eachindex(a);A[k,a[k]]=1;A[k,b[k]]=-1;end
    expected=sqrt.(diag(inv(A[:,free]'*Diagonal(w)*A[:,free]+ridge*I)))
    sd,ok=DDSync.estimateThetaStdExact(10,a,b,w,1,free,ridge,1.0,0.0,false,3)
    @test ok
    @test sd[free]≈expected rtol=1e-12
    bound=DDSync.estimateThetaStdGraphBound(10,a,b,w,1,1.0)
    @test all(bound[free].>=expected.-1e-12)
    @test all(isfinite,bound) && bound[1]==0
    tree=DDSync.estimateThetaStdGraphBound(4,[1,2,3],[2,3,4],[2.0,4.0,0.5],1,1.0)
    @test tree≈sqrt.([0.0,0.5,0.75,2.75])
    @test_throws ErrorException DDSync.estimateThetaStdGraphBound(3,[1],[2],[1.0],1,1.0)
    @test_throws ErrorException DDSync.estimateThetaStdGraphBound(10,a,b,w,1,NaN)
    for dist in ("rademacher","gaussian")
        Random.seed!(42)
        sd,info=DDSync.estimateThetaStdHutch(10,a,b,w,1,free,ridge,1.0,0.0,false,8192,dist,128,0,1e-12;rel_tol=0.05)
        @test info.reason=="converged"
        @test 128<=info.probes<8192
        @test info.max_rel_se<=0.05
        @test maximum(abs.(sd[free]./expected.-1))<0.12
    end
    @test DDSync.hutchDiagnostic([-1.0,1.0],[1.0,1.0],128,128,0.1)[1]=="nonpositive_estimate"
    @test DDSync.hutchDiagnostic([NaN,1.0],[1.0,1.0],128,128,0.1)[1]=="nonfinite_estimate"
    @test DDSync.hutchDiagnostic([1.0],[1.0],1,128,0.1)[1]=="insufficient_probes"
    @test DDSync.hutchDiagnostic([1.0],[1e8],128,128,0.1)[1]=="not_converged"
    # Singular numerical system is a failure, not a positive precision floor.
    sd,ok=DDSync.estimateThetaStdExact(3,[1],[2],[1.0],1,[2,3],0.0,1.0,0.0,false,2)
    @test !ok && isnan(sd[2]) && isnan(sd[3])
    @test isnan(DDSync.sigmaForThetaStd(NaN,0.0005,true))
    for (key,value) in ((:batch,0),(:min_probes,1),(:rel_tol,NaN),(:probes,-1),(:mode,"typo"))
        cfg=DDSync.config_default()[:std];cfg[key]=value
        @test_throws Exception DDSync.validateStdConfig(cfg)
    end
end
