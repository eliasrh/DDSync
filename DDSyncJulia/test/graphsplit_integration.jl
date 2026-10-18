# Optional downstream integration test; GraphSplit is not a DDSync dependency.
# julia DDSyncJulia/test/graphsplit_integration.jl /path/to/GraphSplit.jl [MATLAB-test-output]
using Test, Random, DelimitedFiles, LinearAlgebra
length(ARGS)>=1 || error("Pass the path to a GraphSplit.jl checkout")
include(joinpath(abspath(ARGS[1]),"src","GraphSplit.jl"))
include(joinpath(@__DIR__,"..","src","DDSync.jl"))
include("fixture_config.jl")

function check_graphsplit(theta_dir,std_dir,fixture)
    groups=GraphSplit.load_theta_folder(theta_dir,std_dir)
    catalog=GraphSplit.read_catalog(joinpath(fixture,"catalog.txt"))
    mktemp() do path,io
        println(io,"TEST 64 -20 0");close(io)
        stations=GraphSplit.read_stations(path)
        cfg=GraphSplit.default_config()
        obs=GraphSplit.build_star_observations(groups,stations,catalog,cfg)
        @test length(obs)==36 # 18 non-reference nodes × 2 phases
        @test all(isfinite,obs.dt)
        @test all(s->isfinite(s)&&s>0,obs.sigma)
        @test all(k->(obs.i[k]<=10 && obs.j[k]==3) || (11<=obs.i[k]<=20 && obs.j[k]==13),eachindex(obs.i))
        # Exactly the two intended gauge freedoms; no missing-sigma disconnection.
        B=zeros(length(obs),20)
        for k in eachindex(obs.i)
            B[k,obs.i[k]]=1/obs.sigma[k];B[k,obs.j[k]]=-1/obs.sigma[k]
        end
        free=setdiff(1:20,[3,13])
        @test rank(B)==18
        @test isposdef(Symmetric(B[:,free]'*B[:,free]))
        for group in groups, id in group.event_id
            entry=group.entries[id]
            @test isfinite(entry.sigma)
            @test entry.sigma>0 || id==entry.ref_id
        end
    end
end

@testset "DDSync failure fallback through GraphSplit" begin
    mktempdir() do root
        for fallback in ("exact","graph_bound","nan")
            cfg=split_config(joinpath(root,fallback),mode="hutch")
            cfg[:std][:probes]=1 # deliberately cannot converge
            cfg[:std][:exact_max_nred]=0 # force finite graph bound
            cfg[:std][:fallback]=fallback
            Random.seed!(1234)
            redirect_stdout(devnull) do; DDSync.run(cfg); end
            check_graphsplit(cfg[:io][:thetadir],cfg[:io][:thetastd_dir],dirname(cfg[:io][:infile_dt]))
        end
    end
    if length(ARGS)>=2
        root=joinpath(abspath(ARGS[2]),"bound")
        fixture=joinpath(@__DIR__,"..","..","test","fixtures")
        check_graphsplit(joinpath(root,"theta"),joinpath(root,"thetastd"),fixture)
    end
end
