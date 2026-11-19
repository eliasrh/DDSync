using Test, Random, DelimitedFiles
include(joinpath(@__DIR__, "..", "src", "DDSync.jl"))

function fixture(directory, ids)
    mkpath(directory)
    catalog=joinpath(directory,"catalog.txt")
    dt=joinpath(directory,"dt.cc")
    open(catalog,"w") do io
        println(io,"# Test catalog; whitespace and other columns must survive")
        for (i,id) in enumerate(ids)
            println(io,"2026 1 1 0 0 0   64.00",i," -20.00",i," 3.0 1.0\t",id)
        end
    end
    open(dt,"w") do io
        for i in 1:length(ids)-1, j in i+1:length(ids)
            println(io,"  # ",ids[i],"\t",ids[j],"  0.125")
            for s in 1:2, p in ("P","S")
                delta=(i-j)*0.013*s + sin(i*j)*0.0001
                println(io,"STA",s," ",delta," 0.9 ",p)
            end
        end
    end
    cfg=DDSync.config_default()
    for (key,name) in ((:catalog_file,"catalog.txt"),(:infile_dt,"dt.cc"),(:out_dt_sync,"dt_sync.cc"),
            (:metrics_file,"metrics.txt"),(:tmpdir,"tmp"),(:thetadir,"theta"),(:thetastd_dir,"thetastd"))
        cfg[:io][key]=joinpath(directory,name)
    end
    cfg[:robust][:MIN_EDGES]=3
    cfg[:std][:mode]="pseudo_degree"
    cfg[:std][:probes]=32
    cfg[:std][:batch]=16
    cfg
end

function execute(cfg)
    Random.seed!(12345)
    redirect_stdout(devnull) do
        DDSync.run(cfg)
    end
end

@testset "DDSyncJulia event IDs" begin
    @test DDSync.canonicalEventID("+00012","test")=="12"
    @test DDSync.canonicalEventID("-00012","test")=="-12"
    @test DDSync.canonicalEventID("-000","test")=="0"
    for id in ("1.0","1e9","NaN","a","+")
        @test_throws ErrorException DDSync.canonicalEventID(id,"test")
    end
    mktempdir() do root
        ids=["123456789012345678901234567890123456789", "9007199254740993", "16777217", "0", "-9007199254740993"]
        for (mode,scale,dense,pruned,stdmode) in (("base","fixed",false,false,"pseudo_degree"),
                ("robust","fixed",false,false,"pseudo_degree"), ("combined","fixed",true,true,"pseudo_degree"),
                ("thetaStd","fixed",false,false,"pseudo_degree"), ("thetaStd","median",false,false,"pseudo_degree"),
                ("thetaStd","median",true,true,"hutch"))
            label=join((mode,scale,dense,pruned,stdmode),'_')
            a=fixture(joinpath(root,label,"compact"),string.(1:5))
            b=fixture(joinpath(root,label,"large"),ids)
            for cfg in (a,b)
                cfg[:output][:weight_mode]=mode
                cfg[:output][:thetastd_scale_mode]=scale
                cfg[:output][:write_pruned_edges]=pruned
                cfg[:runtime][:store_dense_group_arrays]=dense
                cfg[:std][:mode]=stdmode
                execute(cfg)
            end
            adir=dirname(a[:io][:catalog_file]);bdir=dirname(b[:io][:catalog_file])
            @test !isfile(joinpath(adir,"event_id_map.csv"))
            @test !isfile(joinpath(adir,"catalog_seq.txt"))
            @test !isfile(joinpath(adir,"dt_sync_seq.cc"))
            @test read(joinpath(bdir,"event_id_map.csv"),String)=="original_id,sequential_id\n"*join(["$(ids[i]),$i\n" for i in 1:5])
            @test read(joinpath(bdir,"catalog_seq.txt"),String)==read(a[:io][:catalog_file],String)
            @test read(joinpath(bdir,"dt_sync_seq.cc"),String)==read(a[:io][:out_dt_sync],String)
            original=readlines(b[:io][:out_dt_sync]);sequential=readlines(joinpath(bdir,"dt_sync_seq.cc"))
            @test length(original)==length(sequential)>0
            for (x,y) in zip(original,sequential)
                if startswith(strip(x),"#")
                    fx,fy=split(x),split(y)
                    @test fx[2]==ids[parse(Int,fy[2])]
                    @test fx[3]==ids[parse(Int,fy[3])]
                    @test fx[4:end]==fy[4:end]
                else
                    @test x==y
                end
            end
            for name in readdir(a[:io][:thetadir])
                @test read(joinpath(a[:io][:thetadir],name))==read(joinpath(b[:io][:thetadir],name))
                @test length(readlines(joinpath(b[:io][:thetadir],name)))==5
            end
            for name in readdir(a[:io][:thetastd_dir])
                @test read(joinpath(a[:io][:thetastd_dir],name))==read(joinpath(b[:io][:thetastd_dir],name))
            end
            metrics=read(b[:io][:metrics_file],String)
            @test occursin("# event_ids_remapped=true",metrics)
            @test occursin("# dt_sync_event_ids=original",metrics)
        end
        # Shuffled 1:N identifiers remain unchanged; gaps map in catalog row order.
        for (name,values,expected) in (("shuffled",["3","1","5","2","4"],false),
                ("gapped",["99","7","123","21","81"],true))
            cfg=fixture(joinpath(root,name),values)
            index=DDSync.readEventIndex(cfg[:io][:catalog_file])
            @test index.remapped==expected
            @test [index.internal[v] for v in values]==(expected ? collect(1:5) : parse.(Int,values))
            execute(cfg)
        end
        # Catalog IDs are checked before computation; aliases represent the same integer.
        for values in (["1","1"],["01","+1"],["1","1.5"],["1","1e10"])
            cfg=fixture(joinpath(root,"invalid_catalog"),values)
            @test_throws ErrorException execute(cfg)
        end
        cfg=fixture(joinpath(root,"missing_id"),["10","20","30"])
        write(cfg[:io][:infile_dt],"# 10 40 0\nSTA1 0.1 0.9 P\n")
        @test_throws ErrorException execute(cfg)
        @test !isfile(joinpath(root,"missing_id","event_id_map.csv"))
        for text in ("# 10 10 0\nSTA1 0.1 0.9 P\n", "STA1 0.1 0.9 P\n", "# invalid\n")
            write(cfg[:io][:infile_dt],text)
            @test_throws ErrorException execute(cfg)
        end
        # Theta-only runs still provide the catalog and correspondence required by GraphSplit.
        cfg=fixture(joinpath(root,"theta_only"),ids)
        cfg[:output][:write_dt_sync]=false
        execute(cfg)
        @test isfile(joinpath(root,"theta_only","catalog_seq.txt"))
        @test isfile(joinpath(root,"theta_only","event_id_map.csv"))
        @test !isfile(joinpath(root,"theta_only","dt_sync.cc"))
        @test !isfile(joinpath(root,"theta_only","dt_sync_seq.cc"))
        # Retained zero-weight rows and entirely omitted blocks use matching IDs.
        for keep_pruned in (false,true)
            cfg=fixture(joinpath(root,"pruned_$keep_pruned"),ids)
            cfg[:robust][:MIN_EDGES]=100
            cfg[:output][:write_pruned_edges]=keep_pruned
            execute(cfg)
            original=readlines(cfg[:io][:out_dt_sync])
            sequential=readlines(joinpath(root,"pruned_$keep_pruned","dt_sync_seq.cc"))
            @test length(original)==length(sequential)
            @test isempty(original)==!keep_pruned
            for (a,b) in zip(original,sequential)
                if !startswith(strip(a),"#")
                    @test a==b
                    @test parse(Float64,split(a)[3])==0.0
                end
            end
        end
        # Custom output names, absent events, and output/input collision protection.
        cfg=fixture(joinpath(root,"custom"),ids)
        cfg[:io][:out_dt_sync]=joinpath(root,"custom","exports","cleaned.dat")
        open(cfg[:io][:catalog_file],"a") do io
            println(io,"2026 1 1 0 0 0 64 -20 3 1 888888888888888888888888888")
        end
        execute(cfg)
        @test isfile(joinpath(root,"custom","exports","cleaned_seq.dat"))
        theta=readdlm(joinpath(cfg[:io][:thetadir],"theta_STA1_P.txt"),Float64)
        @test size(theta,1)==6
        @test isnan(theta[end,2])
        cfg=fixture(joinpath(root,"collision"),ids)
        cfg[:io][:out_dt_sync]=cfg[:io][:infile_dt]
        before=read(cfg[:io][:infile_dt])
        @test_throws ErrorException execute(cfg)
        @test read(cfg[:io][:infile_dt])==before
    end
end

include("numerical_regressions.jl")

include("per_event_uncertainty.jl")

include("metrics_reporting.jl")

include("cholesky_uncertainty.jl")
