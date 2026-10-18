function split_config(root; mode="exact")
    cfg=DDSync.config_default()
    fixture=joinpath(@__DIR__,"..","..","test","fixtures")
    cfg[:io][:infile_dt]=joinpath(fixture,"split.cc")
    cfg[:io][:catalog_file]=joinpath(fixture,"catalog.txt")
    for (key,name) in ((:out_dt_sync,"dt_sync.cc"),(:metrics_file,"metrics.txt"),(:tmpdir,"tmp"),(:thetadir,"theta"),(:thetastd_dir,"thetastd"))
        cfg[:io][key]=joinpath(root,name)
    end
    cfg[:robust][:K_SIGMA]=4.0
    cfg[:std][:mode]=mode
    cfg[:std][:report_every_batch]=0
    cfg[:output][:dt_decimals]=9
    cfg[:output][:theta_decimals]=12
    cfg[:output][:thetastd_decimals]=12
    cfg
end
