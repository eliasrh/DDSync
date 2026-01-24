#!/usr/bin/env julia
# run_ddsync.jl
#
# Runner for DDSync.jl with optional TOML configuration.
#
# Usage examples:
#   julia run_ddsync.jl
#   julia run_ddsync.jl dt.cc catalog.txt
#   julia run_ddsync.jl ddsync_config.toml
#   julia run_ddsync.jl ddsync_config.toml dt.cc catalog.txt
#
# Behavior:
#   - Always starts from DDSync.config_default().
#   - If a TOML file is provided (first argument ending in .toml) it is loaded
#     and deep-merged into the defaults.
#   - If no TOML file is provided but a file named "ddsync_config.toml" exists
#     in the current directory, it is loaded automatically.
#   - Remaining positional args override infile_dt and catalog_file.

include(joinpath(@__DIR__, "src", "DDSync.jl"))
using .DDSync
using TOML  # stdlib

"Recursively merge TOML Dicts into a Dict{Symbol,Any}."
function deep_merge!(dst::Dict{Symbol,Any}, src)
    # src is usually Dict{String,Any} (from TOML.parsefile)
    for (k,v) in src
        ks = Symbol(k)
        if v isa AbstractDict
            if haskey(dst, ks) && (dst[ks] isa Dict{Symbol,Any})
                deep_merge!(dst[ks], v)
            else
                # create new nested dict
                nd = Dict{Symbol,Any}()
                deep_merge!(nd, v)
                dst[ks] = nd
            end
        else
            dst[ks] = v
        end
    end
    return dst
end

cfg = DDSync.config_default()

# --- config file loading ---
argpos = 1
cfgfile = ""
if length(ARGS) >= 1 && endswith(lowercase(ARGS[1]), ".toml")
    cfgfile = ARGS[1]
    argpos = 2
elseif isfile("ddsync_config.toml")
    cfgfile = "ddsync_config.toml"
end

if !isempty(cfgfile)
    @info("Loading config", file=cfgfile)
    toml = TOML.parsefile(cfgfile)
    deep_merge!(cfg, toml)
end

# --- positional overrides for dt + catalog ---
if length(ARGS) >= argpos
    cfg[:io][:infile_dt] = ARGS[argpos]
    argpos += 1
end
if length(ARGS) >= argpos
    cfg[:io][:catalog_file] = ARGS[argpos]
    argpos += 1
end

DDSync.run(cfg)
