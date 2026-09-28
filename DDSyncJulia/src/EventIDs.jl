# External IDs remain decimal strings; only compact internal IDs enter numeric arrays.
struct EventIndex
    original::Vector{String}          # catalog row order
    internal::Dict{String,Int}        # canonical external ID -> internal ID
    remapped::Bool
end

function canonicalEventID(token::AbstractString, location::AbstractString)
    occursin(r"^[+-]?[0-9]+$", token) ||
        error("$location: event ID '$token' must be a decimal integer (not scientific notation).")
    negative = startswith(token, "-")
    digits = first(token) in ('+', '-') ? token[2:end] : token
    digits = lstrip(digits, '0')
    isempty(digits) && return "0"
    return negative ? "-" * digits : String(digits)
end

catalogComment(line) = isempty(strip(line)) || first(strip(line)) in ('#', '%', '*', '!')

function readEventIndex(catalog_file::String)
    original = String[]
    internal = Dict{String,Int}()
    for (line, text) in enumerate(eachline(catalog_file))
        catalogComment(text) && continue
        id = canonicalEventID(last(split(text)), "$catalog_file:$line")
        haskey(internal, id) && error("$catalog_file:$line: duplicate event ID $id.")
        # The existing group spool stores compact indices as Float32.
        length(original) < 2^24 || error("DDSyncJulia supports at most 16,777,216 catalog events in its current spool format.")
        push!(original, id)
        internal[id] = length(original)
    end
    isempty(original) && error("No event IDs found in $catalog_file.")
    n = length(original)
    sequential = all(original) do id
        value = tryparse(Int, id)
        value !== nothing && 1 <= value <= n
    end
    # Unique IDs all in 1:N already form the complete range, even if rows are shuffled.
    if sequential
        for id in original
            internal[id] = parse(Int, id)
        end
    end
    return EventIndex(original, internal, !sequential)
end

function mappedPair(header::AbstractString, index::EventIndex, location::AbstractString)
    fields = split(strip(header[2:end]))
    length(fields) >= 2 || error("$location: expected '# event1 event2 [origin-time term]'.")
    a = canonicalEventID(fields[1], location)
    b = canonicalEventID(fields[2], location)
    i = get(index.internal, a, 0)
    j = get(index.internal, b, 0)
    i > 0 || error("$location: event ID $a is absent from the catalog.")
    j > 0 || error("$location: event ID $b is absent from the catalog.")
    i != j || error("$location: event pair $a $b refers to the same event twice.")
    return i, j
end

function sequentialHeader(header::AbstractString, i::Int, j::Int)
    m = match(r"^(\s*#\s*)\S+(\s+)\S+(.*)$", header)
    m === nothing && error("Invalid differential-time header: $header")
    return string(m[1], i, m[2], j, m[3])
end

# Retain the standalone spool helper for compact, already-reindexed inputs.
function mappedPair(header::AbstractString, ::Nothing, location::AbstractString)
    fields = split(strip(header[2:end]))
    length(fields) >= 2 || error("$location: invalid event-pair header.")
    ids = map(fields[1:2]) do token
        id = tryparse(Int, canonicalEventID(token, location))
        id !== nothing && 1 <= id <= 2^24 || error("$location: use DDSync.run to map external IDs before spooling.")
        id
    end
    ids[1] != ids[2] || error("$location: self-pair.")
    return ids[1], ids[2]
end

function indexOutputPaths(out_dt_sync::String)
    stem, ext = splitext(out_dt_sync)
    return (catalog=joinpath(dirname(out_dt_sync), "catalog_seq.txt"),
            mapping=joinpath(dirname(out_dt_sync), "event_id_map.csv"),
            dt=stem * "_seq" * ext)
end

function checkIndexOutputPaths(paths, inputs::Vector{String})
    # Resolve existing symlinks too: generated files must never replace an input.
    resolved(p) = ispath(p) ? realpath(p) : joinpath(realpath(dirname(abspath(p))), basename(p))
    destinations = [resolved(p) for p in paths]
    length(unique(destinations)) == length(destinations) || error("DDSync output paths overlap.")
    for source in inputs
        resolved(source) in destinations && error("An output would overwrite input $source; choose different output paths.")
    end
end

function writeIndexFiles(catalog_file::String, index::EventIndex, paths)
    open(paths.mapping, "w") do io
        println(io, "original_id,sequential_id")
        for id in index.original
            println(io, id, ',', index.internal[id])
        end
    end
    open(paths.catalog, "w") do out
        for (line, text) in enumerate(eachline(catalog_file; keep=true))
            if catalogComment(text)
                write(out, text)
            else
                m = match(r"^(.*?)(\S+)(\s*)$"s, text)
                id = canonicalEventID(m[2], "$catalog_file:$line")
                print(out, m[1], index.internal[id], m[3])
            end
        end
    end
end

function writeIndexMetrics(io, index::EventIndex, paths, write_dt::Bool)
    println(io, "# event_count=", length(index.original))
    println(io, "# event_ids_remapped=", index.remapped)
    println(io, "# theta_event_ids=", index.remapped ? "sequential" : "original (already 1:N)")
    println(io, "# dt_sync_event_ids=original")
    if index.remapped
        println(io, "# mapping_order=catalog row order")
        println(io, "# event_id_map=", abspath(paths.mapping))
        println(io, "# graphsplit_catalog=", abspath(paths.catalog))
        write_dt && println(io, "# sequential_dt_sync=", abspath(paths.dt))
    end
end
