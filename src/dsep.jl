# Blocking and d-separation


"""
    dsep(g::AbstractGraph, u, v, s; verbose = false)

Check  whether `u` and `v` are d-separated given set `s`.
Algorithm: unrolled https://arxiv.org/abs/1304.1505

"""
function dsep(g::AbstractGraph, U, V, S; verbose = false)
    T = eltype(g)
    in_seen = falses(nv(g)) # nodes reached earlier backwards
    out_seen = falses(nv(g)) # nodes reached earlier forwards
    descendant = falses(nv(g)) # descendant in s
    isv = falses(nv(g))
    blocked = falses(nv(g))

    for ve in S
        in_seen[ve] = true
        blocked[ve] = true
    end

    next = Vector{T}()

    # mark vertices with descendants in S
    next = Vector{T}()
    for w in S
        push!(next, w)
        descendant[w] = true
    end

    while !isempty(next)
        for w in inneighbors(g, popfirst!(next))
            if !descendant[w]
                push!(next, w) # push onto queue
                descendant[w] = true
            end
        end
    end

    verbose && println(descendant)

    in_next = Vector{T}()
    out_next = Vector{T}()

    for u in U
        push!(in_next, u) # treat u as vertex reached backwards
        in_seen[u] && throw(ArgumentError("U and S not disjoint."))
        in_seen[u] = true
    end
    if V isa Integer
        in_seen[V] && throw(ArgumentError("U, V and S not disjoint."))
        return dsep_inner!(g, in_next, out_next, descendant, ==(V), blocked, out_seen, in_seen; verbose)
    else
        isv = falses(nv(g))
        for v in V
            in_seen[v] && throw(ArgumentError("U, V and S not disjoint."))
            isv[v] = true
        end
        return dsep_inner!(g, in_next, out_next, descendant, w->isv[w], blocked, out_seen, in_seen; verbose)
    end
end

function dsep_inner!(g, in_next, out_next, descendant, found, blocked, out_seen, in_seen; verbose=false)
    while true
        sin = isempty(in_next)
        sout = isempty(out_next)
        sin && sout && return true # no vertices in the queue

        if !sin # some vertex reach backwards in the queue
            src = popfirst!(in_next)
            for w in outneighbors(g, src) # possible collider at destination
                if !out_seen[w] && (!blocked[w] || descendant[w])
                    verbose && println("<- $src -> $w")
                    found(w) && return false
                    push!(out_next, w)
                    out_seen[w] = true
                end
            end
            for w in inneighbors(g, src)
                if !in_seen[w]
                    verbose && println("<- $src <- $w")
                    found(w) && return false
                    push!(in_next, w)
                    in_seen[w] = true
                end
            end
        end
        if !sout # some vertex reach forwards in the queue
            src = popfirst!(out_next)
            for w in outneighbors(g, src) # possible collider at destination
                if !out_seen[w] && !blocked[src] && (!blocked[w] || descendant[w])
                    verbose && println("-> $src -> $w")
                    found(w) && return false
                    push!(out_next, w)
                    out_seen[w] = true
                end
            end
            for w in inneighbors(g, src) # collider at source
                if !in_seen[w] && descendant[src] # shielded collider
                    verbose && println("-> $src <- $w")
                    found(w) && return false
                    push!(out_next, w)
                    in_seen[w] = true
                end
            end
        end
    end
end
