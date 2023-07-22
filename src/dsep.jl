# Blocking and d-seperation

isblocked(g, x, y, nodesRemoved) = !has_a_path(g, [x], y, nodesRemoved)

#Do the nodesRemoved block all semi-directed paths between src and dest?
#= """
    isblocked(g, src, dest, nodesRemoved)

Return `true` if there is no semi-directed path between `src` and `dest` in the graph `g`.
A set of vertices (`nodesRemoved`) can be removed from the graph before searching for a semi-directed path.
A semi-directed path between `src` and `dest` is a list of edges in `g` where every edge is either undirected or points toward `dest`. 
    src → x₁ - x₂ → dest ✓
    src → x₁ ← x₂ - dest ✖
"""
function isblocked(g, x, y, nodesRemoved)

    # Keep track of all the nodes visited
    visited = zeros(Bool, nv(g))

    # Mark excluded vertices as visited
    for vᵢ in nodesRemoved 
        visited[vᵢ] = true
    end

    # If src or dest were in nodesRemoved, the path is blocked
    (visited[x] || visited[y]) && return true

    # If the src and dest are the same, path is itself
    x == y && return false

    # Check if scr or dest have no neighbors
    if hasnoneighbors(g, x) || hasnoneighbors(g, y)
        return true
    end

    queue = [x]
    visited[x] = true
    while !isempty(queue)
        current = popfirst!(queue) # get new element from queue
        for vᵢ in descendent(g, current)
            vᵢ == y && return false
            if !visited[vᵢ]
                push!(queue, vᵢ) # push onto queue
                visited[vᵢ] = true
            end
        end
    end
    return true
end
=#

"""
    dsep(g::AbstractGraph, u, v, s; verbose = false)

Check  whether `u` and `v` are d-separated given set `s`.
Algorithm: unrolled https://arxiv.org/abs/1304.1505

"""
function dsep(g::AbstractGraph, u::Integer, v::Integer, S; verbose = false)
    T = eltype(g)
    in_seen = falses(nv(g)) # nodes reached earlier backwards
    out_seen = falses(nv(g)) # nodes reached earlier forwards
    descendant = falses(nv(g)) # descendant in s
    blocked = falses(nv(g))

    for ve in S
        in_seen[ve] = true
        blocked[ve] = true
    end

    (in_seen[u] || in_seen[v]) && throw(ArgumentError("S should not contain u or v"))

    u == v && throw(ArgumentError("u == v"))

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

    push!(in_next, u) # treat u as vertex reached backwards
    in_seen[u] = true

    while true
        sin = isempty(in_next)
        sout = isempty(out_next)
        sin && sout && return true # no vertices in the queue

        if !sin # some vertex reach backwards in the queue
            src = popfirst!(in_next)
            for w in outneighbors(g, src) # possible collider at destination
                if !out_seen[w] && (!blocked[w] || descendant[w])
                    verbose && println("<- $src -> $w")
                    w == v && return false
                    push!(out_next, w)
                    out_seen[w] = true
                end
            end
            for w in inneighbors(g, src)
                if !in_seen[w]
                    verbose && println("<- $src <- $w")
                    w == v && return false
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
                    w == v && return false
                    push!(out_next, w)
                    out_seen[w] = true
                end
            end
            for w in inneighbors(g, src) # collider at source
                if !in_seen[w] && descendant[src] # shielded collider
                    verbose && println("-> $src <- $w")
                    w == v && return false
                    push!(out_next, w)
                    in_seen[w] = true
                end
            end
        end
    end
end
