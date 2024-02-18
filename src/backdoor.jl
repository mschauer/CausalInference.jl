export backdoor_criterion
export backdoors

# TODOs
# - start with implementing a general Bayes-Ball type graph search
# - then implement efficient algorithms for finding backdoor adjustment sets


function backdoors(g::AbstractGraph{T}, U::Vector) where {T}
    seen = zeros(Bool, nv(g))
    bd = zeros(Bool, nv(g))
    next = Vector{T}()
    for u in U
        seen[u] = true
    end
    for u in U
        for vertex in outneighbors(g, u)
            if !seen[vertex]
                bd[vertex] = true
                seen[vertex] = true
            end
        end
    end
    findall(bd)
end
"""
    backdoor_criterion(g::AbstractGraph, u::Integer, v::Integer, Z; verbose = false)

Test that given directed graph `g`, no node in `Z` is descendant of `u` and `Z` d-separates `u` from `v` 
in the subgraph that only has the backdoors of `u` left (outgoing edges of `u` removed).

If so, the causal effect of `u` on `v` is identifiable and is given by the formula:

    ∑{z∈Z} p(v | u, z)p(z)

In the linear Gaussian model, find `E[Y | X = x, Z = z] = α + βx + γ'z` and obtain `E[Y | do(x)] = α + βx + γ'E[Z].
"""
function backdoor_criterion(g::AbstractGraph{T}, u::Integer, v::Integer, S = T[];
                            verbose = false) where {T}
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
        w == u && return false # u has descendant in S
    end

    while !isempty(next)
        for w in inneighbors(g, popfirst!(next))
            if !descendant[w]
                push!(next, w) # push onto queue
                descendant[w] = true
                w == u && return false # u has descendant in S
            end
        end
    end

    verbose && println(descendant)

    in_next = Vector{T}()
    out_next = Vector{T}()

    push!(in_next, u) # reached u backwards
    in_seen[u] = true
    out_seen[u] = true

    while true
        sin = isempty(in_next)
        sout = isempty(out_next)
        sin && sout && return true # no vertices in the queue

        if !sin # some vertex reach backwards in the queue
            src = popfirst!(in_next)
            src != u && for w in outneighbors(g, src) # possible collider at destination
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
