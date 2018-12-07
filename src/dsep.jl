
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
