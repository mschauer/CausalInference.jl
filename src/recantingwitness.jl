export has_recanting_witness

"""
    has_recanting_witness(g::AbstractGraph, u, v,  blocked_edges::AbstractGraph) -> Bool

In a causal DAG with edges `g`, determine whether path-specific causal effect from vertex `u` to `v`
with edges in `blocked_edges` blocked can be can be computed uniquely from the data available to the investigator (assuming complete observations),
which is the case if there is no "recanting witness".
Essentially this means that `blocked_edges` could equivalently be replaced by a blocking only
outgoing edges of `u`. 

See Alvin, Shpitser, Pearl (2005): "Identifiability of Path-Specific Effects", 
https://ftp.cs.ucla.edu/pub/stat_ser/r321-L.pdf.    
"""
function has_recanting_witness(g::AbstractGraph{T}, u::Integer, v::Integer,
                               blocked_edges::AbstractGraph{T}) where {T}
    on_blocked = zeros(Bool, nv(g))
    on_open = zeros(Bool, nv(g))
    may_blocked = zeros(Bool, nv(g))
    for e in edges(blocked_edges) # mark sources of blocked edges
        may_blocked[src(e)] = true
    end
    u == v && error("u == v")

    for nbr in outneighbors(g, u) # if there is a recanting witness, a neighbour of u is one too
        on_blocked .= false
        on_open .= false
        may_blocked[v] && has_edge(blocked_edges, v, nbr) && continue
        nbr == v && continue # direct edge

        next_open = Vector{T}()
        next_blocked = Vector{T}()

        push!(next_open, nbr)
        on_open[nbr] = true
        while !isempty(next_open) # find all reachable via direct neighbour
            src = popfirst!(next_open) # get new element from queue
            for vertex in outneighbors(g, src)
                if may_blocked[src] && has_edge(blocked_edges, src, vertex)
                    if !on_blocked[vertex]
                        on_blocked[vertex] = true
                        if vertex != v
                            push!(next_blocked, vertex)  # push onto blocked queue
                        end
                    end
                elseif !on_open[vertex]
                    on_open[vertex] = true
                    if vertex != v
                        push!(next_open, vertex) # push onto queue
                    end
                end
            end
        end
        while !isempty(next_blocked) # find all reachable through the blocked edge
            src = popfirst!(next_blocked) # get new element from queue
            for vertex in outneighbors(g, src)
                if !on_blocked[vertex]
                    on_blocked[vertex] = true
                    if vertex != v
                        push!(next_blocked, vertex)  # push onto blocked queue
                    end
                end
            end
        end
        if on_blocked[v] && on_open[v]
            return true
        end
    end
    return false
end
