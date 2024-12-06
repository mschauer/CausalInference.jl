
""""
    do!(g, v)

Graphical do operator, removes all 
incoming edges to vertex `v` in a DiGraph `g`.
Returns the modified graph.
"""
function do!(g::DiGraph{T}, v::T) where {T}
    for u in collect(inneighbors(g, v))
        rem_edge!(g, u, v)
    end
    return g
end
