export has_a_path


"""
    has_a_path(g::AbstractGraph, U::Vector, V::Vector, exclude_vertices::AbstractVector = T[], nbs=Graphs.outneighbors)

Find if there is a path connecting U with V not passing `exclude_vertices`, where ` nbs=Graphs.outneighbors`
determines the direction of traversal. 
"""
function has_a_path(g::AbstractGraph{T}, U::Vector, V::Vector, exclude_vertices::AbstractVector = T[], nbs=Graphs.outneighbors) where {T}
    seen = zeros(Bool, nv(g))
    target = zeros(Bool, nv(g))
    for ve in V
        target[ve] = true
    end
    for ve in exclude_vertices # mark excluded vertices as seen
        seen[ve] = true
        target[ve] = false # excluded target
    end
    any(target) || return false # no vertex to go to
    next = Vector{T}()
    for u in U
        target[u] && return true
        push!(next, u)
        seen[u] = true
    end
    while !isempty(next)
        src = popfirst!(next) # get new element from queue
        for vertex in outneighbors(g, src)
            if !seen[vertex]
                target[vertex] && return true
                push!(next, vertex) # push onto queue
                seen[vertex] = true
            end
        end
    end
    return false
end 