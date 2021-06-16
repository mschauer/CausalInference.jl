export has_recanting_witness
function has_recanting_witness(g::AbstractGraph{T}, u::Integer, v::Integer,  blocked_edges::AbstractGraph{T}) where T
    on_blocked = zeros(Bool, nv(g))
    on_open = zeros(Bool, nv(g))
    may_blocked = zeros(Bool, nv(g))
    for e in edges(blocked_edges) # mark excluded vertices as seen
        may_blocked[src(e)] = true
    end
    u == v && error("u == v")

    next_open = Vector{T}()
    next_blocked = Vector{T}()
    
    push!(next_open, u)
    on_open[u] = true
    while !isempty(next_open)
        src = popfirst!(next_open) # get new element from queue
        for vertex in outneighbors(g, src)
            if !on_open[vertex]
                if may_blocked[src] && has_edge(blocked_edges, src, vertex)
                    on_blocked[vertex] = true
                    if vertex != v 
                        push!(next_blocked, vertex)  # push onto blocked queue
                    end
                else
                    on_open[vertex] = true
                    if vertex != v 
                        push!(next_open, vertex) # push onto queue
                    end
                end
            end
        end
    end
    while !isempty(next_blocked)
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
    @show on_blocked
    @show on_open
    if on_blocked[v] && on_open[v]
        return true
    end
    return false
end