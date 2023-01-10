using TabularDisplay

export has_a_path, graph_to_text

"""
    has_a_path(g::AbstractGraph, U::Vector, V::Vector, exclude_vertices::AbstractVector = T[], nbs=Graphs.outneighbors)

Find if there is a path connecting U with V not passing `exclude_vertices`, where ` nbs=Graphs.outneighbors`
determines the direction of traversal. 
"""
function has_a_path(g::AbstractGraph{T}, U::Vector, V::Vector,
                    exclude_vertices::AbstractVector = T[],
                    nbs = Graphs.outneighbors) where {T}
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

"""
    graph_to_text(g::AbstractGraph, node_labels::AbstractVector{<:AbstractString}=[]; edge_styles::AbstractDict=Dict())

Print out a graph `g` as a table of edges labeled by `node_labels`.

# Arguments
- `g::AbstractGraph`: a graph to print
- `node_labels::AbstractVector{<:AbstractString}=[]`: labels for nodes (same order as indices of nodes in `g`)
- `edge_styles::AbstractDict=Dict()`: dictionary of edge styles (e.g. `Dict((1, 2) => "->", (2, 3) => "<->")`)

# Example
```
g = DiGraph(4)
for (i, j) in [(1, 2), (2, 3), (2, 4)]
    add_edge!(g, i, j)
end
graph_to_text(g)
```
"""
function graph_to_text end;
function graph_to_text(g::AbstractGraph,
                       node_labels::AbstractVector{<:AbstractString} = String[];
                       edge_styles::AbstractDict = Dict())
    graph_to_text(stdout, g, node_labels, edge_styles = edge_styles)
end
function graph_to_text(io::IO, g::AbstractGraph,
                       node_labels::AbstractVector{<:AbstractString} = String[];
                       edge_styles::AbstractDict = Dict())
    if length(node_labels) == 0
        node_labels = map(string, 1:nv(g))
    end
    @assert length(node_labels)==nv(g) "node_labels must be the same length as number of nodes (provided: $(length(node_labels)), expected: $(nv(g)))"

    edges_array = String[]
    for e in edges(g)
        arrow_symbol = get(edge_styles, (e.src, e.dst), "->")
        push!(edges_array, string(node_labels[e.src], arrow_symbol, node_labels[e.dst]))
    end
    displaytable(io, edges_array)
end
