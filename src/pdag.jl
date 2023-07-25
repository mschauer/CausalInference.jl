# Partially directed graphs

#= Terminology and representation

PDAGs are implemented as DAG with undirected edges x --- y represented 
as the pair x --> y and y --> x.

neighbors_undirected(g, x) - Neighbors of x in g are vertices y such that there is an undirected edge y --- x,
thus different from `neighbors(g, x)` from Graphs.

Unshielded collider: triple of vertices x, y, z with directed edges
   x --> y <-- z, and x and z are not adjacent

Semi-directed path: A path of directed and undirected edges, where no
directed edge is going in the opposite direction. `outneighbors(g, x)` gives 
possible continuations of semi-directed path. A directed path has only directed edges.
`children(g, x)` gives possible continuations.
=#

x → y = Edge(x, y)
# x => y # Pair as directed edge
struct Undirected{T<:Integer}
    x::T
    y::T
Undirected(x::T, y::T) where {T} = new{T}(minmax(x, y)...)
end
struct Directed{T<:Integer}
    x::T
    y::T
Directed(x::T, y::T) where {T} = new{T}(x, y)
end

x ⟹ y = Directed(x, y)
x ⟺ y = Undirected(x, y)

isolated(g, v) = isempty(inneighbors(g, v)) && isempty(outneighbors(g, v))

"""
    isadjacent(g, x, y)

Test if `x` and `y` are connected by a any edge in the graph `g`
(i.e. x --- y, x --> y, or x <-- y .)
"""
isadjacent(dg, v, w) = has_edge(dg, v, w) || has_edge(dg, w, v)

"""
    remove!(dg::DiGraph, x => y)

Removes directed edge.
"""
remove!(dg::DiGraph, e::Edge) = rem_edge!(dg, e)
remove!(dg::Graph, e::Tuple) = rem_edge!(dg, Edge(e))

"""
    isundirected(g, edge::Edge)
    isundirected(g, x, y)

Test if `x` and `y` are connected by a undirected edge in the graph `g`.
"""
isundirected(g, x, y) = has_edge(g, x, y) && has_edge(g, y, x)
isundirected(g, edge) = has_edge(g, edge) && has_edge(g, reverse(edge))
const has_both = isundirected

"""
    isparent(g, x, y)

Test if `x` is a parent of `y` in the graph `g`, meaning x→y.
"""
isparent(g, x, y) = has_edge(g, x, y) && !has_edge(g, y, x)


"""
    ischild(g, x, y)

Test if `x` is a child of `y` in the graph `g`, meaning x←y.
"""
ischild(g, x, y) = !has_edge(g, x, y) && has_edge(g, y, x)

"""
    isoriented(g, edge::Edge)
    isoriented(g, x, y)

Test if `x` and `y` are connected by a directed edge in the graph `g`, either x←y OR x→y.
Can also perform the same test given an `edge`.
"""
isoriented(g, edge) = has_edge(g,edge) ⊻ has_edge(g, reverse(edge)) # xor
isoriented(g, x, y) = has_edge(g,x,y) ⊻ has_edge(g,y,x)


"""
    isclique(g, nodes)

Return `true` if all vertices in `nodes` are adjacent to each other in the graph `g`.
Chickering, "Learning equivalence classes" (2002).
Note that a 3-clique with already two undirected edges is also a clique in the neighbor sense.
"""
function isclique(g, nodes)
    for (u, v) in allpairs(nodes)
        if !isadjacent(g, u, v)
            return false
        end
    end
    return true
end


"""
    orientedge!(g, x, y)

Update the edge `x`-`y` to `x`→`y` in the graph `g`. 
Throws error if not `x - y`
"""
orientedge!(g, edge) = @assert has_edge(g, edge) && rem_edge!(g, reverse(edge))
orientedge!(g, x, y) = @assert has_edge(g, x, y) && rem_edge!(g, y, x)

"""
    neighbors_undirected(g, x)

All vertices in `g` connected to `x` by an undirected edge.
Returns sorted array.
"""
neighbors_undirected(g, x) = inneighbors(g, x) ∩ outneighbors(g, x)

"""
    neighbors_adjacent(g, x)

All vertices in `g` connected to `x` by an any edge.
Returns sorted array.
"""
neighbors_adjacent(g, x) = sort(outneighbors(g, x) ∪ inneighbors(g, x))

"""
    parents(g, x)

Parents are vertices y such that there is a directed edge y --> x.
Returns sorted array.
"""
parents(g, x) = setdiff(inneighbors(g, x), outneighbors(g, x))

"""
    children(g, x)

Children of x in g are vertices y such that there is a directed edge y <-- x.
Returns sorted array.
"""
children(g, x) = setdiff(outneighbors(g, x), inneighbors(g, x))
