# Edges in partially directed graphs

"""
    isadjacent(g, edge::Edge)
    isadjacent(g, x, y)

Test if `x` and `y` are connected by a any edge in the graph `g`.
Can also perform the same test given an `edge`.
"""
isadjacent(g, edge) = has_edge(g, edge) || has_edge(g, reverse(edge))
isadjacent(dg, v, w) = has_edge(dg, v, w) || has_edge(dg, w, v)

remove!(dg::DiGraph, e::Pair) = rem_edge!(dg, Edge(e))
remove!(dg::Graph, e::Tuple) = rem_edge!(dg, Edge(e))


"""
    isundirected(g, edge::Edge)
    isundirected(g, x, y)

Test if `x` and `y` are connected by a undirected edge in the graph `g`.
Can also perform the same test given an `edge`.
"""
isundirected(g, x, y) = has_edge(g, x, y) && has_edge(g, y, x)
isundirected(g, edge) = has_edge(g, edge) && has_edge(g, reverse(edge))
has_both(dg, v, w) = has_edge(dg, v, w) && has_edge(dg, w, v)


"""
    isparent(g, edge::Edge)
    isparent(g, x, y)

Test if `x` is a parent of `y` in the graph `g`, meaning x→y.
Can also perform the same test given an `edge`.
"""
isparent(g, edge) = has_edge(g, edge) && !has_edge(g, reverse(edge))
isparent(g, x, y) = has_edge(g, x, y) && !has_edge(g, y, x)


"""
    ischild(g, edge::Edge)
    ischild(g, x, y)

Test if `x` is a child of `y` in the graph `g`, meaning x←y.
Can also perform the same test given an `edge`.
"""
ischild(g, edge) = !has_edge(g, edge) && has_edge(g, reverse(edge))
ischild(g, x, y) = !has_edge(g, x, y) && has_edge(g, y, x)


"""
    isdescendent(g, edge::Edge)
    isdescendent(g, x, y)

Return `true` if `x`←`y` OR `x`-`y` in the graph `g`.
"""
isdescendent(g, edge) = has_edge(g, reverse(edge))
isdescendent(g, x, y) = has_edge(g, y, x)

"""
    isoriented(g, edge::Edge)
    isoriented(g, x, y)

Test if `x` and `y` are connected by a directed edge in the graph `g`, either x←y OR x→y.
Can also perform the same test given an `edge`.
"""
isoriented(g, edge) = has_edge(g,edge) ⊻ has_edge(g,reverse(edge)) # xor
isoriented(g, x, y) = has_edge(g,x,y) ⊻ has_edge(g,y,x)


"""
    isclique(g, nodes)

Return `true` if all vertices in `nodes` are connected to each other in the graph `g`.
"""
function isclique(g, nodes)

    for nodePair in allpairs(nodes)
        if !isadjacent(g, nodePair)
            return false
        end
    end

    return true
end




"""
    orientedge!(g, x, y)

Update the edge `x`-`y` to `x`→`y` in the graph `g`. 
"""
orientedge!(g, edge) = rem_edge!(g, reverse(edge))
orientedge!(g, x, y) = rem_edge!(g, y, x)

"""
    neighbors_undirected(g, x)

All vertices in `g` connected to `x` by an undirected edge.
"""
neighbors_undirected(g, x) = intersect(inneighbors(g, x), outneighbors(g, x))
#neighbors_undirect(g, x) = neighborsGeneral(isundirected, g, x)
# TODO: check if neighbors are sorted.


#these are just aliases for the functions above
parents(g, x) = inneighbors(g, x)
children(g, x) = outneighbors(g, x)
