import Base: iterate, length

# REMARKS:
# - implemented own topological sort temporarily as topological_sort_by_dfs from Julia Graphs appears to have an issue (quadratic run-time for sparse graphs)
# - cpdag(g) potentially has O(m * sqrt(m)) worst-case run-time due to iterating over all w -> x, not only compelled ones. However, in contrast to topological_sort_by_dfs, quadratic run-time behaviour for sparse graphs (e.g. O(n*n*)/O(m*m)) does not seem to appear.

function dfs(g, u, visited, to)
    visited[u] = true
    for v in outneighbors(g, u)
        !visited[v] && dfs(g, v, visited, to)
    end
    push!(to, u)
end

function tsort(g)
    visited = falses(nv(g))
    to = Vector{Int64}()
    for u in vertices(g)
        !visited[u] && dfs(g, u, visited, to)
    end
    return reverse(to)
end

function nvertexdigraphfromedgelist(n, edgelist)
    g = SimpleDiGraph(edgelist)
    while nv(g) < n
        add_vertex!(g)
    end
    return g
end

"""
    alt_cpdag(g::DiGraph)

Computes CPDAG from a DAG using a simpler adaption of Chickering's conversion algorithm without ordering all edges
(equivalent to the PC algorithm with `d`-separation as independence test, see `pc_oracle`.)

Reference: M. Chickering: Learning equivalence classes of Bayesian
network structures. Journal of Machine Learning Research 2 (2002).
M. Chickering: A Transformational Characterization of Equivalent Bayesian Network Structures. (1995).
"""
function alt_cpdag(g)
    compelledingoing = [Vector{Int64}() for _ = 1:nv(g)]
    edgelist = Vector{Edge{Int64}}()
    to = tsort(g)
    invto = invperm(to)
    iscompelled = falses(nv(g))
    for y in to
        indegree(g, y) == 0 && continue
        x, xidx = -1, -1
        for z in inneighbors(g, y)
            invto[z] > xidx && ((x, xidx) = (z, invto[z]))
        end
        iscompelled[inneighbors(g, y)] .= false
        for w in compelledingoing[x]
            if !has_edge(g, w, y)
                iscompelled[inneighbors(g, y)] .= true
                @goto add_edges
            else
                iscompelled[w] = true
            end
        end
        for z in inneighbors(g, y)
            z == x && continue
            if !has_edge(g, z, x)
                iscompelled[inneighbors(g, y)] .= true
                @goto add_edges
            end
        end
        @label add_edges
        for v in inneighbors(g, y)
            if iscompelled[v]
                push!(compelledingoing[y], v)
                push!(edgelist, Edge(v, y))
            else 
                push!(edgelist, Edge(v, y))
                push!(edgelist, Edge(y, v))
            end
        end
    end
    return nvertexdigraphfromedgelist(nv(g), edgelist)
end

# TODOs:
# - add function for checking whether a graph is a CDPAG
# - add function for checking whether a PDAG has a DAG extension
function hasdirectedcycle(G)
    # TODO: implement
    return true
end

function ismaximallyoriented(G)
    # TODO: implement
    return true
end

function checkundirectedcomps(G)
    # TODO: implement
    return true, true
end

function isstronglyprotected(G)
    # TODO: implement
    return true
end

"""
    classifygraph(G)
Classifies a graph efficiently into the following classes:
- "cyclic": has a directed cycle (the graph can also satisfy some of the other criteria)
- "not maximally oriented": there is an undirected edge whose direction would follow from one of the Meek rules (Meek 1995)
- "not extendable": there is no consistent DAG extension of this graph
- "CPDAG": the graph is a CPDAG
- "MPDAG, counting works" and "MPDAG, counting not yet implemented": in both cases the graph is an MPDAG, however, in the second case it has not enough structure for the counting algorithm to work (TODO: add more details)
"""

function classifygraph(G)
    hasdirectedcycle(G) && return "cyclic"
    ismaximallyoriented(G) && return "not maximally oriented"
    areinduced, arechordal = checkundirectedcomps(G)
    !arechordal && return "not extendable"
    if areinduced
        if isstronglyprotected(G)
            return "CPDAG"
        else
            return "MPDAG, counting works"
        end
    else
        return "MPDAG, counting not yet implemented"
    end
end

function consistent_extension(G)
    # TODO: implement Dor-Tarsi algorithm
end


"""
    ordered_edges(dag)

Iterator of edges of a dag, ordered in Chickering order:

    Perform a topological sort on the NODES
    while there are unordered EDGES in g
        Let y be the lowest ordered NODE that has an unordered EDGE incident into it
        Let x be the highest ordered NODE for which x => y is not ordered
        return x => y 
    end
"""
ordered_edges(g) = OrderedEdges(g)

struct OrderedEdges
    g::Any
    ys::Vector{Int}
    yp::Vector{Int}
    n::Int
    m::Int

    function OrderedEdges(g::DiGraph)
        ys = tsort(g)
        yp = sortperm(ys)
        new(g, ys, yp, nv(g), ne(g))
    end
end

length(iter::OrderedEdges) = iter.m

function iterate(iter::OrderedEdges, state = (0, 0, 0, Int[]))
    i, j, k, xs = state
    k >= iter.m && return nothing
    g, ys, yp = iter.g, iter.ys, iter.yp
    while i == 0
        j = j + 1
        xs = sort(inneighbors(g, ys[j]), by = x -> yp[x])
        i = length(xs)
    end
    Edge(xs[i], ys[j]), (i - 1, j, k + 1, xs)
end

# came in handy while prototyping
function chickering_order(g)
    outp = Pair{Int, Int}[]
    ys = tsort(g)
    yp = sortperm(ys)
    n = nv(g)
    i = 0 # src index
    j = 0 # dst index
    xs = Int[]
    while true
        while i == 0
            j = j + 1
            if j > n
                @goto ende
            end
            xs = sort(inneighbors(g, ys[j]), by = x -> yp[x])
            i = length(xs)
        end
        x, y = xs[i], ys[j]
        push!(outp, x => y)
        i -= 1
    end
    @label ende
    outp
end

@enum Status reversible=0 compelled=1 unknown=2
"""
    cpdag(skel::DiGraph)

Computes CPDAG from a DAG using Chickering's conversion algorithm 
(equivalent to the PC algorithm with `d`-seperation as independence test, see `pc_oracle`.)

Reference: M. Chickering: Learning equivalence classes of Bayesian
network structures. Journal of Machine Learning Research 2 (2002).
M. Chickering: A Transformational Characterization of Equivalent Bayesian Network Structures. (1995).


Note that the edge order defined there is already partly encoded into the
representation of a DiGraph.
"""
function cpdag(skel::DiGraph)
    g = copy(skel)
    ys = tsort(g)
    yp = sortperm(ys)
    n = nv(g)
    m = ne(g)
    label = Dict(zip(map(Pair, edges(g)), repeated(unknown)))

    i = 0 # src index
    j = 0 # dst index
    xs = Int[]
    while true
        while i == 0
            j = j + 1
            if j > n
                @goto ende
            end
            xs = sort(inneighbors(g, ys[j]), by = x -> yp[x])
            i = length(xs)
        end

        x, y = xs[i], ys[j]

        if label[x => y] == unknown
            for w in inneighbors(g, x)
                label[w => x] == compelled || continue
                if !has_edge(g, w => y)
                    for z in inneighbors(g, y)
                        label[z => y] = compelled
                    end
                    @goto next
                else
                    label[w => y] = compelled
                end
            end

            lbl = reversible
            for z in inneighbors(g, y)
                if x == z || has_edge(g, z => x)
                    continue
                end
                lbl = compelled
            end
            for z in inneighbors(g, y)
                if label[z => y] == unknown
                    label[z => y] = lbl
                end
            end
        end
        @label next

        i = i - 1
    end
    @label ende
    for e in collect(edges(g))
        x, y = Tuple(e)
        if label[x => y] == reversible
            add_edge!(g, y => x)
        end
    end
    g
end
"""
    vskel(g)

Reduce a (P)DAG to skeleton and `v`-structures. 
"""
function vskel(g::DiGraph)
    g2 = DiGraph(Graph(g))
    for v in vertices(g)
        ns = parents(g, v) 
        n = length(ns)
        protected = falses(n) # mark parents which are v structures
        for (i, j) in combinations(1:n, 2) 
            if !isadjacent(g, ns[i], ns[j]) 
                protected[i] = protected[j] = true
            end
        end
        for i in 1:n
            if protected[i]
                remove!(g2, v → ns[i])
            end
        end
    end
    g2
end
function vskel!(g::DiGraph)
    es = Pair{Int,Int}[]
    for v in vertices(g)
        ns = parents(g, v)
        n = length(ns)
        protected = falses(n) # mark parents which are v structures
        for (i, j) in combinations(1:n, 2) 
            if !isadjacent(g, ns[i], ns[j]) 
                protected[i] = protected[j] = true
            end
        end
        for i in 1:n
            if !protected[i]
                push!(es, v=>ns[i]) # make undirected
            end
        end
    end
    for e in es
        add_edge!(g, e)
    end
    g
end