import Base: iterate, length

"""
    alt_cpdag(g::DiGraph)

Computes CPDAG from a DAG using a simpler adaption of Chickering's conversion algorithm without ordering all edges
(equivalent to the PC algorithm with `d`-separation as independence test, see `pc_oracle`.)

Reference: M. Chickering: Learning equivalence classes of Bayesian
network structures. Journal of Machine Learning Research 2 (2002).
M. Chickering: A Transformational Characterization of Equivalent Bayesian Network Structures. (1995).
"""
function alt_cpdag(g)
    compelled = SimpleDiGraph(nv(g))
    reversible = SimpleDiGraph(nv(g))
    to = topological_sort_by_dfs(g)
    invto = invperm(to)
    iscompelled = falses(nv(g))
    for y in to
        indegree(g, y) == 0 && continue
        x = argmax(x -> invto[x], inneighbors(g, y))
        iscompelled[inneighbors(g, y)] .= false
        for w in inneighbors(compelled, x)
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
                add_edge!(compelled, v, y) 
            else
                add_edge!(reversible, v, y)
                add_edge!(reversible, y, v)
            end
        end
    end
    return union(compelled, reversible)
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
        ys = topological_sort_by_dfs(g)
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
    ys = topological_sort_by_dfs(g)
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
    ys = topological_sort_by_dfs(g)
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
