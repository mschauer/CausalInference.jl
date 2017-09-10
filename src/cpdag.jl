import Base: start, next, done, length

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

type OrderedEdges
    g
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
start(iter::OrderedEdges) = 0, 0, 0, Int[] #i, j
function done(iter::OrderedEdges, state) 
    i, j, k, xs  = state
    k >= iter.m
end

function next(iter::OrderedEdges, state)
    i, j, k, xs  = state
    g, ys, yp = iter.g, iter.ys, iter.yp
    while i == 0
        j = j + 1
        xs = sort(in_neighbors(g, ys[j]), by=x->yp[x])
        i = length(xs)
    end
    Edge(xs[i], ys[j]), (i - 1, j, k + 1, xs)
end

# came in handy while prototyping
function print_chickering_order(g)
    ys = topological_sort_by_dfs(g)
    yp = sortperm(ys)
    println(ys)
    n = nv(g)
    i = 0 # src index
    j = 0 # dst index
    k = 1 # edge number

    while true
        while i == 0
            j = j + 1
            if j > n
                @goto ende 
            end
            xs = sort(in_neighbors(g, ys[j]), by=x->yp[x])
            i = length(xs)
        end

        x, y = xs[i], ys[j]
        println(k, " ", x=>y)
        k += 1
        i -= 1
    end
    @label ende
end


@enum Status reversible=0 compelled=1 unknown=2
"""
    cpdag(skel::DiGraph)

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
            xs = sort(in_neighbors(g, ys[j]), by=x->yp[x])
            i = length(xs)
        end

        x, y = xs[i], ys[j]

        if label[x => y] == unknown
            for w in in_neighbors(g, x)
                label[w => x] == compelled || continue
                if !has_edge(g, w => y)
                    for z in in_neighbors(g, y)
                        label[z => y] = compelled
                    end
                    @goto next
                else
                    label[w => y] = compelled
                end
            end

            lbl = reversible
            for z in in_neighbors(g, y)
                if x == z || has_edge(g, z => x) 
                    continue
                end
                lbl = compelled
            end
            for z in in_neighbors(g, y)
                if label[z => y] == unknown
                    label[z => y] = lbl
                end
            end
        end
        @label next

        i = i - 1
    end
    @label ende
    for (x, y) in collect(edges(g))
        if label[x=>y] == reversible
            add_edge!(g, y => x)
        end
    end
    g
end