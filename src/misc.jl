using Random
"""
    digraph(E)

Create `DiGraph` from edge-list.
"""
function digraph(E)
    d = maximum(flatten(E))
    g = DiGraph(d)
    for (i, j) in E
        add_edge!(g, i, j)
    end
    g
end

"""
    graph(E)

Create `Graph` from edge-list.
"""
function graph(E)
    d = maximum(flatten(E))
    g = Graph(d)
    for (i, j) in E
        add_edge!(g, i, j)
    end
    g
end


"""
    vpairs(g)

Return the edge-list as `Pair`s.
"""
vpairs(g) = map(Pair, collect(edges(g)))

"""
    skel_oracle(g)

Compute the `skeleton` using the `dseporacle` for the DAG `g`.
"""
skel_oracle(g) = skeleton(nv(g), dseporacle, g)


"""
    pc_oracle(g)

Compute CPDAG using the PC algorithm using the `dseporacle` on the DAG `g`. 
"""
pc_oracle(g) = pcalg(nv(g), dseporacle, g)

"""
    randdag(n, alpha = 0.1)

Create random DAG from randomly permuted random triangular matrix with
edge probability `alpha`.
"""
function randdag(n, alpha = 0.1)
    g = DiGraph(n)
    p = randperm(n)
    for i in 1:n
        for j in i+1:n
            rand() < alpha && add_edge!(g, p[i], p[j])
        end
    end
    g
end

