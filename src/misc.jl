function digraph(E)
    d = maximum(flatten(E))
    g = DiGraph(d)
    for (i, j) in E
        add_edge!(g, i, j)
    end
    g
end

pairs(g) = map(Pair, collect(edges(g)))

oracle(g) = skeleton(nv(g), dseporacle, g)
pc_oracle(g) = pcalg(nv(g), dseporacle, g)

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

