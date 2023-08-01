using CausalInference
using Graphs
using Test
using Random
using Graphs
using Distributions
using CausalInference: disjoint_sorted

@test disjoint_sorted([], [1, 2])
@test disjoint_sorted([1, 2], [])
@test disjoint_sorted([], [])
@test disjoint_sorted([1, 2], [3, 4])
@test disjoint_sorted([1, 2], [3])
@test disjoint_sorted([1], [3, 4])
@test !disjoint_sorted([1, 3], [2, 3, 4])
@test !disjoint_sorted([1, 2, 3], [3, 4, 5])

# tests modelled after examples discussed in chapter 2 of
# Pearl, Judea. Causality. Cambridge University Press, 2009.

g = DiGraph(5)
d = nv(g)
for (i, j) in [(1, 2), (1, 3), (2, 4), (3, 4), (4, 5)]
    add_edge!(g, i, j)
end

h, s = skeleton(d, dseporacle, g)
Z = orientable_unshielded(h, s)
@testset "unshielded" begin
    @test Graph(g) == h
    for z in Z
        u, v, w = z
        n = inneighbors(h, v)
        @test u in n
        @test w in n
        @test !(u in neighbors(h, w))
    end
end

dg = pcalg(d, dseporacle, g)

Random.seed!(123)
N = 1000
p = 0.01
x = randn(N)
v = x + randn(N) * 0.25
w = x + randn(N) * 0.25
z = v + w + randn(N) * 0.25
s = z + randn(N) * 0.25

X = [x v w z s]
C = cor(X)
df = (x = x, v = v, w = w, z = z, s = s)

println("Timing Gaussian pcalg")
@time gaussci_g = pcalg(df, p, gausscitest)

println("Timing CMI pcalg")
@time cmi_g = pcalg(map(x->getindex(x, 1:100), df), 0.1, cmitest) # use less observations for test performance

@testset "pcalg_edgde_test" begin
    @test collect(Graphs.edges(gaussci_g)) == collect(Graphs.edges(dg))
    @test_broken collect(Graphs.edges(cmi_g)) == collect(Graphs.edges(dg)) # not enough observations
end

@testset "PC alg plotting utils" begin
    g = DiGraph(5)
    for (i, j) in [(1, 2), (2, 3), (2, 4), (4, 2), (4, 5)]
        add_edge!(g, i, j)
    end
    objs = CausalInference.prepare_pc_graph(g)
    @test keys(objs) == (:plot_g, :node_labels, :edge_styles, :node_style, :options)
    @test nv(objs.plot_g) == nv(g)
    @test ne(objs.plot_g) == (ne(g) - 1) # bi-directional edge removed
    @test objs.edge_styles[(1, 2)] == "->" # default 
    @test objs.edge_styles[(2, 4)] == "--" # recognize bi-directional edge
    @test objs.node_labels == string.(1:nv(g))

    objs = CausalInference.prepare_pc_graph(g, collect(string.('a':'e')))
    @test keys(objs) == (:plot_g, :node_labels, :edge_styles, :node_style, :options)
    @test nv(objs.plot_g) == nv(g)
    @test ne(objs.plot_g) == (ne(g) - 1) # bi-directional edge removed
    @test objs.edge_styles[(1, 2)] == "->" # default 
    @test objs.edge_styles[(2, 4)] == "--" # recognize bi-directional edge
    @test objs.node_labels == collect(string.('a':'e'))
end
