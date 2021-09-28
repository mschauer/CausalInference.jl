using CausalInference
using LightGraphs
using Test
using Random
using LightGraphs
using Distributions
using CategoricalArrays
using CausalInference: disjoint_sorted

@test disjoint_sorted([],[1,2])
@test disjoint_sorted([1,2],[])
@test disjoint_sorted([],[])
@test disjoint_sorted([1,2],[3,4])
@test disjoint_sorted([1,2],[3])
@test disjoint_sorted([1],[3,4])
@test !disjoint_sorted([1,3],[2,3,4])
@test !disjoint_sorted([1,2,3],[3,4,5])

# tests modelled after examples discussed in chapter 2 of
# Pearl, Judea. Causality. Cambridge University Press, 2009.

g = DiGraph(5)
d = nv(g)
for (i,j) in [(1,2), (1,3), (2,4), (3,4), (4,5)]
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
v = x + randn(N)*0.25
w = x + randn(N)*0.25
z = v + w + randn(N)*0.25
s = z + randn(N)*0.25

X = [x v w z s]
C = cor(X)
df = (x=x, v=v, w=w, z=z, s=s)

println("Running Gaussian tests")
@time gaussci_g = pcalg(df, p, gausscitest)

println("Running CMI tests")
@time cmi_g = pcalg(df, 0.1, cmitest)

@testset "continuous pcalg_edge_test" begin
    @test collect(LightGraphs.edges(cmi_g)) == collect(LightGraphs.edges(dg))
    @test collect(LightGraphs.edges(gaussci_g)) == collect(LightGraphs.edges(dg))
end

x = rand(1:5, N)
v = x + rand(1:2, N)
w = x + rand(1:2, N)
z = ((v + w)/2) + rand(1:2, N)
s = z + rand(1:2, N)

df = (x=CategoricalArray(x), v=CategoricalArray(v),
      w=CategoricalArray(w), z=CategoricalArray(z),
      s=CategoricalArray(s))

@time cmi_cat_g = pcalg(df, 0.1, cmitest)
@testset "categorical pcalg_edge_test" begin
    @test collect(LightGraphs.edges(cmi_cat_g)) == collect(LightGraphs.edges(dg))
end
