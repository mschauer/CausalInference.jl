using CausalInference
using Graphs
using Test
using Random
using Statistics
using TikzGraphs
using TikzPictures # not in dependencies!

# deprecated in favour of explicit Tikz-backed function `plot_pc_graph_tikz`
# include("plotdag.jl") # required for `plot_dag`

# After examples discussed in chapter 2 of
# Pearl, Judea. Causality. Cambridge University Press, 2009.

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

println("Running Gaussian tests")
@time est_g = pcalg(df, p, gausscitest)

variables = [String(k) for k in keys(df)]
tp = plot_pc_graph_tikz(est_g, variables)
save(PDF("estdag"), tp)

# This should be the true DAG

g = DiGraph(5)
d = nv(g)
for (i, j) in [(1, 2), (1, 3), (2, 4), (3, 4), (4, 5)]
    add_edge!(g, i, j)
end

tp = plot_pc_graph_tikz(g, variables)
save(PDF("truedag"), tp)

# Use the oracle to tell about a graph what can be learned from observations
dg = pcalg(d, dseporacle, g)

tp = plot_pc_graph_tikz(dg, variables)
save(PDF("equivalencedag"), tp)

@testset "pctest" begin @test collect(Graphs.edges(est_g)) == collect(Graphs.edges(dg)) end
