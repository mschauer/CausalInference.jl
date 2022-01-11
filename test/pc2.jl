using Distributions, CausalInference, Graphs
using Test, Random
# Example from Pearl, Judea. Causality. Cambridge University Press, 2009.

Random.seed!(100)
N = 2000
d = Normal()
p = 0.01

x = randn(N)
v = x + randn(N)*0.5
w = x + randn(N)*0.5
z = v + w + randn(N)*0.5
s = z + randn(N)*0.5
X = [x v w z s]
C = cor(X)
g = pcalg(size(X,2), gausscitest, (C, N), quantile(Normal(), 1-p/2))

@test collect(Graphs.edges(g)) == map(x->Edge(x...), [1 => 2
 1 => 3
 2 => 1
 2 => 4
 3 => 1
 3 => 4
 4 => 5])
