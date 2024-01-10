using Random, CausalInference, Statistics, Test, Graphs
Random.seed!(1)

# Generate some sample data to use with the GES algorithm

N = 2000 # number of data points

# define simple linear model with added noise
x = randn(N)
v = x + randn(N)*0.25
w = x + randn(N)*0.25
z = v + w + randn(N)*0.25
s = z + randn(N)*0.25

df = (x=x, v=v, w=w, z=z, s=s)
iterations = 5_000
n = length(df) # vertices
κ = n - 1 # max degree
penalty = 2.0 # increase to get more edges in truth
Random.seed!(101)
C = cor(CausalInference.Tables.matrix(df))
score = GaussianScore(C, N, penalty)
gs = @time causalzigzag(n; score, κ, iterations)
graphs, graph_pairs, hs, τs, ws, ts, scores = CausalInference.unzipgs(gs)
posterior = sort(keyedreduce(+, graph_pairs, ws); byvalue=true, rev=true)

@test first(posterior).first == [1=>2, 1=>3, 2=>1, 2=>4, 3=>1, 3=>4, 4=>5] 
@test score_dag(pdag2dag!(copy(graphs[end])), score) ≈ scores[end] + score_dag(DiGraph(n), score)