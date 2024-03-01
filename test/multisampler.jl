using Random, CausalInference, Statistics, Test, Graphs, LinearAlgebra
@testset "MultiSampler" begin
    Random.seed!(1)

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
    Random.seed!(101)
    C = cor(CausalInference.Tables.matrix(df))
    score = GaussianScore(C, N, penalty)
    global bestgraph, samplers = multisampler(n; score, iterations)
    #posterior = sort(keyedreduce(+, graph_pairs, ws); byvalue=true, rev=true)

    # maximum aposteriori estimate
    MAP = [1=>2, 1=>3, 2=>1, 2=>4, 3=>1, 3=>4, 4=>5]
    @test bestgraph == digraph(MAP, n)
    @test abs(-(extrema(getfield.(last.(samplers), :τ))...)) < 0.1
    cm = sort(countmap(vpairs.(getfield.(last.(samplers), :g))), byvalue=true, rev=true)
    @test first(cm).first == MAP
end #testset