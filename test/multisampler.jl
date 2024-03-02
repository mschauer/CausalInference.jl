using Random, CausalInference, StatsBase, Statistics, Test, Graphs, LinearAlgebra
@testset "MultiSampler" begin
    Random.seed!(1)

    N = 500 # number of data points

    # define simple linear model with added noise
    x = randn(N)
    v = x + randn(N)*0.25
    w = x + randn(N)*0.25
    z = v + w + randn(N)*0.25
    s = z + randn(N)*0.25
    
    df = (x=x, v=v, w=w, z=z, s=s)
    iterations = 5_000
    penalty = 2.0 # increase to get more edges in truth
    n = length(df) # vertices
    Random.seed!(101)
    C = cor(CausalInference.Tables.matrix(df))
    score = GaussianScore(C, N, penalty)
    bestgraph, samplers = multisampler(n; score, σ=2.0, iterations)
    #posterior = sort(keyedreduce(+, graph_pairs, ws); byvalue=true, rev=true)

    # maximum aposteriori estimate
    MAP = [1=>2, 1=>3, 2=>1, 2=>4, 3=>1, 3=>4, 4=>5]
    @test bestgraph == digraph(MAP, n)
    cm = sort(countmap(vpairs.(getfield.(last.(samplers), :g))), byvalue=true, rev=true)
    @test first(cm).first == MAP
end #testset

@testset "MultiSampler" begin
    Random.seed!(1)

    N = 200 # number of data points

    # define simple linear model with added noise
    x = randn(N)
    v = x + randn(N)*0.25
    w = x + randn(N)*0.25
    z = v + w + randn(N)*0.25
    s = z + randn(N)*0.25
    
    df = (x=x, v=v, w=w, z=z, s=s)
    iterations = 500
    penalty = 2.0 # increase to get more edges in truth
    n = length(df) # vertices
    Random.seed!(101)
    C = cor(CausalInference.Tables.matrix(df))
    score = GaussianScore(C, N, penalty)
    M = 2000
    bestgraph, samplers = multisampler(n; M, σ=2.0, score, iterations)
    coldness = CausalInference.expcoldness(minimum(getfield.(last.(samplers), :τ)))

    gs = causalzigzag(n; score, κ=n-1, coldness, iterations)
    graphs, graph_pairs, hs, τs, ws, ts, scores = CausalInference.unzipgs(gs)
    posterior = sort(keyedreduce(+, graph_pairs, ws); byvalue=true, rev=true)


    # maximum aposteriori estimate
    MAP = [1=>2, 1=>3, 2=>1, 2=>4, 3=>1, 3=>4, 4=>5]
    @test bestgraph == digraph(MAP, n)
    cm = sort((proportionmap(vpairs.(getfield.(last.(samplers), :g)))), byvalue=true, rev=true)
    @test first(cm).first == MAP
    logΠ = map(g->score_dag(pdag2dag!(digraph(g, n)), score), collect(keys(cm)))
    Π = normalize(exp.(coldness*(logΠ .- maximum(logΠ))), 1)
    Πhat = normalize(collect(values(cm)), 1)
    @show coldness
    display([Π Πhat])
    s = 0.0
    for (i, k) in enumerate(keys(cm))
        s += posterior[k]
        #@show cm[k] Π[i] 
    end
    @show s
    @test s > 0.98
    @test norm(collect(values(cm)) - Π) < 0.02
end #testset
