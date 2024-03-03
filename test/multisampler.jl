using Random, CausalInference, StatsBase, Statistics, Test, Graphs, LinearAlgebra
@testset "MultiSampler" begin
    Random.seed!(1)

    N = 400 # number of data points

    # define simple linear model with added noise
    x = randn(N)
    v = x + randn(N)*0.25
    w = x + randn(N)*0.25
    z = v + w + randn(N)*0.25
    s = z + randn(N)*0.25
    
    df = (x=x, v=v, w=w, z=z, s=s)
    iterations = 1_000
    penalty = 2.0 # increase to get more edges in truth
    n = length(df) # vertices
    Random.seed!(101)
    C = cor(CausalInference.Tables.matrix(df))
    score = GaussianScore(C, N, penalty)
    decay = 1e-5
    schedule = (τ -> 1.0 + τ*decay, τ -> decay) # linear
    M = 20
    baseline = 0.0
    bestgraph, samplers = multisampler(n; M, ρ = 1.0, score, baseline, schedule, iterations)
    #posterior = sort(keyedreduce(+, graph_pairs, ws); byvalue=true, rev=true)

    # maximum aposteriori estimate
    MAP = [1=>2, 1=>3, 2=>1, 2=>4, 3=>1, 3=>4, 4=>5]
    @test bestgraph == digraph(MAP, n)
    cm = sort(countmap(vpairs.(getfield.(samplers, :g))), byvalue=true, rev=true)
    Tmin, T = extrema(getfield.(samplers, :τ))
    @show Tmin T schedule[1](T)
    @test first(cm).first == MAP
end #testset

@testset "MultiSampler" begin
    Random.seed!(1)
    decay = 1e-5
    schedule = (τ -> 1.0 + τ*decay, τ -> decay) # linear
  
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
    M = 100
    bestgraph, samplers = multisampler(n; M, score, schedule, iterations)
    Tmin, T = extrema(getfield.(samplers, :τ))
    coldness = schedule[1](T)
    @show Tmin T coldness

    gs = causalzigzag(n; score, κ=n-1, coldness, iterations=iterations*100)
    graphs, graph_pairs, hs, τs, ws, ts, scores = CausalInference.unzipgs(gs)
    posterior = sort(keyedreduce(+, graph_pairs, ws); byvalue=true, rev=true)

   
    # maximum aposteriori estimate
    MAP = [1=>2, 1=>3, 2=>1, 2=>4, 3=>1, 3=>4, 4=>5]
    @test bestgraph == digraph(MAP, n)
    cm = sort((proportionmap(vpairs.(getfield.(samplers, :g)))), byvalue=true, rev=true)
    @test first(cm).first == MAP
    logΠ = map(g->score_dag(pdag2dag!(digraph(g, n)), score), collect(keys(cm)))
    Π = normalize(exp.(coldness*(logΠ .- maximum(logΠ))), 1)
    Πhat = normalize(collect(values(cm)), 1)

    display([Π Πhat])
    s = 0.0
    for (i, k) in enumerate(keys(cm))
        s += get(posterior, k, 0.0)
        #@show cm[k] Π[i] 
    end
    @show s
    @test s > 0.98
    @test norm(collect(values(cm)) - Π) < 0.04
end #testset
