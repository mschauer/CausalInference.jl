using Random, CausalInference, Statistics, Test, Graphs
@testset "Zig-Zag" begin
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
    κ = n - 1 # max degree
    penalty = 2.0 # increase to get more edges in truth
    Random.seed!(101)
    C = cor(CausalInference.Tables.matrix(df))
    score = GaussianScore(C, N, penalty)
    gs = causalzigzag(n; score, κ, iterations)
    graphs, graph_pairs, hs, τs, ws, ts, scores = CausalInference.unzipgs(gs)
    posterior = sort(keyedreduce(+, graph_pairs, ws); byvalue=true, rev=true)

    # maximum aposteriori estimate
    @test first(posterior).first == [1=>2, 1=>3, 2=>1, 2=>4, 3=>1, 3=>4, 4=>5] 
    # score of last sample
    @test score_dag(pdag2dag!(copy(graphs[end])), score) ≈ scores[end] + score_dag(DiGraph(n), score)
    
    gs = causalzigzag(n; balance=CausalInference.sqrt_balance, score, κ, iterations, coldness = 100.0)
    g, s = ges(df; penalty=penalty, parallel=true)
    @test gs[end][1] == g
    @test gs[end][end] ≈ s
    @test gs[1][end] == 0
    @test ne(gs[1][1]) == 0
end #testset


@testset "CPDAG-Zig-Zag" begin
    Random.seed!(1)
      
    iterations = 4_000
    n = 3
    gs = causalzigzag(n; iterations);
    graphs, graph_pairs, hs, τs, ws, ts, scores = CausalInference.unzipgs(gs);
    posterior = sort(keyedreduce(+, graph_pairs, ws); byvalue=true, rev=true)
    @test length(posterior) == 11
    @test maximum(values(posterior)) < 1/11 + 0.02
    @test minimum(values(posterior)) > 1/11 - 0.02
    @show T_skew = sum(τs)

    gs = causalzigzag(n; iterations, σ=1.0, ρ=0.0);
    graphs, graph_pairs, hs, τs, ws, ts, scores = CausalInference.unzipgs(gs);
    posterior = sort(keyedreduce(+, graph_pairs, ws); byvalue=true, rev=true)
    @test length(posterior) == 11
    @test maximum(values(posterior)) < 1/11 + 0.02
    @test minimum(values(posterior)) > 1/11 - 0.02
    
    @show T = sum(τs)

    @test 1.98 < T_skew/T < 2.02
end 