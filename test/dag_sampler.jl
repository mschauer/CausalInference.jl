using Random, CausalInference, Statistics, Test, Graphs, LinearAlgebra
@testset "Dag-Zig-Zag" begin

    Random.seed!(1)

    N = 100 # number of data points

    # define simple linear model with added noise
    x = randn(N)
    v = x + randn(N)*0.25
    w = x + randn(N)*0.25
    z = v + w + randn(N)*0.25
    s = z + randn(N)*0.25

    df = (x=x, v=v, w=w, z=z, s=s)
    iterations = 10_000
    n = length(df) # vertices
    κ = n - 1 # max degree
    penalty = 2.0 # increase to get more edges in truth
    Random.seed!(101)
    C = cor(CausalInference.Tables.matrix(df))
    score = GaussianScore(C, N, penalty)
    for balance in (CausalInference.sqrt_balance, CausalInference.metropolis_balance)
        @testset "$balance" begin
            gs = dagzigzag(n; score, balance, κ, iterations, verbose=false)
            @test ne.(first.(gs)) == getindex.(gs, 4)
            # score of last sample
            graphs, graph_pairs, hs, τs, ws, ts, scores = CausalInference.unzipgs(gs)
            @test score_dag(graphs[end], score) ≈ scores[end] + score_dag(DiGraph(n), score)
            @test !any(Graphs.is_cyclic.(graphs))
            posterior = sort(keyedreduce(+, graph_pairs, ws); byvalue=true, rev=true)
            logΠ = map(g->score_dag(digraph(g, n), score), collect(keys(posterior)))
            Π = normalize(exp.(logΠ .- maximum(logΠ) ), 1)
            @test norm(collect(values(posterior)) - Π, 1) < 30/sqrt(iterations)

            posterior = sort(keyedreduce(+, vpairs.(cpdag.(graphs)), ws); byvalue=true, rev=true)
            @test first(posterior).first == [1=>2, 1=>3, 2=>1, 2=>4, 3=>1, 3=>4, 4=>5] 
        end
    end
end #testset





@testset "Dag-Zig-Zag" begin
    Random.seed!(1)
      
    iterations = 8_000
    n = 3
    gs = dagzigzag(n; iterations);
    graphs, graph_pairs, hs, τs, ws, ts, scores = CausalInference.unzipgs(gs);
    posterior = sort(keyedreduce(+, graph_pairs, ws); byvalue=true, rev=true)
    @test length(posterior) == 25
    @test maximum(values(posterior)) < 1/25 + 0.01
    @test minimum(values(posterior)) > 1/25 - 0.01
    

end #testset

