using Random, CausalInference, Statistics, Test, Graphs
@testset "Dag-Zig-Zag" begin

    Random.seed!(1)

    N = 20 # number of data points

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
    gs = dagzigzag(n; score, balance=CausalInference.sqrt_balance, κ, iterations, verbose=false)
    @test ne.(first.(gs)) == getindex.(gs, 4)
    # score of last sample
    graphs, graph_pairs, hs, τs, ws, ts, scores = CausalInference.unzipgs(gs)
    @test score_dag(graphs[end], score) ≈ scores[end] + score_dag(DiGraph(n), score)
    @test !any(Graphs.is_cyclic.(graphs))
    posterior = sort(keyedreduce(+, graph_pairs, ws); byvalue=true, rev=true)
    
    posterior = sort(keyedreduce(+, vpairs.(cpdag.(graphs)), ws); byvalue=true, rev=true)
    @test first(posterior).first == [1=>2, 1=>3, 2=>1, 2=>4, 3=>1, 3=>4, 4=>5] 
  
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

