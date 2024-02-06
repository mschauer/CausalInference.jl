using Random, CausalInference, Statistics, Test, Graphs
@testset "Dag-Zig-Zag" begin
    Random.seed!(1)

    iterations = 500_000
    n = 3
    gs = dagzigzag(n; iterations);
    graphs, graph_pairs, hs, τs, ws, ts, scores = CausalInference.unzipgs(gs);
    posterior = sort(keyedreduce(+, graph_pairs, ws); byvalue=true, rev=true)


end #testset