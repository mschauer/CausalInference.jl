using CausalInference, Graphs, MetaGraphs, Random

@testset "Plotting to text" begin
    g = DiGraph(5)
    for (i, j) in [(1, 2), (2, 3), (2, 4), (4, 2), (4, 5)]
        add_edge!(g, i, j)
    end

    io = IOBuffer()
    graph_to_text(io, g, collect(string.('a':'e')))
    @test strip(String(take!(io))) == "a->b  b->c  b->d  d->b  d->e"

    graph_to_text(io, g)
    @test strip(String(take!(io))) == "1->2  2->3  2->4  4->2  4->5"

    graph_to_text(io, g, edge_styles=Dict((2, 4) => "<->"))
    @test strip(String(take!(io))) == "1->2  2->3  2<->4  4->2  4->5"
end