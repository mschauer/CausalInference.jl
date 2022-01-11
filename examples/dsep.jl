using CausalInference
using Graphs
using Test

include("plotdag.jl")

g1 = g = DiGraph(7)
d = nv(g)
for (i,j) in [(1,2), (2,3), (2,4), (4,5), (3,5), (5,6), (7,5)]
    add_edge!(g, i, j)
end

tp = plot_dag(g)
save(PDF("dag1"), tp)

@testset "dsep g1" begin
    @test !dsep(g, 1, 2, [])
    @test !dsep(g, 2, 5, [])

    @test !dsep(g, 1, 5, [3])
    @test !dsep(g, 1, 5, [4])
    @test dsep(g, 1, 5, [3, 4])

    @test dsep(g, 3, 4, [2])
    @test dsep(g, 3, 4, [2, 7])

    @test !dsep(g, 3, 4, [5])
    @test !dsep(g, 3, 4, [2, 6])


    @test !dsep(g, 3, 4, [2, 5])
    @test !dsep(g, 3, 4, [2, 6])

    @test !dsep(g, 3, 5, [6])
    @test !dsep(g, 3, 5, [7])

end



g2 = g = DiGraph(7)
d = nv(g)
for (i,j) in [(1,3), (2,3), (3,4),(3,5), (4,6), (6, 7)]
    add_edge!(g, i, j)
end

tp = plot_dag(g)
save(PDF("dag2"), tp)

@testset "dsep g2" begin
    @test !dsep(g, 1, 2, [4, 5])
    @test dsep(g, 1, 2, [])
    @test !dsep(g, 1, 2, [3])
    @test dsep(g, 4, 5, [3])
    @test !dsep(g, 4, 5, [])
    @test !dsep(g, 4, 5, [1, 2])
end
