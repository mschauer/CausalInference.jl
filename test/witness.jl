using CausalInference
using LightGraphs
using Test
@testset "recanting witness" begin

    g = DiGraph(7)
    g_blocked = DiGraph(7)
    d = nv(g)
    for (i,j) in [(1,2), (2,3), (2,4), (4,5), (3,5), (5,6), (7,1)]
        add_edge!(g, i, j)
    end

    for (i,j) in [(3,5)]
        add_edge!(g_blocked, i, j)
    end
    @test has_recanting_witness(g, 1, 3, g_blocked) == false
    @test has_recanting_witness(g, 1, 4, g_blocked) == false
    @test has_recanting_witness(g, 1, 5, g_blocked) == true
    @test has_recanting_witness(g, 1, 6, g_blocked) == true

    for (i,j) in [(2,3),(3,5)]
        add_edge!(g_blocked, i, j)
    end
    @test has_recanting_witness(g, 1, 3, g_blocked) == false
    @test has_recanting_witness(g, 1, 5, g_blocked) == true
    @test has_recanting_witness(g, 1, 6, g_blocked) == true
end

