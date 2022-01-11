using CausalInference
using Graphs
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

@testset "recanting witness azt example" begin
    g = DiGraph(6)
    d = nv(g)
    for (i,j) in [(1,2), (1,3), (1,6), (2,4), (3,5), (3,6), (4,6), (5,6)]
        add_edge!(g, i, j)
    end

    g_blocked = DiGraph(6) # 2a
    for (i,j) in [(1,2),(1,3)]
        add_edge!(g_blocked, i, j)
    end
    @test has_recanting_witness(g, 1, 6, g_blocked) == false

    g_blocked = DiGraph(6) # 2b
    for (i,j) in [(1,6)]
        add_edge!(g_blocked, i, j)
    end
    @test has_recanting_witness(g, 1, 6, g_blocked) == false

    g_blocked = DiGraph(6) # 3a
    for (i,j) in [(5,6)]
        add_edge!(g_blocked, i, j)
    end
    @test has_recanting_witness(g, 1, 6, g_blocked) == true

    g_blocked = DiGraph(6) # 3b
    for (i,j) in [(4,6)]
        add_edge!(g_blocked, i, j)
    end
    @test has_recanting_witness(g, 1, 6, g_blocked) == false

    
end

@testset "recanting witness boundary" begin
    g = DiGraph(4)
    d = nv(g)
    for (i,j) in [(1,2), (2,4), (1,3), (3,4)]
        add_edge!(g, i, j)
    end

    g_blocked = DiGraph(4)
    for (i,j) in [(1,2)]
        add_edge!(g_blocked, i, j)
    end
    @test has_recanting_witness(g, 1, 4, g_blocked) == false

    g_blocked = DiGraph(4)
    for (i,j) in [(2,4)]
        add_edge!(g_blocked, i, j)
    end
    @test has_recanting_witness(g, 1, 4, g_blocked) == false

end


@testset "recanting witness kite example" begin
    g = DiGraph(5)
    d = nv(g)
    for (i,j) in [(1,2),(2,3),(2,4),(3,5),(4,5)]
        add_edge!(g, i, j)
    end

    g_blocked = DiGraph(5)
    for (i,j) in [(2,3)]
        add_edge!(g_blocked, i, j)
    end
    @test has_recanting_witness(g, 1, 5, g_blocked) == true

    g_blocked = DiGraph(5)
    for (i,j) in [(2,4)]
        add_edge!(g_blocked, i, j)
    end
    @test has_recanting_witness(g, 1, 5, g_blocked) == true

    g_blocked = DiGraph(5)
    for (i,j) in [(1,2)]
        add_edge!(g_blocked, i, j)
    end
    @test has_recanting_witness(g, 1, 5, g_blocked) == false
end