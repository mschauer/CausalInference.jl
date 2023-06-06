using Graphs
using CausalInference
using Random
using Test
Random.seed!(1)

g = SimpleDiGraph(Edge.([(1, 2), (2, 3), (4, 1), (4, 3)]))
X = Set{Integer}(1)
Y = Set{Integer}(3)

@testset "gensearch in1" begin
    @test issetequal(ancestors(g, X), Set{Integer}([1,4]))
    @test issetequal(descendants(g, X), Set{Integer}([1,2,3]))
    @test alt_test_dsep(g, X, Y, Set{Integer}([2,4]))
    @test !alt_test_dsep(g, X, Y, Set{Integer}(2))
    @test !alt_test_dsep(g, X, Y, Set{Integer}(4))
end
