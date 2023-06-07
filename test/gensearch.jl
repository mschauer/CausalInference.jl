using Graphs
using CausalInference
using Random
using Test
Random.seed!(1)

g = SimpleDiGraph(Edge.([(1, 2), (2, 3), (3, 4), (5, 1), (6, 5), (6, 4), (1, 7)]))
X = Set(1)
Y = Set(4)

@testset "gensearch in1" begin
    @test ancestors(g, X) == Set([1,5,6])
    @test descendants(g, X) == Set([1,2,3,4,7])
    @test alt_test_dsep(g, X, Y, Set([2,5]))
    @test alt_test_dsep(g, X, Y, Set([3,5]))
    @test alt_test_dsep(g, X, Y, Set([3,6]))
    @test !alt_test_dsep(g, X, Y, Set(2))
    @test !alt_test_dsep(g, X, Y, Set(5))
    @test !alt_test_backdoor(g, X, Y, Set())
    @test !alt_test_backdoor(g, X, Y, Set(2))
    @test alt_test_backdoor(g, X, Y, Set(5))
    @test find_dsep(g, X, Y) == Set([2,3,5,6])
    @test find_dsep(g, X, Set(2)) == false
    mindsepZ = find_min_dsep(g, X, Y)
    @test alt_test_dsep(g, X, Y, mindsepZ)
    for S in powerset(collect(mindsepZ))
        length(S) == length(result) && continue
        @test !alt_test_dsep(g, X, Y, S)
    end
    @test find_covariate_adjustment(g, X, Y) == Set([5,6])
    @test find_backdoor_adjustment(g, X, Y) == Set([5,6])
    @test find_min_covariate_adjustment(g, X, Y) ==  Set(5)
    @test find_frontdoor_adjustment(g, X, Y) == Set([2,3,7])
    @test  
end
