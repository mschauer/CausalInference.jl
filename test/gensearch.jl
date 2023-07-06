using Graphs
using CausalInference
using Test

g1 = SimpleDiGraph(Edge.([(1, 2), (2, 3), (3, 4), (5, 1), (6, 5), (6, 4), (1, 7)]))
X = Set(1)
Y = Set(4)

# for some of the tests other results would not be wrong per se 
# e.g. for the finding methods, other valid sets could be returned
# however, the algorithms should be deterministic, so we think the
# tests are okay for now

@testset "gensearch in1" begin
    @test ancestors(g1, X) == Set([1,5,6])
    @test descendants(g1, X) == Set([1,2,3,4,7])
    @test alt_test_dsep(g1, X, Y, Set([2,5]))
    @test alt_test_dsep(g1, X, Y, Set([3,5]))
    @test alt_test_dsep(g1, X, Y, Set([3,6]))
    @test !alt_test_dsep(g1, X, Y, Set(2))
    @test !alt_test_dsep(g1, X, Y, Set(5))
    @test !alt_test_backdoor(g1, X, Y, Set([5,7]))
    @test !alt_test_backdoor(g1, X, Y, Set([2,5,7]))
    @test alt_test_backdoor(g1, X, Y, Set(5))
    @test test_covariate_adjustment(g1, X, Y, Set([5,7]))
    @test !test_covariate_adjustment(g1, X, Y, Set([2,5,7]))
    @test test_covariate_adjustment(g1, X, Y, Set(5))
    @test find_dsep(g1, X, Y) == Set([2,3,5,6])
    @test find_dsep(g1, X, Set(2)) == false
    @test find_min_dsep(g1, X, Y) == Set([2,5])
    @test find_covariate_adjustment(g1, X, Y) == Set([5,6])
    @test find_backdoor_adjustment(g1, X, Y) == Set([5,6])
    @test find_min_covariate_adjustment(g1, X, Y) ==  Set(5)
    @test find_min_backdoor_adjustment(g1, X, Y) == Set(5)
    @test find_frontdoor_adjustment(g1, X, Y) == Set([2,3,7])
    @test Set(list_dseps(g1, X, Y)) == Set([Set([3,6]), Set([2,6]), Set([2,3,6]), Set([3,6,7]), Set([2,6,7]), Set([2,3,6,7]), Set([3,5]), Set([2,5]), Set([2,3,5]), Set([3,5,7]), Set([2,5,7]), Set([2,3,5,7]), Set([3,5,6]), Set([2,5,6]), Set([2,3,5,6]), Set([3,5,6,7]), Set([2,5,6,7]), Set([2,3,5,6,7])])
    @test Set(list_covariate_adjustment(g1, X, Y)) == Set([Set([5,7]), Set([5]), Set([6]), Set([6,7]), Set([5,6]), Set([5,6,7])])
    @test Set(list_backdoor_adjustment(g1, X, Y)) == Set([Set([5]), Set([6]), Set([5,6])])
    @test Set(list_frontdoor_adjustment(g1, X, Y)) == Set([Set([3]), Set([7,2]), Set([7,2,3]), Set([7,3]), Set([2]), Set([2,3])]) 
end

g2 = SimpleDiGraph(Edge.([(1, 3), (3, 6), (2, 5), (5, 8), (6, 7), (7, 8), (1, 4), (2, 4), (4, 6), (4, 8)]))
X = Set(6)
Y = Set(8)
@testset "gensearch in2" begin
    @test Set(list_covariate_adjustment(g2, Set([6]), Set([8]), Set(Int[]), setdiff(Set(1:8), [1,2]))) == Set([Set([3,4]), Set([4,5]), Set([3,4,5])])
end
