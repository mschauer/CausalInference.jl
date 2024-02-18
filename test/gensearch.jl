using Graphs
using CausalInference
using Test
using CausalInference: graph
using Combinatorics: combinations

# for some of the tests other results would not be wrong per se 
# e.g. for the finding methods, other valid sets could be returned
# however, the algorithms should be deterministic, so we think the
# tests are okay for now



g1 = SimpleDiGraph(Edge.([(1, 2), (2, 3), (3, 4), (5, 1), (6, 5), (6, 4), (1, 7)]))
X = Set(1)
Xint = 1
Xvec = [1]
Y = Set(4)
Yint = 4
Yvec = [4]

@testset "gensearch in1" begin
    @test ancestors(g1, X) == Set([1,5,6])
    @test descendants(g1, Xint) == Set([1,2,3,4,7])
    @test alt_test_dsep(g1, Xvec, Yint, Set([2,5]))
    @test alt_test_dsep(g1, X, Yint, Set([3,5]))
    @test alt_test_dsep(g1, Xvec, Y, Set([3,6]))
    @test !alt_test_dsep(g1, X, Y, Set(2))
    @test !alt_test_dsep(g1, Xint, Y, Set(5))
    @test !alt_test_backdoor(g1, Xint, Yint, Set([5,7]))
    @test !alt_test_backdoor(g1, Xvec, Yvec, Set([2,5,7]))
    @test alt_test_backdoor(g1, X, Y, Set(5))
    @test test_covariate_adjustment(g1, Xvec, Yint, Set([5,7]))
    @test !test_covariate_adjustment(g1, X, Y, Set([2,5,7]))
    @test test_covariate_adjustment(g1, Xint, Yvec, Set(5))
    @test find_dsep(g1, Xint, Yint) == Set([2,3,5,6])
    @test find_dsep(g1, X, Set(2)) == false
    @test find_min_dsep(g1, X, Y) == Set([2,5])
    @test find_covariate_adjustment(g1, Xint, Yint) == Set([5,6])
    @test find_backdoor_adjustment(g1, X, Y) == Set([5,6])
    @test find_frontdoor_adjustment(g1, X, Y) == Set([2,3,7])
    @test find_min_covariate_adjustment(g1, Xvec, Yvec) ==  Set(5)
    @test find_min_backdoor_adjustment(g1, Xint, Yvec) == Set(5)
    @test find_min_frontdoor_adjustment(g1, X, Y) == Set(3)
    @test Set(list_dseps(g1, Xvec, Y)) == Set([Set([3,6]), Set([2,6]), Set([2,3,6]), Set([3,6,7]), Set([2,6,7]), Set([2,3,6,7]), Set([3,5]), Set([2,5]), Set([2,3,5]), Set([3,5,7]), Set([2,5,7]), Set([2,3,5,7]), Set([3,5,6]), Set([2,5,6]), Set([2,3,5,6]), Set([3,5,6,7]), Set([2,5,6,7]), Set([2,3,5,6,7])])
    @test Set(list_covariate_adjustment(g1, X, Yint)) == Set([Set([5,7]), Set([5]), Set([6]), Set([6,7]), Set([5,6]), Set([5,6,7])])
    @test Set(list_backdoor_adjustment(g1, Xvec, Yint)) == Set([Set([5]), Set([6]), Set([5,6])])
    @test Set(list_frontdoor_adjustment(g1, X, Y)) == Set([Set([3]), Set([7,2]), Set([7,2,3]), Set([7,3]), Set([2]), Set([2,3])]) 
end

g2 = SimpleDiGraph(Edge.([(1, 3), (3, 6), (2, 5), (5, 8), (6, 7), (7, 8), (1, 4), (2, 4), (4, 6), (4, 8)]))
X = Set(6)
Y = Set(8)
@testset "gensearch in2" begin
    @test ancestors(g2, X) == Set([1,2,3,4,6])
    @test descendants(g2, X) == Set([6,7,8])
    @test alt_test_dsep(g2, X, Y, Set([3,4,7]))
    @test alt_test_dsep(g2, X, Y, Set([4,5,7]))
    @test alt_test_dsep(g2, X, Y, Set([1,4,7]))
    @test !alt_test_dsep(g2, X, Y, Set(7))
    @test !alt_test_dsep(g2, X, Y, Set([4,7]))
    @test !alt_test_backdoor(g2, X, Y, Set([3,4,7]))
    @test !alt_test_backdoor(g2, X, Y, Set([3,5]))
    @test alt_test_backdoor(g2, X, Y, Set([4,2]))
    @test test_covariate_adjustment(g2, X, Y, Set([3,4]))
    @test !test_covariate_adjustment(g2, X, Y, Set([2,4,7]))
    @test test_covariate_adjustment(g2, X, Y, Set([5,4]))
    @test find_dsep(g2, X, Y) == Set([1,2,3,4,5,7])
    @test find_dsep(g2, X, Y, Set{Int64}(), setdiff(Set(1:8), [4,6,8])) == false
    @test find_dsep(g2, Set([1,6]), Set(2)) == false
    @test find_min_dsep(g2, X, Y) == Set([3,4,7])
    @test find_covariate_adjustment(g2, X, Y, Set(7), Set([1,2,3,4,5])) == false
    @test find_covariate_adjustment(g2, X, Y, Set{Int64}(), Set([3,4,5,7])) == Set([3,4,5]) 
    @test find_backdoor_adjustment(g2, X, Y) == Set([1,2,3,4,5])
    @test find_frontdoor_adjustment(g2, X, Y) == Set(7)
    @test find_min_covariate_adjustment(g2, X, Y) ==  Set([3,4])
    @test find_min_backdoor_adjustment(g2, X, Y) == Set([3,4])
    @test find_min_frontdoor_adjustment(g2, X, Y) == Set(7)
    @test Set(list_dseps(g2, X, Y, Set{Int64}(), Set{Int64}([3,4,5,7]))) == Set([Set([4, 7, 3]), Set([5, 4, 7]), Set([5, 4, 7, 3])])
    @test Set(list_covariate_adjustment(g2, Set([6]), Set([8]), Set(Int[]), setdiff(Set(1:8), [1,2]))) == Set([Set([3,4]), Set([4,5]), Set([3,4,5])])
    @test Set(list_backdoor_adjustment(g2, Set([6]), Set([8]), Set(Int[]), setdiff(Set(1:8), [1,2]))) == Set([Set([3,4]), Set([4,5]), Set([3,4,5])])
    @test Set(list_frontdoor_adjustment(g2, X, Y)) == Set([Set(7)]) 
end

using Test, CausalInference, Combinatorics
function test_dsep(g) 
    n = nv(g)
    for (_, v, w, z) in partitions(1:n, 4)
        @test dsep(g, v, w, z) == alt_test_dsep(g, v, w, z)
    end
end
@testset "dsep vs alt_test_dsep" begin
    test_dsep(g1)
    test_dsep(g2)
    @test_throws ArgumentError dsep(g1, [1,2], [2,3], [4,5])
    @test_throws ArgumentError dsep(g1, [1,2], [3,4], [4,5])
    @test_throws ArgumentError dsep(g1, [1,2], [3,4], [5,1])
end

@testset "bayesball_graph" begin
    g = digraph([1=>3, 2=>3, 3=>4, 2=>4, 1=>4])
    g2 = CausalInference.bayesball_graph(g, 2, [3], back=true)
    @test g2 == digraph([2 => 5, 2 => 7, 4 => 5, 4 => 7, 5 => 2, 5 => 4], 8)


    for d in 2:8
        g = randdag(d, 0.3)
        for S in combinations(1:d)
            Sᶜ = setdiff(1:d, S)
            for v in Sᶜ
                for back in (true, false)
                    g2 = CausalInference.bayesball_graph(g, v, S)
                    for w in Sᶜ
                        w == v && continue
                        @test !dsep(g, v, w, S) == has_path(g2, 2v, 2w-1) || has_path(g2, 2v, 2w) 
                    end
                end
            end
        end
    end
end