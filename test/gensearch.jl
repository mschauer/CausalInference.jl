using Graphs
using CausalInference
using Random
using Test
Random.seed!(1)

@testset "adjustment" begin
    g = digraph([1 => 3, 2 => 1, 2 => 3])
    X = Set{Integer}(1)
    Y = Set{Integer}(3)
    I = Set{Integer}()
    R = Set{Integer}(2)
    
#    @test TODO
end
