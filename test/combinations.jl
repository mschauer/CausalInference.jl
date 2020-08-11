using Combinatorics
using CausalInference
using Test
using Random

Random.seed!(2)

@testset "combinations" begin
    for n in 1:10
        r = shuffle(1:n)
        for k in 0:n+1
            w = rand(1:n)

            c1 = collect(copy(c) for c in combinations(r, k) if r[w] ∉ c)
            comb = CausalInference.combinations_without(r, k, w)
            c2 = collect(copy(c) for c in comb)

            @test length(c1) == length(comb)
            @test sort(c1) == sort(c2)
        end
    end
end
