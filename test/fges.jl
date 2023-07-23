using Distributions, CausalInference, Graphs
using Test, Random

# we use intersection of sorted arrays is sorted
@test issorted(intersect((1:1000)[rand(Bool, 1000)], (1:1000)[rand(Bool, 1000)]))
@test issorted(intersect((1:1000)[rand(Bool, 1000)], intersect((1:1000)[rand(Bool, 1000)], (1:1000)[rand(Bool, 1000)])))

Random.seed!(100)
@testset "GES " begin
    wrong = 0
    for n in 0:10
        alpha = 0.1
        @testset "randdag($n)" begin for k in 1:100
            global g = randdag(n, alpha)
            h2 = cpdag(g)
            h1 = CausalInference.fges_internal(g, Float64, h2)
        
            #h1 == h2 || println(vpairs(g))
            @test vpairs(h1) ⊆ vpairs(h2)
            @test_skip vpairs(h2) ⊆ vpairs(h1)
            wrong += !(vpairs(h2) ⊆ vpairs(h1))
        end end
    end
    println("Wrong: $wrong")
end

Random.seed!(100)
N = 2000
d = Normal()
p = 0.01

x = randn(N)
v = x + randn(N)*0.5
w = x + randn(N)*0.5
z = v + w + randn(N)*0.5
s = z + randn(N)*0.5
X = [x v w z s]


g = fges(X)

@test collect(Graphs.edges(g)) == map(x -> Edge(x...), [1 => 2
                            1 => 3
                            2 => 1
                            2 => 4
                            3 => 1
                            3 => 4
                            4 => 5])


