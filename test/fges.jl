using Distributions, CausalInference, Graphs
using Test, Random, Statistics

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

@test sort(collect(Graphs.edges(g))) == sort(Edge.([1 => 2
                            1 => 3
                            2 => 1
                            2 => 4
                            3 => 1
                            3 => 4
                            4 => 5]))

using LinearAlgebra
# https://cran.r-project.org/web/packages/BCDAG/vignettes/bcdag_generatedata.html
@testset "GES" begin
    #seed = reinterpret(UInt, time())
    seed = 0x41d92f73c5dfc9a7
    @show seed
    Random.seed!(seed)
    qu(x) = x * x'
    d = 6
    n = 40000
    alpha = 0.5
    p = 0.05
    c = 0.5

    #E = LowerTriangular([i > j ? 1 * (rand() < alpha) : 0 for i in 1:d, j in 1:d])
    E = Matrix(adjacency_matrix(randdag(d, 3, true)))
    L = E .* (rand(d, d) .- 0.5)
    println("\nVertices: $d, Edges: ", sum(E))

    X = (I - L) \ randn(d, n)
    C = cor(X, dims = 2)
    Σtrue = inv(qu((I - L)'))
    di = sqrt.(diag(Σtrue))
    Ctrue = (Σtrue) ./ (di * di')

    @show cond(I - L)
    g = DiGraph(E)
    cg = cpdag(g)
 
    @time h1, _ = skeleton(d, gausscitest, (C, n), c)
    @time h2, _ = skeleton(d, dseporacle, g)
    println("skel edges ", ne(h2))
    println("skel missing edges ", length(setdiff(vpairs(h1), vpairs(h2))))
    println("skel wrong edges ", length(setdiff(vpairs(h2), vpairs(h1))))

    println("Vertices $(nv(g)) Edges $(ne(g))")

    g1 = @time fges(Matrix(X'))
    g2 = pcalg(d, gausscitest, (C, n), c)

    println("ges not in skel ", length(setdiff(vpairs(g1), vpairs(DiGraph(h2)))))
    println("pc not in skel ", length(setdiff(vpairs(g2), vpairs(DiGraph(h2)))))

    println("cpdag dir ", ne(g) - (2ne(h2) - ne(cg)) , " undir ", 2ne(h2) - ne(cg))
    println("ges missing edges ", length(setdiff(vpairs(cg), vpairs(g1))))
    println("ges wrong edges ", length(setdiff(vpairs(g1), vpairs(cg))))
    println("pc missing edges ", length(setdiff(vpairs(cg), vpairs(g2))))
    println("pc wrong edges ", length(setdiff(vpairs(g2), vpairs(cg))))

end
                            
