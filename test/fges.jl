using Distributions, CausalInference, Graphs
using Test, Random, Statistics

# we use intersection of sorted arrays is sorted
@test issorted(intersect((1:1000)[rand(Bool, 1000)], (1:1000)[rand(Bool, 1000)]))
@test issorted(intersect((1:1000)[rand(Bool, 1000)], intersect((1:1000)[rand(Bool, 1000)], (1:1000)[rand(Bool, 1000)])))

Random.seed!(100)
@testset "GES " begin
    total = 0
    wrong0 = wrong1 = wrong2 = 0
    for n in 0:10
        alpha = 0.1
        @testset "randdag($n)" begin for k in 1:100
            global g = randdag(n, alpha)
            skel = DiGraph(Graph(g))
            h2 = cpdag(g)
            h1 = CausalInference.fges_internal(n, Float64, h2)
            wrong0 += !(vpairs(h1) ⊆ vpairs(skel))
            #h1 == h2 || println(vpairs(g))
            #@test vpairs(h1) ⊆ vpairs(h2)
            #@test vpairs(h2) ⊆ vpairs(h1)

            wrong1 += !(vpairs(h1) ⊆ vpairs(h2))
            wrong2 += !(vpairs(h2) ⊆ vpairs(h1))
            total += 1
        end end
    end
    println("$(wrong0/total) too large $(wrong1/total) too many $(wrong2/total) missing")
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
    seed = reinterpret(UInt, time())
    seed = 0x41d92f901f6981c3
    @show seed
    Random.seed!(seed)
    qu(x) = x * x'
    d = 5
    n = 2000
    alpha = 3/d
    c1 = 1.0
    c2 = 1.0
    penalty = 1.0

    E = Matrix(adjacency_matrix(randdag(d, alpha)))
    Random.seed!(7)
    L = E .* (0.9rand(d, d) .+ 0.1)
    println("\nVertices: $d, Edges: ", sum(E))
    
    X = (I - L) \ randn(d, n)
    C = cor(X, dims = 2)
    Σtrue = inv(qu((I - L)'))
    di = sqrt.(diag(Σtrue))
    Ctrue = (Σtrue) ./ (di * di')

    @show cond(I - L)
    g = DiGraph(E)
    cg = cpdag(g)
 
    @time h1, _ = skeleton(d, gausscitest, (C, n), c1)
    @time h2, _ = skeleton(d, dseporacle, g)
    println("skel edges ", ne(h2))
    println("skel missing edges ", length(setdiff(vpairs(h1), vpairs(h2))))
    println("skel wrong edges ", length(setdiff(vpairs(h2), vpairs(h1))))

    println("Vertices $(nv(g)) Edges $(ne(g))")

    g1 = @time fges(Matrix(X'), penalty=penalty)
    g2 = pcalg(d, gausscitest, (C, n), c2)

    println("ges not in skel ", length(setdiff(vpairs(g1), vpairs(DiGraph(h2)))))
    println("pc not in skel ", length(setdiff(vpairs(g2), vpairs(DiGraph(h2)))))

    println("cpdag dir ", ne(g) - (2ne(h2) - ne(cg)) , " undir ", 2ne(h2) - ne(cg))
    println("ges missing edges ", length(setdiff(vpairs(cg), vpairs(g1))))
    println("ges wrong edges ", length(setdiff(vpairs(g1), vpairs(cg))))
    println("pc missing edges ", length(setdiff(vpairs(cg), vpairs(g2))))
    println("pc wrong edges ", length(setdiff(vpairs(g2), vpairs(cg))))

end
                            
