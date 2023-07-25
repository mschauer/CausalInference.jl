using Distributions, CausalInference, Graphs
using Test, Random, Statistics, Combinatorics

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
            h1 = CausalInference.ges_internal(n, Float64, h2)
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


g = ges(X)

@test sort(collect(Graphs.edges(g))) == sort(Edge.([1 => 2
                            1 => 3
                            2 => 1
                            2 => 4
                            3 => 1
                            3 => 4
                            4 => 5]))

using LinearAlgebra, StatsBase
# https://cran.r-project.org/web/packages/BCDAG/vignettes/bcdag_generatedata.html
#@testset "GES" begin
begin
    seed = reinterpret(UInt, time())
  #  seed = 0x41d92fae13925cd9
    @show seed
    Random.seed!(seed)
    qu(x) = x * x'
    d = 6
    n = 20000
    alpha = 2.5/d
    c1 = 0.001
    penalty = 0.01
    #g = digraph([1 => 2, 1 => 4, 2 => 4, 3 => 2, 3 => 4])
    g = randdag(d, alpha)
    E = Matrix(adjacency_matrix(g)) # Markov style
    L = E .* (0.3rand(d, d) .+ 0.3)
    println("\nVertices: $d, Edges: ", sum(E))
    
    # Do not actually sample
    #X = (I - L)' \ randn(d, n)
    #C = cor(X, dims = 2)
    Σtrue = Float64.(inv(big.(qu((I - L)))))
    di = sqrt.(diag(Σtrue))
    Ctrue = (Σtrue) ./ (di * di')
    println("Cond ", cond(Ctrue))
    @assert g == DiGraph(E)

    for i in 1:d
        for j in 1:d
            i == j && continue
            #I = setdiff(sample(1:d, rand(0:d), replace=false), [i,j])
            for I in powerset(setdiff(1:d, [i,j]))
                t = dseporacle(i, j, I, g) == gausscitest(i, j, I, (Ctrue, n), c1)
                if t == false
                    println(vpairs(g))
                    println("$i $j $I")
                    println(dsep(g, i, j, I), " ", CausalInference.partialcor(i, j, I, Ctrue))
                    error()
                end
            end
        end
    end

    cg = cpdag(g)
 
    h1, _ = skeleton(d, gausscitest, (Ctrue, n), c1)
    h2, _ = skeleton(d, dseporacle, g)
    println("skel edges ", ne(h2))
    println("skel missing edges ", length(setdiff(vpairs(h1), vpairs(h2))))
    println("skel wrong edges ", length(setdiff(vpairs(h2), vpairs(h1))))

    println("Vertices $(nv(g)) Edges $(ne(g))")

    #g1 = @time ges(Matrix(X'), penalty=penalty)
    g1 = @time CausalInference.ges_internal(d, Float64, GaussianScore(Ctrue, n, penalty))
    g2 = pcalg(d, gausscitest, (Ctrue, n), c1)

    println("ges not in skel ", length(setdiff(vpairs(g1), vpairs(DiGraph(h2)))))
    println("pc not in skel ", length(setdiff(vpairs(g2), vpairs(DiGraph(h2)))))

    println("cpdag undir ", ne(g) - (2ne(h2) - ne(cg)) , " dir ", 2ne(h2) - ne(cg))
    println("ges missing edges ", length(setdiff(vpairs(cg), vpairs(g1))))
    println("ges wrong edges ", length(setdiff(vpairs(g1), vpairs(cg))))
    println("pc missing edges ", length(setdiff(vpairs(cg), vpairs(g2))))
    println("pc wrong edges ", length(setdiff(vpairs(g2), vpairs(cg))))

end
                            
