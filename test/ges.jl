using Distributions, CausalInference, Graphs
using Test, Random, Statistics, Combinatorics
using LinearAlgebra, StatsBase

@testset "Sorted intersection" begin
    # we use intersection of sorted arrays is sorted
    @test issorted(intersect((1:1000)[rand(Bool, 1000)], (1:1000)[rand(Bool, 1000)]))
    @test issorted(intersect((1:1000)[rand(Bool, 1000)], intersect((1:1000)[rand(Bool, 1000)], (1:1000)[rand(Bool, 1000)])))
end

@testset "GES example" begin
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


    g, _ = ges(X; penalty=2.0)

    @test sort(vpairs(g)) == sort(([1 => 2
                            1 => 3
                            2 => 1
                            2 => 4
                            3 => 1
                            3 => 4
                            4 => 5]))
end

if Threads.nthreads() == 1
    println("Skipping parallel tests.")
end
# https://cran.r-project.org/web/packages/BCDAG/vignettes/bcdag_generatedata.html
@testset "GES randdag" begin
  #  seed = reinterpret(UInt, time())
    seed = 123
    @show seed
    Random.seed!(seed)
    t0 = t1 = t2 = t3 = 0.0
    K = 50
    for k in 1:K
        qu(x) = x * x'
        d = 8
        n = 20000
        alpha = 2.5/d
        c1 = 0.0001 # take small as no sampling error, only numerical errors
        penalty = 0.00005 # same
        #g = digraph([1 => 2, 1 => 4, 2 => 4, 3 => 2, 3 => 4])
        g = randdag(d, alpha)
        E = Matrix(adjacency_matrix(g)) # Markov operator multiplies from right 
        L = E .* (0.3rand(d, d) .+ 0.3)
        
        # Do not actually sample but compute true correlation
        #X = (I - L)' \ randn(d, n)
        #C = cor(X, dims = 2)
        Σtrue = Float64.(inv(big.(qu((I - L)))))
        di = sqrt.(diag(Σtrue))
        Ctrue = (Σtrue) ./ (di * di')
        @assert g == DiGraph(E)
        if k < 5
            t0 += @elapsed for i in 1:d
                for j in 1:d
                    i == j && continue
                    #I = setdiff(sample(1:d, rand(0:d), replace=false), [i,j])
                    for I in powerset(setdiff(1:d, [i,j]))
                        test_dsep = dseporacle(i, j, I, g) == gausscitest(i, j, I, (Ctrue, n), c1)
                        if !test_dsep
                            println(vpairs(g))
                            println("$i $j $I")
                            println(dsep(g, i, j, I), " ", CausalInference.partialcor(i, j, I, Ctrue))
                        end
                        @test test_dsep
                    end
                end
            end
        end

        cg = cpdag(g)
    
        h1, _ = skeleton(d, gausscitest, (Ctrue, n), c1)
        h2, _ = skeleton(d, dseporacle, g)

        t1 += @elapsed g1, _ = ges(d, GaussianScore(Ctrue, n, penalty))
        t2 += @elapsed g2 = pcalg(d, gausscitest, (Ctrue, n), c1)
        if Threads.nthreads() > 1
            t3 += @elapsed g3, _ = ges(d, GaussianScore(Ctrue, n, penalty), parallel=true)
            @test g1 == g3
        else
            @test_skip g1 == g3
        end
        test_pcges = g1 == g2
        if !test_pcges
            println("cond(C) ", cond(Ctrue))
    
            println("skel edges ", ne(h2))
            println("skel missing edges ", length(setdiff(vpairs(h1), vpairs(h2))))
            println("skel wrong edges ", length(setdiff(vpairs(h2), vpairs(h1))))
    
            println("vertices $(nv(g)) edges $(ne(g))")

            println("ges not in skel ", length(setdiff(vpairs(g1), vpairs(DiGraph(h2)))))
            println("pc not in skel ", length(setdiff(vpairs(g2), vpairs(DiGraph(h2)))))

            println("cpdag undir ", ne(g) - (2ne(h2) - ne(cg)) , " dir ", 2ne(h2) - ne(cg))
            println("ges missing edges ", length(setdiff(vpairs(cg), vpairs(g1))))
            println("ges wrong edges ", length(setdiff(vpairs(g1), vpairs(cg))))
            println("pc missing edges ", length(setdiff(vpairs(cg), vpairs(g2))))
            println("pc wrong edges ", length(setdiff(vpairs(g2), vpairs(cg))))
        end
        @test test_pcges
    end
    println("timing: PC ", t2/K, " GES ", t1/K, " GES-P ", t3/K, " (", t0/K, ")")
end
                            
