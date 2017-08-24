
using CausalInference
using Base.Test

function skeleton2(n::V, I, par...) where {V}
    g = CompleteGraph(n)
    S = Dict{edgetype(g),Vector{V}}()
    d = 0 # depth
    while true
        isdone = true
        for e in collect(edges(g)) # cannot remove edges while iterating
            n = copy(neighbors(g, src(e)))::Vector{V}
            append!(n, neighbors(g, dst(e)))
            sort!(n)
            if length(n) > d  # i.e. |n\{dst(e)}| >= d 
                removesorted!(n, dst(e))
                removesorted!(n, src(e))

                isdone = false
                for s in combinations(n, d)
                    
                    if I(src(e), dst(e), s, par...) 
                        rem_edge!(g, e)
                        if !(e in keys(S))
                            S[e] = s
                        end
                        break 
                    end
                end
            end
        end 
        d = d + 1
        if isdone
            return g, S
        end
    end    
end

qu(x) = x*x'
 
@testset "d10" begin
    d = 10 # 40 disconnected
    n = 10000
    alpha = 0.3
    srand(5)
    E = LowerTriangular([i > j ? 1*(rand() < alpha) : 0 for i in 1:d, j in 1:d]) 
    L = E .* rand(d,d)
    println("\nVertices: $d, Edges: ", sum(E))

    X = (I - L)\randn(d, n)
    Σ = cov(X,2)
    Σtrue = inv(qu((I-L)'))

    @test norm(Σ -  Σtrue)  < .2*d*d/sqrt(n)
    di = sqrt.(diag(Σtrue));
    Ctrue = (Σtrue)./(di*di');
    C = cor(X,2)

    gd = DiGraph(E)
    g = Graph(Symmetric(E+E'))
    
    @time h, s = skeleton(d, dseporacle, gd)
    @test g == h


    println("Using true correlation")
    @time h, s = skeleton(d, gausscitest, (Ctrue, n*n), 0.05)
    O = full.(adjacency_matrix.(g)) 
    a = O +  2*full.(adjacency_matrix.(h))
    println("num edges found ", div(sum(a .== 3),2), " of ", ne(g), ", false edges ", div(sum(a .== 2),2) )

    println("Using data (n = $n)")
    @time h, s = skeleton(d, gausscitest, (C,n), 2.5)

    a = O +  2*full.(adjacency_matrix.(h))
    println("num edges found ", div(sum(a .== 3),2), " of ", ne(g), ", false edges ", div(sum(a .== 2),2) )
end

@testset "d23" begin
    d = 23 # 40 disconnected
    n = 100000
    alpha = 0.2
    srand(5) 
    E = LowerTriangular([i > j ? 1*(rand() < alpha) : 0 for i in 1:d, j in 1:d]) 
    L = E .* rand(d,d)
    println("\nVertices: $d, Edges: ", sum(E))

    X = (I - L)\randn(d, n)
    Σ = cov(X,2)
    Σtrue = inv(qu((I-L)'))

    @test norm(Σ -  Σtrue)  < .2*d*d/sqrt(n)
    di = sqrt.(diag(Σtrue));
    Ctrue = (Σtrue)./(di*di');
    C = cor(X,2)
    g = Graph(Symmetric(E+E'))
    gd = DiGraph(E)
    
    @time h, s = skeleton(d, dseporacle, gd)
    @test g == h

    println("Using true correlation")
    @time h, s = skeleton(d, gausscitest, (Ctrue, n*n), 0.05)
    O = full.(adjacency_matrix.(g)) 
    a = O +  2*full.(adjacency_matrix.(h))
    println("num edges found ", div(sum(a .== 3),2), " of ", ne(g), ", false edges ", div(sum(a .== 2),2) )

    println("Using data (n = $n)")
    @time h, s = skeleton(d, gausscitest, (C,n), 1.96)

    a = O +  2*full.(adjacency_matrix.(h))
    println("num edges found ", div(sum(a .== 3),2), " of ", ne(g), ", false edges ", div(sum(a .== 2),2) )

end


@testset "d5" begin
     g = DiGraph(5)
    d = nv(g)
    for (i,j) in [(1,2), (2,3), (2,4),(4,5), (3,5)]
        add_edge!(g,i,j)
    end
    
    h, s = skeleton(d, dseporacle, g)
    @test Graph(g) == h
    h, s = skeleton(d, truetest)
    @test ne(h) == 0
    
end    

@testset "truefalse" begin 
    d = 15
    @time h, s = skeleton(d, falsetest)
    @test ne(h) == div(d*(d-1),2)
     d = 100
    @time h, s = skeleton(d, truetest)
    @test ne(h) == 0
end    
