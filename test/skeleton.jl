

srand(5)
let 
    d = 20 # 40 disconnected
    n = 100000
    alpha = 0.2
    qu(x) = x*x'
    E = LowerTriangular([i > j ? 1*(rand() < alpha) : 0 for i in 1:d, j in 1:d]) 
    L = E .* rand(d,d)
    println("Edges: ", sum(E))


    X = (I - L)\randn(d, n)
    Σ = cov(X,2)
    Σtrue = inv(qu((I-L)'))

    @test norm(Σ -  Σtrue)  < .2*d*d/sqrt(n)
    di = sqrt.(diag(Σtrue));
    Ctrue = (Σtrue)./(di*di');
    C = cor(X,2)
    g = Graph(Symmetric(E+E'))
    #@time h, s = skeleton(d, oracle, g)
    #@test g == h

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

