qu(x) = x*x'
function make_gaussian_model(n, N = 50, # increase to get more concentrated posterior
                            alpha = 0.12, # increase to get more edges in truth
                            penalty = 2.0 # increase to get less edges in sample
                            )
    g = randdag(n, alpha)
    E = Matrix(adjacency_matrix(g)) # Markov operator multiplies from right 
    L = E .* (0.3rand(n, n) .+ 0.3)
    Σtrue = Float64.(inv(big.(qu((I - L)))))
    di = sqrt.(diag(Σtrue))
    Ctrue = (Σtrue) ./ (di * di')
    cpdag(g), GaussianScore(Symmetric(Ctrue), N, penalty)
end