using CausalInference, Distributions, Random

@testset "Entropy" begin
    Random.seed!(123)
    
    @test isapprox(n_ball(2), π)
    @test isapprox(n_ball(3), 4/3. * π)
    N = 250000

    d = Normal(0,1)
    samples = rand(d, N)
    est = kl_entropy(collect(transpose(samples)))
    ent = entropy(d)
    @test (est-ent)/ent < 0.01

    d = Gamma(2,3)
    samples = rand(d, N)
    est = kl_entropy(collect(transpose(samples)))
    ent = entropy(d)
    @test (est-ent)/ent < 0.01

    d = Exponential(2.5)
    samples = rand(d, N)
    est = kl_entropy(collect(transpose(samples)))
    ent = entropy(d)
    @test (est-ent)/ent < 0.01

    r = 0.75
    d = MvNormal([1 r; r 1])
    samples = rand(d, N)
    X = collect(transpose(samples[1,:]))
    Y = collect(transpose(samples[2,:]))
    est = kl_mutual_information(X,Y)
    MI = -1/2. * log(1-r^2)
    @test (est-MI)/MI < 0.01
end

@testset "Permutation tests" begin
    Random.seed!(123)
    N = 10000
    p = 0.05
    x = randn(N)
    v = x + randn(N)*0.25
    w = x + randn(N)*0.25
    z = v + w +randn(N)*0.25

    @test kl_perm_mi_test(collect(transpose(v)), collect(transpose(w)))<p
    @test kl_perm_cond_mi_test(collect(transpose(v)),
                               collect(transpose(w)),
                               collect(transpose(z)))<p
    @test kl_perm_cond_mi_test(collect(transpose(v)),
                               collect(transpose(w)),
                               collect(transpose(x)))>p
    
end

@testset "Categorical Variables" begin
    Random.seed!(123)
    N = 1000
    p = 0.05
    
    x = rand(1:5, N)
    @test abs((cat_H(x) - log(5))/log(5)) < 0.01

    x = rand(1:4, N)
    v = map(d->floor(d/2), x) .+ rand(1:2, N)
    w = map(d->floor(d/2), x) .+ rand(1:2, N)
    @test perm_cat_MI_test(v,w) < p
    @test perm_cat_CMI_test(v,w,x) > p
end
