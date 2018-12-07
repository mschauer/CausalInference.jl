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
