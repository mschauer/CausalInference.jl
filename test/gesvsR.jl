using Distributions
using Graphs
using DelimitedFiles
using LinearAlgebra
using Test
using CausalInference
using Random
@testset "ges vs GES from pcalg (R)" begin 
    Random.seed!(1)
    data = [:nci60, :gmG8][2]

    score_R = Dict([:nci60 => -27741.897919271432, :gmG8 => -20106.287947087083])[data]

    if data == :nci60
        X = Matrix(Float64.(readdlm(joinpath(@__DIR__, "..", "nci60.csv"), ',')[2:end, 2:end]))
    elseif data == :gmG8
        X = Matrix(Float64.(readdlm(joinpath(@__DIR__,  "..", "gmG8.csv"), ',')[2:end, 2:end]))
    end
    n, d = size(X)
    C = Symmetric(cov(X, dims = 1, corrected=false))

    if data == :nci60
        g3 = DiGraph(map(x->x=="TRUE", readdlm(joinpath(@__DIR__,  "..", "nci60adj.csv"), ',')[2:end, 2:end]))
    else
        g3 = DiGraph(map(x->x=="TRUE", readdlm(joinpath(@__DIR__,  "..", "gmG8adj.csv"), ',')[2:end, 2:end]))
    end
    penalty = 1.0

    g2, s = ges(X; penalty, method=:gaussian_bic_raw)
    g2b, sb, (t1b, t2b) = ges(X; penalty)
    @test g2 == g2b
    @test s ≈ sb
    #g2c, sc, (t1c, t2c) = ges(X; penalty, parallel=true)
    @test score_R ≈ score_dag(DiGraph(d), GaussianScore(C, n, penalty)) + s
    @show score_R ≈ score_dag(pdag2dag!(copy(g2)), GaussianScore(C, n, penalty))
    @show score_R ≈ score_dag(pdag2dag!(copy(g3)), GaussianScore(C, n, penalty))

    @test isempty(symdiff(vpairs(g2), vpairs(g2b)))

    @test isempty(symdiff(vpairs(g2), vpairs(g3)))

end