using Distributions
using Graphs
using DelimitedFiles
using LinearAlgebra
using Revise
using CausalInference
data = [:nci60, :gmG8][1]
#=
if !@isdefined(pathNCI)
#    pathNCI = download("http://nugget.unisa.edu.au/ParallelPC/data/real/NCI-60.csv")
    pathNCI = download("https://raw.githubusercontent.com/vincentarelbundock/Rdatasets/master/csv/ISLR/NCI60.csv")
    data, header = readdlm(pathNCI, ',', header=true)
    X = Matrix(data[1:end, 2:end-1]')    
end=#
if data == :nci60
    X = Matrix(Float64.(readdlm("nci60.csv", ',')[2:end, 2:end]))
else
    X = Matrix(Float64.(readdlm("gmG8.csv", ',')[2:end, 2:end]))
end
n, d = size(X)
C = Symmetric(cov(X, dims = 1, corrected=false))
#h, s = @time skeleton(d, gausscitest, (C, n), penalty)
#g1 = @time pcalg(d, gausscitest, (C, n), penalty)
if data == :nci60
    g3 = DiGraph(map(x->x=="TRUE", readdlm("nci60adj.csv", ',')[2:end, 2:end]))
else
    g3 = DiGraph(map(x->x=="TRUE", readdlm("gmG8adj.csv", ',')[2:end, 2:end]))
end
penalty = 1.0
#g2b, sb, (t1b, t2b) = @time ges(X; penalty, verbose=true)
#error("done")
ges(X; penalty)

#g2, s = @time ges(X; penalty, method=:gaussian_bic_raw)
CausalInference.Memoization.empty_all_caches!()
g2b, sb, (t1b, t2b) = @time ges(X; penalty)
CausalInference.Memoization.empty_all_caches!()
g2c, sc, (t1c, t2c) = @time ges(X; penalty, parallel=true)
@show score_dag(DiGraph(d), GaussianScore(C, n, penalty)) + s
@show score_dag(pdag2dag!(copy(g2)), GaussianScore(C, n, penalty))
@show score_dag(pdag2dag!(copy(g3)), GaussianScore(C, n, penalty))

@show length(symdiff(vpairs(g2), vpairs(g2b)))

@show length(symdiff(vpairs(g2), vpairs(g3)))

error("Done.")
g2, s = ges(d, GaussianScore(C, n, penalty), parallel=true)
g2, s = @time ges(d, GaussianScore(C, n, penalty))
g2, s = @time ges(d, GaussianScore(C, n, penalty), parallel=true)

#println("inferred edges ", ne(h))
using ProfileView
@time g2 = @profview ges(d, GaussianScore(C, n, penalty))
#@time g2 = @profview ges(d, GaussianScore(C, n, penalty), parallel=true)

