using Distributions
using CausalInference
using Graphs
using DelimitedFiles
using LinearAlgebra
using Revise

if !@isdefined(pathNCI)
#    pathNCI = download("http://nugget.unisa.edu.au/ParallelPC/data/real/NCI-60.csv")
    pathNCI = download("https://raw.githubusercontent.com/vincentarelbundock/Rdatasets/master/csv/ISLR/NCI60.csv")
end
data, header = readdlm(pathNCI, ',', header=true)
X = data[1:end, 2:end-1]
d, n = size(X)
C = Symmetric(cor(X, dims = 2))
#h, s = @time skeleton(d, gausscitest, (C, n), 2.0)
#g1 = @time pcalg(d, gausscitest, (C, n), 2.0)
g2 = ges(d, GaussianScore(C, n, 2.0))
g2 = ges(d, GaussianScore(C, n, 2.0), parallel=true)
g2 = @time ges(d, GaussianScore(C, n, 2.0))
g2 = @time ges(d, GaussianScore(C, n, 2.0), parallel=true)

#println("inferred edges ", ne(h))
using ProfileView
@time g2 = @profview ges(d, GaussianScore(C, n, 2.0))
#@time g2 = @profview ges(d, GaussianScore(C, n, 2.0), parallel=true)

