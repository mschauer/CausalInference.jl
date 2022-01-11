using Distributions
using CausalInference
using Graphs
using DelimitedFiles
using LinearAlgebra


if !@isdefined(pathNCI)
    pathNCI = download("http://nugget.unisa.edu.au/ParallelPC/data/real/NCI-60.csv")
end
data = readdlm(pathNCI, ',')
d = size(data,2)
X = data[:,1:d]'
d, n = size(X)
C = Symmetric(cor(X, dims=2))
p = 0.01
@time h, s = skeleton(d, gausscitest, (C,n), quantile(Normal(), 1-p/2))
println("inferred edges ", ne(h))
