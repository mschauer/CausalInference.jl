using Distributions
using CausalInference
using LightGraphs

quantile(Normal(), 0.95)

if !isfile(joinpath("data","NCI-60.csv"))
    cd("data") do
        run(`wget http://nugget.unisa.edu.au/ParallelPC/data/real/NCI-60.csv`)
    end
end  

data = readdlm("data/NCI-60.csv")
d = size(data,2)
X = data[:,1:d]'
d, n = size(X)
C = Symmetric(cor(X, 2))
p = 0.01
@time h, s = skeleton(d, gausscitest, (C,n), quantile(Normal(), 1-p/2)) 
println("inferred edges ", ne(h))
writedlm("data/adjac.csv", adjacency_matrix(h))

