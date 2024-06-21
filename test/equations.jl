using CausalInference
using TikzGraphs
using Random
Random.seed!(1)

# Generate some sample data to use with the GES algorithm

N = 2000 # number of data points

# define simple linear model with added noise
x = randn(N)
v = x + randn(N)*0.25
w = x + randn(N)*0.25
z = v + w + randn(N)*0.25
s = z + randn(N)*0.25

df = (x=x, v=v, w=w, z=z, s=s)

est_g, score = ges(df; penalty=1.0, parallel=true)


est_dag= pdag2dag!(est_g)

scm= estimate_equations(df,est_dag)

display(scm)

#println(CI.SCM)

df_generated= generate_data(scm, 2000)

println("df: ")

display(df)

println("df_generated: ")



display(df_generated)