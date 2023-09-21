using GLMakie
using Colors
using GraphMakie
using CausalInference
using NetworkLayout

G = [:A=>:T, :T=>:E, :E=>:X, :L=>:E, :S=>:L, :B=>:D, :E=>:D]
#V = unique(Iterators.flatten(G))
V = [:B, :A, :T, :S, :L,  :E, :X, :D]
V⁻ = Dict(v=>i for (i,v) in enumerate(V))
g = digraph(map(e-> V⁻[e[1]]=>V⁻[e[2]], G))
m = cpdag(g)

gp(g; other...) = graphplot(g; kwargs_pdag_graphmakie(g)..., other...)
fig0, ax0, pl0 = gp(m, layout=Shell(), ilabels=V,  ilabels_fontsize = 50, node_color= colorant"#c7e1ca")
hidedecorations!.(ax0)
fig0

iterations = 50_000; verbose = false
n = nv(g) # vertices
κ = n - 1 # max degree
N = 50 # increase to get more concentrated posterior
alpha = 0.12 # increase to get more edges in truth
Random.seed!(101)

E = Matrix(adjacency_matrix(g)) # Markov operator multiplies from right 
L = E .* (0.3rand(n, n) .+ 0.3)
penalty = 2.0 # increase to get less edges in sample
Σtrue = Float64.(inv(big.(qu((I - L)))))
di = sqrt.(diag(Σtrue))
Ctrue = (Σtrue) ./ (di * di')
#score = GaussianScore(Ctrue, N, penalty)
score = UniformScore()

#if !(@isdefined gs)
    gs = @time randcpdag(n; score, ρ=1.0, σ=0.0, wien=true, κ, iterations, verbose)[burnin:end]
#end

graphs = first.(gs)
graph_pairs = as_pairs.(graphs)
hs = map(last, gs)
τs = map(x->getindex(x, 2), gs)
ws = normalize(τs, 1)
ts = cumsum(ws)

print(as_pairs(m))
cm = keyedreduce(+, graph_pairs, ws)
cm = sort(cm; byvalue=true, rev=true)



import CausalInference.kwargs_pdag_graphmakie
function CausalInference.kwargs_pdag_graphmakie(g::Observable; ilabels=1:nv(g[]), arrowsize=25, ilabels_fontsize=25)
    
    kwargs = @lift kwargs_pdag_graphmakie($g; ilabels, arrowsize, ilabels_fontsize)
    arrow_size = @lift ($kwargs).arrow_size
    edge_width = @lift ($kwargs).edge_width
    (; arrow_shift=:end, arrowsize, edge_width)
end
frames = 1000

g = Observable(graphs[1])

fig = Figure()

ax2 = fig[1,1] = Axis(fig)
hidedecorations!.(ax2)

ax1 = fig[2,1] = Axis(fig)

xlims!(ax1, 0.0, frames/length(ts))
tim = Observable(copy(ts))
hei = Observable(copy(hs))
resize!(tim[],1)
resize!(hei[],1)

stairs!(ax1, tim, hei, step=:post)
ylims!(ax1, 0, n*(n-1)÷2)

#xlims!(ax1, tim[max(1, length(hs) - 1000)], tim[end])
#ax1.yzoomlock = true


graphplot!(ax2, g, kwargs_pdag_graphmakie(g, arrowsize=100)..., curve_distance_usage = false, arrow_shift = :end, 
    layout=pl0[:node_pos][], ilabels=V,  ilabels_fontsize = 30, node_color= colorant"#c7e1ca")
 
xlims!(ax2, -1.3, 1.3)
ylims!(ax2, -1.3, 1.3)

display(fig)

record(fig, "append_animation2.mp4", 2:frames; framerate = 30) do i
    g[] = graphs[i]
    #tim[] =
     push!(tim[], ts[i])
    hei[] = push!(hei[], hs[i])
end

fig1 = fig