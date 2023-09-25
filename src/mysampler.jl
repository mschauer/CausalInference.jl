using Revise
using LinearAlgebra
using CausalInference
using CausalInference.Graphs
using StatsBase
using ProgressMeter
const as_pairs = vpairs
using Random
using Distributions
using OffsetArrays
using GLMakie
include("dags20.jl")
Random.seed!(5)
balance(t) = min(one(t), t)
#balance(t) = sqrt(t)
#balance(t) = t/(1+t)
ncpdags = [1, 2, 11, 185, 8782, 1067825, 312510571, 212133402500, 326266056291213, 1118902054495975181, 8455790399687227104576, 139537050182278289405732939, 4991058955493997577840793161279]

const coldness = 1.0

qu(x) = x*x'


using CausalInference: isadjacent, tails_and_adj_neighbors, adj_neighbors, isclique, Insert!, Delete!
using CausalInference: isundirected, parents, meek_rule1, meek_rule2, meek_rule3, meek_rule4, children, 
    neighbors_adjacent, neighbors_undirected, Δscoreinsert, Δscoredelete, keyedreduce
using CausalInference.Combinatorics
const adjacents = neighbors_adjacent

#g1 = SimpleDiGraph(Edge.([(1, 2), (2, 1), (1, 3), (3, 1), (1, 4), (4, 1), (2, 4), (4, 2), (3, 4), (4, 3)]))
#println("start")
#c = countcliques(g1)
#println(c)
#bigc = Dict{Vector{Int64}, Int64}()
#rep = 1000000
#for i=1:rep
#    ret = sampleclique(g1, rand(1:c))
#    if haskey(bigc, ret)
#        bigc[ret] += 1
#    else
#        bigc[ret] = 1
#    end
#end
#
#println(bigc)
#
#println("end")

iterations = 10000; verbose = false
n = 20 # vertices
κ = n - 1 # max degree
#κ = 4
reversible_too = true # do baseline 
#iterations = 50; verbose = true
#burnin = iterations÷2
burnin = 1
uniform = false
verbose = false

if uniform # sample uniform
    score = UniformScore()
elseif n == 5 # infer example data from https://mschauer.github.io/CausalInference.jl/latest/examples/ges_basic_examples/
    score = let N = 200
        Random.seed!(100)
        x = randn(N)
        v = x + randn(N)*0.5
        w = x + randn(N)*0.5
        z = v + w + randn(N)*0.5
        s = z + randn(N)*0.5
        X = [x v w z s]
        penalty = 0.5
        C = Symmetric(cov(X, dims = 1, corrected=false))
        GaussianScore(C, N, penalty)
    end
    true_cpdag = [1 => 2, 1 => 3, 2 => 1, 2 => 4, 3 => 1, 3 => 4, 4 => 5]
elseif n == 4
    score = let N = 200
        Random.seed!(100)
        v = randn(N)*0.5
        w = randn(N)*0.5
        z = v + w + randn(N)*0.5
        s = z + randn(N)*0.5
        X = [v w z s]
        penalty = 0.5
        C = Symmetric(cov(X, dims = 1, corrected=false))
        GaussianScore(C, N, penalty)
    end
    true_cpdag = [1 => 3, 2 => 3, 3 => 4]
elseif n == 3
    score = let N = 200
        Random.seed!(100)
        v = randn(N)*0.5
        w = randn(N)*0.5
        z = v + w + randn(N)*0.5
        X = [v w z]
        penalty = 0.5
        C = Symmetric(cov(X, dims = 1, corrected=false))
        GaussianScore(C, N, penalty)
    end
    true_cpdag = [1 => 2, 1 => 3]
else # 
    score, true_cpdag = let N = 50 # increase to get more concentrated posterior
        alpha = 0.12 # increase to get more edges in truth
        Random.seed!(101)
        g = randdag(n, alpha)
        E = Matrix(adjacency_matrix(g)) # Markov operator multiplies from right 
        L = E .* (0.3rand(n, n) .+ 0.3)
        penalty = 2.0 # increase to get less edges in sample
        Σtrue = Float64.(inv(big.(qu((I - L)))))
        di = sqrt.(diag(Σtrue))
        Ctrue = (Σtrue) ./ (di * di')
        GaussianScore(Symmetric(Ctrue), N, penalty), as_pairs(cpdag(g))
    end
end 
#G = (complete_digraph(n), n*κ÷2) 
G = DiGraph(n), 0
gs = @time randcpdag(n, G; score, ρ=1.0, σ=0.0, naive=false, κ, iterations, verbose)[burnin:end]

graphs = first.(gs)
graph_pairs = as_pairs.(graphs)
hs = hsnonrev = map(last, gs)
undirecteds = ne.(graphs) - hs


i = rand(eachindex(hs))
g, h = graphs[i], hs[i]

τs = map(x->getindex(x, 2), gs)
ws = wsnonrev = normalize(τs, 1)

println("Average nr. of undirected edges: " , sum(undirecteds .* ws))


@show sum(τs)
cm = keyedreduce(+, graph_pairs, ws)
cm = sort(cm; byvalue=true, rev=true)

if reversible_too
    gsrev = @time randcpdag(n, G; score, ρ=0.0, σ=1.0, naive=false, κ, iterations, verbose)[burnin:end]
    hsrev = map(last, gsrev)
    τsrev = map(x->getindex(x, 2), gsrev)
    wsrev = normalize(τsrev, 1)
    graph_pairs_rev = as_pairs.(first.(gsrev))

    cmrev = keyedreduce(+, graph_pairs_rev, wsrev)
    cmrev = sort(cmrev; byvalue=true, rev=true)

end


 
println("# cpdags: $( n ≤ length(ncpdags) ? ncpdags[n] : "NA") (true), $(length(cm)) (sampled) " )
reversible_too && println("# cpdags: $( n ≤ length(ncpdags) ? ncpdags[n] : "NA") (true), $(length(cmrev)) (sampled) " )

println("prob: ", n ≤ length(ncpdags) ? round(1/(ncpdags[n]), sigdigits=3) : "NA"," (true), ", extrema(values(cm)), "(estimates) " )

if n == κ - 1 && n ≤ length(ncpdags) 
    println("rmse ", norm(values(cm) .- 1/(ncpdags[n])))
else

    println("extrema", extrema(values(cm)))
    reversible_too && println("extrema", extrema(values(cmrev)))
end
println("mean repetitions ", mean(rle(hash.(graph_pairs))[2]))


function figure(;autocor=false) 
    
    fig = Figure()
    ax1 = fig[1,1] = Axis(fig)
    if autocor 
        ax2 = fig[2,1] = Axis(fig)
    end
    tim = cumsum(wsnonrev)
    stairs!(ax1, tim, hsnonrev, step=:post)
    autocor && lines!(ax2, autocor(hsnonrev))

    if @isdefined hsrev 
        stairs!(ax1, cumsum(wsrev), hsrev, color=:orange, step=:post)
        autocor && lines!(ax2, autocor(hsrev))
    end

    ylims!(ax1, 0, n*(κ)÷2)
    #ylims!(ax1, 0, 250)

   # xlims!(ax1, tim[max(1, length(hsnonrev) - 1000)], tim[end])
    ax1.yzoomlock = true
    fig
end
# using GLMakie
@isdefined(Figure) && (fig = figure(); display(fig))

if score != UniformScore()
    println("Maximum a posteriory estimate: ", argmax(cm))
    println("True CPDAG:                    ", true_cpdag)
end


if !uniform
    using GraphMakie, GLMakie
    function graphdiff(g1, g2)
        @assert nv(g1) == nv(g2)
        fig, ax, pl = graphplot(g1; edge_color=(:darkorange, 0.3),  kwargs_pdag_graphmakie(g1)...)
        graphplot!(ax, g2; edge_color=(:blue, 0.2), layout=pl[:node_pos][], kwargs_pdag_graphmakie(g2)...)
        fig
    end
    fig2 = graphdiff(digraph(true_cpdag,n), digraph(first(keys(cm)),n))
end

cm

fig2 = lines(log.(values(sort(keyedreduce(+, hs, ws)))), color=:blue)
lines!(log.(values(sort(keyedreduce(+, hsrev, wsrev)))), color=:orange)
fig2
fig
