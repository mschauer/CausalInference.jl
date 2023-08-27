using Revise
using LinearAlgebra
using CausalInference
using CausalInference.Graphs
using StatsBase
using ProgressMeter
const as_pairs = vpairs
using Random
Random.seed!(2)

include("mcs.jl")

using CausalInference: isadjacent, tails_and_adj_neighbors, adj_neighbors, isclique, Insert!, Delete!
using CausalInference: isundirected, parents, meek_rule1, meek_rule2, meek_rule3, meek_rule4, children, 
    neighbors_adjacent, neighbors_undirected
using CausalInference.Combinatorics
const adjacents = neighbors_adjacent

"""
    keyedreduce(op, key::AbstractVector{T}, a, init=0.0) where T

Similar to `countmap` returning a dictionary mapping unique key in `key` to the 
reduction the given collection itr with the given binary operator `op`.

```
julia> keyedreduce(+, [:a, :b, :a], [7, 3, 2])
Dict{Symbol, Float64} with 2 entries:
  :a => 9.0
  :b => 3.0
```
"""
function keyedreduce(op, key::AbstractVector{T}, a, init=0.0) where T
    cm = Dict{T,typeof(init)}()
    for (v, τ) in zip(key, a)
        index = Base.ht_keyindex2!(cm, v)
        if index > 0
            @inbounds cm.vals[index] = op(cm.vals[index] , τ)
        else
            @inbounds Base._setindex!(cm, op(init, τ), v, -index)
        end
    end
    return cm
end

"""
    ne_total(g)

Number of directed or undirected edges in a PDAG represented as DiGraph.
"""
function ne_total(g)
    s = 0
    for x in vertices(g)
        for y in vertices(g)
            y < x && continue
            s += isadjacent(g, x, y)
        end
    end
    s
end


"""
    ne_undirected(g)

Number of undirected edges in a PDAG represented as DiGraph.
"""
function ne_undirected(g)
    s = 0
    for x in vertices(g)
        for y in vertices(g)
            y < x && continue
            s += isundirected(g, x, y)
        end
    end
    s
end

isblocked(g, x, y, nodesRemoved) = !has_a_path(g, [x], y, nodesRemoved)

ncpdags = [1, 2, 11, 185, 8782, 1067825, 312510571, 212133402500, 326266056291213, 1118902054495975181, 8455790399687227104576, 139537050182278289405732939, 4991058955493997577840793161279]

"""
    nup(g, total)

Number of directed edges that can be add to a PDAG `g` with `total` number of edges.
"""
nup(g, total) = nup(nv(g), total)
nup(n::Int, total) = 2*(n*(n - 1)÷2 - total) # twice as number of undirected edges

"""
    ndown(g, total)

Number of edges that can be removed from a `pdag` counting undirected edges twice.
"""
ndown(g, total) = ne(g)



function exact(g, κ, dir=:both) 
    s1, s2, _ = exact2(g, κ, dir)
    -s1 + s2
end

"""
    exact2(g, κ, dir=:both)

Return 
"""
function exact2(g, κ, dir=:both)
    s1 = s2 = 0.0
    x1 = y1 = x2 = y2 = 0
    T1 = Int[]
    H2 = Int[] 
    for x in vertices(g)
        for y in vertices(g)
            x == y && continue
            if !isadjacent(g, x, y) && dir != :down
                if length(neighbors_adjacent(g, x)) < κ && length(neighbors_adjacent(g, y)) < κ 
                    Tyx, NAyx = tails_and_adj_neighbors(g, x, y)
                    @assert length(Tyx) < 127
                    for i in 0:UInt128(2)^length(Tyx) - 1
                        T = Tyx[digits(Bool, i, base=2, pad=length(Tyx))]
                        NAyxT = CausalInference.sorted_union_(NAyx, T)
                        valid = (isclique(g, NAyxT) && isblocked(g, y, x, NAyxT))
                        if valid && rand() > s1/(s1 + 1) # sequentially draw sample
                            x1, y1 = x, y
                            T1 = T
                        end
                        s1 = s1 + valid
                    end
                end
            elseif has_edge(g, x, y) && dir != :up
                Hyx = adj_neighbors(g, x, y)
                @assert length(Hyx) < 127
                for i in 0:UInt128(2)^length(Hyx) - 1
                    mask = digits(Bool, i, base=2, pad=length(Hyx))
                    H = Hyx[mask] 
                    NAyx_H = Hyx[map(~, mask)] 
                    valid = isclique(g, NAyx_H)
                    if valid && rand() > s2/(s2 + 1)
                        x2, y2 = x, y
                        H2 = H
                    end
                    s2 = s2 + valid 
                end
            end 
        end
    end
    s1, s2, (x1, y1, T1), (x2, y2, H2)
end

function countcliques(g)
    n = nv(g)
    preceding, _ = countmcs(g)
    # maybe use BigInt at some point
    cnt = 1 # don't forget "empty" clique
    for i = 1:n
        cnt += 2^preceding[i]
    end
    return cnt
end

function sampleclique(g, r)
    n = nv(g)
    preceding, invmcsorder = countmcs(g)
    cnt = 1 # don't forget "empty" clique
    r <= cnt && return Vector{Int64}()
    for i = 1:n
        cnt += 2^preceding[i]
        if r <= cnt
            p = Vector{Int64}()
            for j in inneighbors(g, i)
                invmcsorder[j] < invmcsorder[i] && push!(p, j) 
            end
            subset = rand(0:2^preceding[i]-1)
            ret = Vector{Int64}()
            for d = 1:length(p)
                if ((1 << d) & subset) > 0
                    push!(ret, p[d])
                end
            end
            return ret
        end
    end
end

function exactup(g, κ)
    n = nv(g)
    ttcnt = 0
    cnts = Vector{Int64}()
    cntids = Vector{Tuple{Int64, Int64}}()
    need = Vector{Set{Int64}}()
    pot = Vector{Vector{Int64}}()
    for y in vertices(g)
        length(neighbors_adjacent(g, y)) == κ && continue
        blocked = falses(n)
        nu = neighbors_undirected(g, y)
        push!(nu, y)
        for v in nu
            blocked[v] = true
        end
        visfrom = Vector{BitVector}()
        for z in nu
            # => find all vertices reachable by semidirected path (with other neighbors of y and y itself blocked) from z and save somewhere
            vis = falses(n)
            q = Vector{Int64}()
            push!(q, z)
            vis[z] = true
            # can optimize further by not looking at undir after first dir edge
            while !isempty(q)
                w = popfirst!(q)
                for v in outneighbors(g, w)
                    vis[v] && continue
                    blocked[v] && continue
                    push!(q, v)
                    vis[v] = true
                end
            end
            push!(visfrom, vis)
        end
        for x in vertices(g)
            x == y && continue
            isadjacent(g, x, y) && continue
            last(visfrom)[x] && continue
            length(neighbors_adjacent(g, x)) == κ && continue
            needtotake = Set(adj_neighbors(g, x, y))
            # => get all other undirected neighbors we have to take to close semidirected paths
            for i = 1:length(nu)
                if visfrom[i][x]
                    push!(needtotake, nu[i])
                end
            end
            !isclique(g, needtotake) && continue
            # => keep all other undirected neighbors which are connected to all the must have ones
            potential = Vector{Int64}()
            for v in nu
                v == y && continue
                v in needtotake && continue
                fullyconnected = true
                for w in needtotake
                    if !isadjacent(g, v, w)
                        fullyconnected = false
                        break
                    end
                end
                fullyconnected && push!(potential, v)
            end
            # => count number of cliques in chordal graph
            sg, _ = induced_subgraph(g, potential)
            cnt = countcliques(sg)
            ttcnt += cnt
            push!(cnts, ttcnt)
            push!(cntids, (x,y))
            push!(need, needtotake)
            push!(pot, potential)
            # println(x, " ", y, " ", ttcnt) 
        end
    end
    # store counts for each pair x,y and sample one of them
    # then sample a clique randomly in chordal graph -> by same procedure find# highest index node and then take random subset of prev visited neighbors
    ttcnt == 0 && return ttcnt, (0, 0, Int[]) 
    r = rand(1:ttcnt)
    for i = 1:length(cnts)
        if r <= cnts[i]
            (x, y) = cntids[i]
            cnt = cnts[i]
            i > 1 && (cnt -= cnts[i-1])
            r2 = rand(1:cnt)
            #println(NAyx)
            sg, _ = induced_subgraph(g, pot[i])
            cl = sampleclique(sg, r2)
            append!(cl, need[i])
            #println(cl)
            return ttcnt, (x, y, cl) 
        end
    end
end

function exactdown(g)
    ttcnt = 0
    cnts = Vector{Int64}()
    cntids = Vector{Tuple{Int64, Int64}}()
    for x in vertices(g)
        for y in [neighbors_undirected(g, x); children(g, x)]
            NAyx = adj_neighbors(g, x, y)
            if length(NAyx) == 0
                ttcnt += 1
            elseif length(NAyx) == 1
                ttcnt += 2
            else
                sg, _ = induced_subgraph(g, NAyx)
                cnt = countcliques(sg)
                ttcnt += cnt
            end
            push!(cnts, ttcnt)
            push!(cntids, (x,y))
        end
    end 
    ttcnt == 0 && return 0, (0, 0, Int[])
    # sample clique
    r = rand(1:ttcnt)
    for i = 1:length(cnts)
        if r <= cnts[i]
            (x, y) = cntids[i]
            cnt = cnts[i]
            i > 1 && (cnt -= cnts[i-1])
            r2 = rand(1:cnt)
            NAyx = adj_neighbors(g, x, y)
            #println(NAyx)
            sg, _ = induced_subgraph(g, NAyx)
            cl = sampleclique(sg, r2)
            #println(cl)
            return ttcnt, (x, y, setdiff(NAyx, cl)) 
        end
    end
end

function exact3(g, κ)
    s1, (x1, y1, T1) = exactup(g, κ)
    s2, (x2, y2, H2) = exactdown(g)
    return s1, s2, (x1, y1, T1), (x2, y2, H2)
end

function randcpdag(n, G = (DiGraph(n), 0); σ = 0.0, ρ = 1.0, wien=true,
                        κ = min(n - 1, 10), iterations=10, verbose=false)
    g, total = G
    if κ >= n 
        κ = n - 1
        @warn "Truncate κ to $κ"
    end
    total = 0
    gs = Vector{Tuple{typeof(g),Float64,Int,Int}}()
    dir = 1
    traversals = 0
    τ = 0.0
    secs = 0.0
    olddown = 0.0
    newdown = 0.0
    oldup = 0.0
    newup = 0.0
    @showprogress for iter in 1:iterations
        τ = 0.0
        total_old = total
        dir_old = dir

        if isodd(traversals) && total == 0 # count number of traversal from empty to full
            traversals += 1
        elseif iseven(traversals) && total == nv(g)*(nv(g) - 1)÷2 
            traversals += 1
        end
        olddown += @elapsed _, ss2, _, _ = exact2(g, κ, :down)
        newdown += @elapsed t2, _ = exactdown(g) 
        @assert ss2 == t2
       
        oldup += @elapsed ss1, _, _, _ = exact2(g, κ, :up)
        newup += @elapsed t1, _ = exactup(g, κ) 
        @assert ss1 == t1

        
        s1, s2, up1, down1 = exact2(g, κ)
        λbar = max(dir*(-s1 + s2), 0.0)
        λrw = (s1 + s2) 
        λup = s1   
        λ = dir == 1 ? abs(s1) : abs(s2)

        while true 

            Δτ = randexp()/(ρ*λ)
            Δτrw = randexp()/(σ*λrw)
            Δτflip = randexp()/(ρ*λbar)
            τ += min(Δτ, Δτflip, Δτrw)
            if isinf(τ) 
                error("$Δτ, $Δτrw, $Δτflip")
            end
            if Δτflip < Δτ &&  Δτflip < Δτrw 
                @goto flip
            end
            up = rand() < λup/λrw           

     
         
            if  (Δτ < Δτrw  && dir == 1) || (Δτ > Δτrw  && up)

                local x, y
                x, y, T = up1
                @assert x != y
                push!(gs, (g, τ, dir, total))
                total += 1
                secs += @elapsed begin
                    if wien
                        g = next_CPDAG(g, :up, x, y, T)
                    else
                        g = copy(g)
                        Insert!(g, x, y, T)
                        vskel!(g)
                        meek_rules!(g)
                    end
                end
                break
            else
                x, y, H = down1
                @assert x != y
                push!(gs, (g, τ, dir, total))
                total -= 1
                secs += @elapsed begin
                    if wien
    		        g = next_CPDAG(g, :down, x, y, H)
                    else
                        g = copy(g)
                        Delete!(g, x, y, H)
                        vskel!(g)
                        meek_rules!(g)
                    end
                end
                break
            end
            @label flip
            dir *= -1
            push!(gs, (g, τ, dir, total))
            break            
        end # break
        verbose && println(total_old, dir_old == 1 ? "↑" : "↓", total, " ", round(τ, digits=8))
        verbose && println("\t", vpairs(g))
    end
    println("time moves $secs")
    println("nr. traversals $traversals")
    println("cmpdown $olddown $newdown")
    println("cmpup $oldup $newup")
    gs
end

iterations = 10_000; verbose = false
n = 50 # vertices
κ = n - 1 # max degree
reversible_too = false # do baseline 
#iterations = 50; verbose = true
burnin = iterations÷2

gs = @time randcpdag(n; ρ=1.0, σ=0.0, wien=true, κ, iterations, verbose)[burnin:end]

graphs = first.(gs)
graph_pairs = as_pairs.(graphs)
hs = hsnonrev = map(last, gs)
undirecteds = ne.(graphs) - hs


i = rand(eachindex(hs))
g, h = graphs[i], hs[i]

τs = map(x->getindex(x, 2), gs)
ws = wsnonrev = normalize(τs, 1)

println("Average nr. of undirected edges: " , sum(undirecteds .* ws))



if reversible_too
    gsrev = @time randcpdag(n; ρ=0.0, σ=1.0, κ, iterations, verbose)[burnin:end]
    hsrev = map(last, gsrev)
    τsrev = map(x->getindex(x, 2), gsrev)
    wsrev = normalize(τsrev, 1)
end

@show sum(τs)
cm = keyedreduce(+, graph_pairs, ws)
dirs = map(x->getindex(x, 3), gs)
cm2 = keyedreduce(+, graph_pairs, ws .* dirs)

println("# cpdags: $( n ≤ length(ncpdags) ? ncpdags[n] : "NA") (true), $(length(cm)) (sampled) " )
println("prob: ", n ≤ length(ncpdags) ? round(1/(ncpdags[n]), sigdigits=3) : "NA"," (true), ", extrema(values(cm)), "(estimates) " )

n ≤ length(ncpdags) && println("rmse ", norm(values(cm) .- 1/(ncpdags[n])))
println("mean repetitions ", mean(rle(hash.(graph_pairs))[2]))


function figure() 
    
    fig = Figure()
    ax1 = fig[1,1] = Axis(fig)
    ax2 = fig[2,1] = Axis(fig)
    tim = cumsum(wsnonrev)
    stairs!(ax1, tim, hsnonrev, step=:post)
    lines!(ax2, autocor(hsnonrev))

    if @isdefined hsrev 
        stairs!(ax1, cumsum(wsrev), hsrev, color=:orange, step=:post)
        lines!(ax2, autocor(hsrev))
    end

    ylims!(ax1, 0, n*(n-1)÷2)
    xlims!(ax1, tim[max(1, length(hsnonrev) - 1000)], tim[end])
    ax1.yzoomlock = true
    fig
end
# using GLMakie
@isdefined(Figure) && (fig = figure(); display(fig))

cm

