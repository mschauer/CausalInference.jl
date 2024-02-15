using LogarithmicNumbers
# Valid balancing functions

metropolis_balance(t) = min(one(t), t)
sqrt_balance(t) = sqrt(t)
barker_balance(t) = t/(1+t) # softmin

struct UniformScore
end
#import CausalInference.local_score
local_score(::UniformScore, _, _) = 0.0


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


"""
    count_moves(g, κ, balance, prior, score, coldness, dir=:both)

Return 
"""
function count_moves(g, κ, balance, prior, score, coldness, total, dir=:both)
    s1 = s2 = Δscorevalue1 = Δscorevalue2 = 0.0
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
                        if valid
                            PAy = parents(g, y)
                            Δscorevalue = Δscoreinsert(score, NAyxT ∪ PAy, x, y, T)
                            s = balance(prior(total, total+1)*exp(coldness*Δscorevalue))
                        else 
                            s = 0.0
                        end
                        @assert s >= 0
                        if valid && rand() > s1/(s1 + s) # sequentially draw sample
                            x1, y1 = x, y
                            T1 = T
                            Δscorevalue1 = Δscorevalue
                        end
                        s1 = s1 + s
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
                    if valid
                        PAy = parents(g, y)
                        PAy⁻ = setdiff(PAy, x)
                        Δscorevalue = Δscoredelete(score, NAyx_H ∪ PAy⁻, x, y, H)
                        s = balance(prior(total, total-1)*exp(coldness*Δscorevalue))
                    else 
                        s = 0.0
                    end
                    @assert s >= 0
                    if valid && rand() > s2/(s2 + s)
                        x2, y2 = x, y
                        H2 = H
                        Δscorevalue2 = Δscorevalue
                    end
                    s2 = s2 + s
                end
            end 
        end
    end
    s1, s2, Δscorevalue1, Δscorevalue2, (x1, y1, T1), (x2, y2, H2)
end

function count_moves_new(g, κ, balance, prior, score, coldness, total, dir=:both)
    s1 = s2 = ULogarithmic(0.0)
    Δscorevalue = Δscorevalue1 = Δscorevalue2 = 0.0
    x1 = y1 = x2 = y2 = 0
    T1 = Int[]
    H2 = Int[] 
    noinsert = nodelete = true
    for y in vertices(g)
        noinsert = (dir == :down)
        nodelete = (dir == :up)
        length(neighbors_adjacent(g, y)) >= κ && (noinsert = true)
        total == 0 && (nodelete = true)
        
        #total == nv(g)*κ÷2-2 && (noinsert = true)
        dir == :up && noinsert && continue
        !noinsert && (semidirected = precompute_semidirected(g, y))
        PAy = parents(g, y)
        for x in vertices(g)
            x == y && continue
            if !noinsert 
                length(neighbors_adjacent(g, x)) >= κ && continue
                insit = InsertIterator(g, x, y, semidirected)
                for T in insit
                    NAyxT = union(adj_neighbors(g, x, y), T)
                    # maybe we could do Δscoreinsert(score, g, x, y, T)
                    # to hide complexity
                    # or just Δscoreinsert(score, g, op)
                    # and op contains all necessary stuff e.g. NAyxT and so on
                    Δscorevalue = Δscoreinsert(score, NAyxT ∪ PAy, x, y, T)
                    s = balance(prior(total, total+1)*exp(ULogarithmic, coldness*Δscorevalue))
                    if rand() > 1/(1 + s/s1) # sequentially draw sample
                        x1, y1 = x, y
                        T1 = T
                        Δscorevalue1 = Δscorevalue
                    end
                    s1 = s1 + s
                end
            end
            if !nodelete
                delit = DeleteIterator(g, x, y)
                for H in delit
                    PAy⁻ = setdiff(PAy, x)
                    # I would prefer Δscoredelete(score, g, x, y, H) as above
                    NAyx_H = setdiff(adj_neighbors(g, x, y), H)
                    Δscorevalue = Δscoredelete(score, NAyx_H ∪ PAy⁻, x, y, H)
                    s = balance(prior(total, total-1)*exp(ULogarithmic, coldness*Δscorevalue))
                    if rand() >  1/(1 + s/s2) 
                        x2, y2 = x, y
                        H2 = H
                        Δscorevalue2 = Δscorevalue
                    end
                    s2 = s2 + s
                end
            end
        end
    end
    #@show x1 y1 x2 y2 total s1 s2 noinsert nodelete Δscorevalue
    noinsert || @assert x1 != 0 && y1 != 0
    nodelete || @assert x2 != 0 && y2 != 0
   
    s1, s2, Δscorevalue1, Δscorevalue2, (x1, y1, T1), (x2, y2, H2)
end

"""
    causalzigzag(n, G = (DiGraph(n), 0); balance = metropolis_balance, prior = (_,_)->1.0, score=UniformScore(),
                        coldness = 1.0, σ = 0.0, ρ = 1.0, naive=false,
                        κ = min(n - 1, 10), iterations=10, verbose=false, save=true)

Run the causal zigzag algorithm starting in a cpdag `(G, t)` with `t` oriented or unoriented edges,
the balance function `balance ∈ {metropolis_balance, barker_balance, sqrt}`, `score` function (see `ges` algorithm)
coldness parameter for iterations. `σ = 1.0, ρ = 0.0` gives purely diffusive behaviour, `σ = 0.0, ρ = 1.0` gives Zig-Zag behaviour.

Returns a vector of tuples with information, each containing a graph, spent time, current direction, number of edges and the score.
"""
function causalzigzag(n, G = (DiGraph(n), 0); balance = metropolis_balance, prior = (_,_)->1.0, score=UniformScore(),
                        coldness = 1.0, σ = 0.0, ρ = 1.0, naive=false,
                        κ = min(n - 1, 10), iterations=10, verbose=false, save=true)
    g, total = G
    if κ >= n 
        κ = n - 1
        @warn "Truncate κ to $κ"
    end
    gs = Vector{Tuple{typeof(g),Float64,Int,Int,Float64}}()
    dir = 1
    global traversals = 0
    global tempty = 0.0
    τ = 0.0
    secs = 0.0
    emax = n*κ÷2
    scorevalue = 0.0
    @showprogress for iter in 1:iterations
        total_old = total
        if isodd(traversals) && total == 0 # count number of traversal from empty to full
            traversals += 1
        elseif iseven(traversals) && total == emax
            traversals += 1
        end
        
        Δscorevalue = Δscorevalue1 = Δscorevalue2 = 0.0
        if !naive 
            if score isa UniformScore
                s1, s2, up1, down1 = count_moves_uniform(g, κ)
                total < emax && (s1 *= balance(prior(total, total+1)))
                total > 0 && (s2 *= balance(prior(total, total-1)))
                
            else 
                s1, s2, Δscorevalue1, Δscorevalue2, up1, down1 = count_moves_new(g, κ, balance, prior, score, coldness, total)
            end
        else
            s1, s2, Δscorevalue1, Δscorevalue2, up1, down1 = count_moves(g, κ, balance, prior, score, coldness, total)
        end
        @label flipped

        stuck = false
        τ = 0.0
        total_old = total
        dir_old = dir

        λbar = max(dir*(-s1 + s2), 0.0)
        λrw = (s1 + s2) 
        λup = (s1)   
        λ = dir == 1 ? (s1) : (s2)
        
        local x, y
        while true 

            Δτ = randexp()/(ρ*λ)
            Δτrw = randexp()/(σ*λrw)
            Δτflip = randexp()/(ρ*λbar)
            τ += float(min(Δτ, Δτflip, Δτrw))
            if min(Δτ, Δτflip, Δτrw) > 1.0e10
                @warn "Frozen in the cold at iteration $iter" # $Δτ, $Δτrw, $Δτflip"
                stuck = true
                @goto flip
            end
            if Δτflip < Δτ &&  Δτflip < Δτrw 
                @goto flip
            end
            up = rand() < λup/λrw           

            # @show λup λrw  Δτ Δτrw
         
            if  (Δτ <= Δτrw  && dir == 1) || (Δτ > Δτrw  && up)

                
                x, y, T = up1
                @assert x != y
                total == 0 && (tempty += τ)
                save && push!(gs, (g, τ, dir, total, scorevalue))
                Δscorevalue = Δscorevalue1
                scorevalue += Δscorevalue
                total += 1
                secs += @elapsed begin
                    if !naive
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
                total == 0 && (tempty += τ)
                save && push!(gs, (g, τ, dir, total, scorevalue))
                Δscorevalue = Δscorevalue2
                scorevalue += Δscorevalue
                total -= 1
                secs += @elapsed begin
                    if !naive
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
            x = y = 0
            dir *= -1
            total == 0 && (tempty += τ)
            save && push!(gs, (g, τ, dir, total, scorevalue))
            stuck && break
            @goto flipped           
        end # break
        verbose && println(total_old, dir_old == 1 ? "↑" : "↓", total,  " $x => $y ", round(τ, sigdigits=5), " ", round(Δscorevalue, sigdigits=5), " ", round(scorevalue, sigdigits=5))
        #verbose && println("\t", vpairs(g))
        stuck && break
    end
    println("time moves $secs")
    println("nr. traversals $traversals")
    println("time empty $tempty")

    gs
end


function unzipgs(gs)
    graphs = first.(gs)
    graph_pairs = vpairs.(graphs)
    scs = map(last, gs)
    hs = map(x->getindex(x, 4), gs)
    τs = map(x->getindex(x, 2), gs)
    ws = normalize(τs, 1)
    ts = cumsum(ws)
    (;graphs, graph_pairs, hs, τs, ws, ts, scs)
end 

const randcpdag = causalzigzag
