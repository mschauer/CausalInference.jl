using LogarithmicNumbers


"""
    count_dag_moves(g, κ, balance, prior, score, coldness, dir=:both)

Return 
""" 
function count_dag_moves(g, κ, balance, prior, score, coldness, total, dir=:both)
    s1 = s2 = ULogarithmic(0.0)
    Δscorevalue = Δscorevalue1 = Δscorevalue2 = 0.0
    x1 = y1 = x2 = y2 = 0
    noinsert = nodelete = true
    for y in vertices(g)
        noinsert = (dir == :down)
        nodelete = (dir == :up)
        length(neighbors_adjacent(g, y)) >= κ && (noinsert = true)
        total == 0 && (nodelete = true)
        dir == :up && noinsert && continue
        for x in vertices(g)
            x == y && continue
            if !noinsert 
                length(neighbors_adjacent(g, x)) >= κ && continue
                if !isadjacent(g, x, y) && !has_path(g, y, x)  
                    Δscorevalue = Δscore(score, inneighbors(g, y), x, y)
                    s = balance(prior(total, total+1)*exp(ULogarithmic, coldness*Δscorevalue))
                    if rand() > 1/(1 + s/s1) # sequentially draw sample
                        x1, y1 = x, y
                        Δscorevalue1 = Δscorevalue
                    end
                    s1 = s1 + s
                end
            end
        end
        for x in inneighbors(g, y)
            if !nodelete
                Δscorevalue = -Δscore(score, setdiff(inneighbors(g, y), x), x, y)
                s = balance(prior(total, total-1)*exp(ULogarithmic, coldness*Δscorevalue))
                if rand() >  1/(1 + s/s2) 
                    x2, y2 = x, y
                    Δscorevalue2 = Δscorevalue
                end
                s2 = s2 + s     
            end
        end
    end
    #@show x1 y1 x2 y2 total s1 s2 noinsert nodelete Δscorevalue
    noinsert || @assert x1 != 0 && y1 != 0
    nodelete || @assert x2 != 0 && y2 != 0
   
    s1, s2, Δscorevalue1, Δscorevalue2, (x1, y1), (x2, y2)
end

function count_dag_turn_moves(g, κ, balance, prior, score, coldness, total)
    s0 = ULogarithmic(0.0)
    Δscorevalue0 = 0.0
    x0 = y0 = 0
    for y in vertices(g)
        for x in inneighbors(g, y) # remove x -> y, add y -> x
            if has_a_path(g, setdiff(outneighbors(g, x), y), y) 
                continue
            end
            Δscorevalue = -Δscore(score, setdiff(inneighbors(g, y), x), x, y)
            Δscorevalue += Δscore(score, inneighbors(g, x), y, x)
                  
            s = balance(prior(total, total-1)*exp(ULogarithmic, coldness*Δscorevalue))
            if rand() >  1/(1 + s/s0) 
                x0, y0 = x, y
                Δscorevalue0 = Δscorevalue
            end
            s0 = s0 + s     
        end
    end
    s0, Δscorevalue0, (x0, y0)
end

"""
    dagzigzag(n, G = DiGraph(n); balance = metropolis_balance, prior = (_,_)->1.0, score=UniformScore(),
                        coldness = 1.0, σ = 0.0, ρ = 1.0, 
                        κ = min(n - 1, 10), iterations=10, verbose=false, save=true)

Run the causal zigzag algorithm starting in a dag `G` 
the balance function `balance ∈ {metropolis_balance, barker_balance, sqrt}`, `score` function (see `ges` algorithm)
coldness parameter for iterations. `σ = 1.0, ρ = 0.0` gives purely diffusive behaviour, `σ = 0.0, ρ = 1.0` gives Zig-Zag behaviour.

Returns a vector of tuples with information, each containing a graph, spent time, current direction, number of edges and the score.
"""
function dagzigzag(n, G = DiGraph(n); balance = metropolis_balance, prior = (_,_)->1.0, score=UniformScore(),
                        coldness = 1.0, σ = 0.0, ρ = 1.0,
                        κ = min(n - 1, 10), iterations=10, verbose=false, save=true)
    g = G
    if κ >= n 
        κ = n - 1
        @warn "Truncate κ to $κ"
    end
    gs = Vector{Tuple{typeof(g),Float64,Int,Int,Float64}}()
    dir = 1
    global traversals = 0
    global tempty = 0.0
    τ = 0.0
    emax = n*κ÷2
    scorevalue = 0.0
    @showprogress for iter in 1:iterations
        total = total_old = ne(g)
       
        if isodd(traversals) && total == 0 # count number of traversal from empty to full
            traversals += 1
        elseif iseven(traversals) && total == emax
            traversals += 1
        end
        
        Δscorevalue =  Δscorevalue0 = Δscorevalue1 = Δscorevalue2 = 0.0
    
        if false #fixme score isa UniformScore
            s1, s2, up1, down1 = count_dag_moves_uniform(g, κ)
            total < emax && (s1 *= balance(prior(total, total+1)))
            total > 0 && (s2 *= balance(prior(total, total-1)))
            
        else 
            s1, s2, Δscorevalue1, Δscorevalue2, up1, down1 = count_dag_moves(g, κ, balance, prior, score, coldness, total)
            s0, Δscorevalue0, turn0 = count_dag_turn_moves(g, κ, balance, prior, score, coldness, total)   
        end
        @label flipped
        total = ne(g)
       
        stuck = false
        τ = 0.0
        dir_old = dir
        λbar = max(dir*(-s1 + s2), 0.0)
        λrw = s1 + s2
        λup = s1 
        λturn = s0  
        λ = dir == 1 ? s1 : s2
        
        local x, y
        while true 

            Δτ = randexp()/(ρ*λ)
            Δτturn = randexp()/(λturn)
            Δτrw = randexp()/(σ*λrw)
            Δτflip = randexp()/(ρ*λbar)
            τ += float(min(Δτ, Δτturn, Δτflip, Δτrw))
            if min(Δτ, Δτturn, Δτflip, Δτrw) > 1.0e10
                @warn "Frozen in the cold at iteration $iter" # $Δτ, $Δτrw, $Δτflip"
                stuck = true
                @goto flip
            end
            if Δτflip < Δτ &&  Δτflip < Δτrw && Δτflip < Δτturn
                @goto flip
            end
            if Δτturn < Δτ &&  Δτturn < Δτrw 
                x, y = turn0
                @assert x != y
                @assert total > 0
                save && push!(gs, (g, τ, dir, total, scorevalue))
                g = copy(g) # copy on save?
                rem_edge!(g, x, y)
                add_edge!(g, y, x)
                Δscorevalue = Δscorevalue0
                scorevalue += Δscorevalue
                break
            end

            up = rand() < λup/λrw           

            # @show λup λrw  Δτ Δτrw
            if  (Δτ <= Δτrw  && dir == 1) || (Δτ > Δτrw  && up)

                
                x, y = up1
                @assert x != y
                total == 0 && (tempty += τ)
                save && push!(gs, (g, τ, dir, total, scorevalue))
                g = copy(g) # copy on save?
                add_edge!(g, x, y)
                Δscorevalue = Δscorevalue1
                scorevalue += Δscorevalue
                total += 1

                break
            else
                x, y = down1

                @assert x != y
                total == 0 && (tempty += τ)
                save && push!(gs, (g, τ, dir, total, scorevalue))

                g = copy(g) # copy on save?
                rem_edge!(g, x, y)
                Δscorevalue = Δscorevalue2
                scorevalue += Δscorevalue
                total -= 1

                break
            end
            @label flip
            x = y = 0
            dir *= -1
            total == 0 && (tempty += τ)
            save && push!(gs, (g, τ, dir, total, scorevalue))
            stuck && break
            @goto flipped
        end # break here
        verbose && println(total_old, dir_old == 1 ? "↑" : "↓", total,  " $x => $y ", round(τ, sigdigits=5), " ", round(Δscorevalue, sigdigits=5), " ", round(scorevalue, sigdigits=5))
        #verbose && println("\t", vpairs(g))
        stuck && break
    end # break here
    println("nr. traversals $traversals")
    println("time empty $tempty")

    gs
end

