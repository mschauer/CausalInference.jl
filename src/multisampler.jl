struct Sample
    g::DiGraph
    τ::Float64
    dir::Int8
    total::Int32 
    scoreval::Float64
end

struct Action
    τ::Float64 
    apply::Function
    args::Tuple{Vararg{Any}}
end

function expcoldness(τ, k=0.0005)
    return exp(k*τ)
end

function Dexpcoldness(τ, k=0.0005)
    return k*exp(k*τ)
end

function init(_, _, nextτ, g, dir, total, scoreval)
    return Sample(g, nextτ, dir, total, scoreval) 
end

function applyup(samplers, i, nextτ, x, y, T, Δscoreval)
    prevsample = last(samplers[i])
    g = next_CPDAG(prevsample.g, :up, x, y, T)
    return Sample(g, nextτ, prevsample.dir, prevsample.total+1, prevsample.scoreval + Δscoreval)
end

function applydown(samplers, i, nextτ, x, y, H, Δscoreval)
    prevsample = last(samplers[i])
    g = next_CPDAG(prevsample.g, :down, x, y, H)
    return Sample(g, nextτ, prevsample.dir, prevsample.total-1, prevsample.scoreval + Δscoreval) 
end

function applyflip(samplers, i, nextτ)
    prevsample = last(samplers[i])
    return Sample(prevsample.g, nextτ, -1*prevsample.dir, prevsample.total, prevsample.scoreval)
end

function applycopy(samplers, _, nextτ, j)
    copysample = last(samplers[j])
    return Sample(copysample.g, nextτ, copysample.dir, copysample.total, copysample.scoreval)
end

# for starters without turn move
function sampleaction(samplers, i, M, balance, prior, score, σ, ρ, κ, coldness, Dcoldness) 
    # preprocess 
    prevsample = last(samplers[i])
    sup, sdown, Δscorevalup, Δscorevaldown, argsup, argsdown = count_moves_new(prevsample.g, κ, balance, prior, score, coldness(prevsample.τ), prevsample.total)

    # propose moves
    λdir = prevsample.dir == 1 ? sup : sdown 
    λupdown = sup + sdown 
    λflip = max(prevsample.dir*(-sup + sdown), 0.0)
    λterm = -exp(ULogarithmic, 0.0)*Dcoldness(prevsample.τ) * prevsample.scoreval # TODO: prior
    @assert λterm >= 0
    Δτdir = randexp()/(ρ*λdir)
    Δτupdown = randexp()/(σ*λupdown)
    Δτflip = randexp()/(ρ*λflip)
    Δτterm = randexp()/(λterm) 

    Δτmin = min(Δτdir, Δτupdown, Δτflip, Δτterm)

    if Δτdir == Δτmin 
        if prevsample.dir == 1
            return Action(prevsample.τ + Δτdir, applyup, (argsup..., Δscorevalup))
        else 
            return Action(prevsample.τ + Δτdir, applydown, (argsdown..., Δscorevaldown))
        end
    end

    if Δτupdown == Δτmin
        λup = sup
        if rand() < λup/λupdown
            return Action(prevsample.τ + Δτupdown, applyup, (argsup..., Δscorevalup))
        else 
            return Action(prevsample.τ + Δτupdown, applydown, (argsdown..., Δscorevaldown))
        end
    end
    
    if Δτflip == Δτmin 
        return Action(prevsample.τ + Δτflip, applyflip, ())
    end

    if Δτterm == Δτmin 
        return Action(prevsample.τ + Δτterm, applycopy, (rand(1:M),)) 
    end

    @assert false
end

function save!(as, a)
    if length(as) == 0
        push!(as, a)
    else
        as[end] = a
    end
end

# remark: chose κ = n-1 as default
function multisampler(n, G = (DiGraph(n), 0); M = 10, balance = metropolis_balance, prior = (_,_) -> 1.0, score=UniformScore(), σ = 0.0, ρ = 1.0, κ = n - 1, iterations = min(3*n^2, 50000), schedule=(expcoldness, Dexpcoldness)) #, verbose = false, save = true)
    if κ >= n 
        κ = n - 1
        @warn "Truncate κ to $κ"
    end
    coldness, Dcoldness = schedule

    # init M samplers
    samplers = [Vector{Sample}() for _ = 1:M]
    nextaction = Vector{Action}(undef, M)
    queue = PriorityQueue{Int32, Float64}()
    
    for i = 1:M 
        nextaction[i] = Action(0.0, init, (first(G), 1, last(G), 0.0)) # pass correct initial score?!
        enqueue!(queue, i, 0.0)
    end

    # todo: multiply iterations by M to keep passed iteration number indep of M?
    # could also stop if one sampler has more than iterations many samples
    # but then @showprogress does not work so nicely?! 
    iterations *= M
    bestgraph = DiGraph(n)
    bestscore = 0.0 # fix if correct initial score is given above
   
    @showprogress for _ in 1:iterations 
        i = dequeue!(queue)
        nextsample = nextaction[i].apply(samplers, i, nextaction[i].τ, nextaction[i].args...)
        if nextsample.scoreval > bestscore 
            bestgraph = nextsample.g
            bestscore = nextsample.scoreval
        end
        save!(samplers[i], nextsample)
        nextaction[i] = sampleaction(samplers, i, M, balance, prior, score, σ, ρ, κ, coldness, Dcoldness)
        enqueue!(queue, i, nextaction[i].τ)
    end


    return bestgraph, samplers
end