using DataStructures

struct Sample
    g::DiGraph
    τ::Float64
    dir::Int8
    total::Int32 
    scoreval::Float64
end

struct Action
    τ::Float64 
    func::Function
    args::Tuple{Any}
end

function init(_, _, nextτ, g, dir, total, scoreval)
    return Sample(g, nextτ, dir, total, 0.0, scoreval) 
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

function applycopy(samplers, i, nextτ, j)
    @assert i != j 
    return Sample(samplers[j].g, nextτ, samplers[j].dir, samplers[j].total, samplers[j].scoreval)
end

# for starters without turn move
function sampleevent(samplers, i, M, balance, prior, score, coldness, σ, ρ, κ, ) 
    # preprocess 
    prevsample = last(samplers[i])
    sup, sdown, Δscorevalup, Δscorevaldown, argsup, argsdown = count_moves_new(g, κ, balance, prior, score, coldness, prevsample.total)

    # propose moves
    λdir = prevsample.dir == 1 ? sup : sdown 
    λupdown = sup + sdown 
    λflip = max(prevsample.dir*(-sup + sdown), 0.0)
    λterm = 0.0 # TODO

    Δτdir = randexp()/(ρ*λdir)
    Δτupdown = randexp()/(σ*λupdown)
    Δτflip = randexp()/(ρ*λflip)
    Δτterm = randexp()/(dβdτ * β * scorevalue) # TODO # prior?

    Δτmin = min(Δτdir, Δτupdown, Δτflip, Δτterm)

    if Δτdir == Δτmin 
        if prevsample.dir == 1
            return prevsample.τ + Δτdir, applyup, (argsup..., Δscorevalup)
        else 
            return prevsample.τ + Δτdir, applydown, (argsdown..., Δscorevaldown)
        end
    end

    if Δτupdown == Δτmin
        λup = sup
        if rand() < λup/λupdown
            return prevsample.τ + Δτupdown, applyup, (argsup..., Δscorevalup)
        else 
            return prevsample.τ + Δτupdown, applydown, (argsdown..., Δscorevaldown)
        end
    end
    
    if Δτflip == Δτmin 
        return prevsample.τ + Δτflip, applyflip, ()
    end

    if Δτterm == Δτmin 
        j = i
        # not pretty
        while j == i 
            j = rand(1:M)
        end
        return prevsample.τ + Δτterm, applycopy, (j) # maybe sample copied process here directly
    end

    @assert false
end

# todo: κ = n-1 as default
function multisampler(n, G = (DiGraph(n), 0); M = 10, balance = metropolis_balance, prior = (_,_) -> 1.0, score=UniformScore(), coldness = 1.0, σ = 0.0, ρ = 1.0, κ = min(n - 1, 10), iterations = 10, verbose = false, save = true)
    if κ >= n 
        κ = n - 1
        @warn "Truncate κ to $κ"
    end
    
    # init M samplers
    # write init function
    samplers = [Vector{Sample}() for _ = 1:M]
    nextaction = Vector{Sample}(undef, M)
    queue = PriorityQueue{Int32, Float64}()
    
    for i = 1:M 
        nextevent[i] = Action(0.0, init, (G, 1, 0, 0.0)) # pass correct initial score?!
        enqueue!(queue, i, 0.0)
    end

    # todo: multiply iterations by M to keep passed iteration number indep of M?
    # could also stop if one sampler has more than M samples
    # but then @showprogress does not work so nicely?! 
    iterations *= M

    @showprogress for _ in 1:iterations 
        i = dequeue!(queue)
        nextsample = nextaction[i].func(samplers, i, nextaction[i].τ, nextaction[i].args)
        push!(samplers[i], nextsample)
        nextaction[i] = sampleaction(samplers, i, M, balance, prior, score, coldness, σ, ρ, κ)
        enqueue!(queue, i, nextaction[i].τ)
    end

    # postprocess
end
