struct Sample
    g::DiGraph
    τ::Float64
    dir::Int8
    total::Int32 
    scoreval::Float64
    alive::Bool
end
Sample(g, nextτ, dir, total, scoreval) = Sample(g, nextτ, dir, total, scoreval, true) 

struct Action
    i::Int
    τ::Float64 
    apply!::Function
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
    prevsample = samplers[i]
    g = next_CPDAG(prevsample.g, :up, x, y, T)
    return samplers[i] = Sample(g, nextτ, prevsample.dir, prevsample.total+1, prevsample.scoreval + Δscoreval)
end

function applydown(samplers, i, nextτ, x, y, H, Δscoreval)
    prevsample = samplers[i]
    g = next_CPDAG(prevsample.g, :down, x, y, H)
    return samplers[i] = Sample(g, nextτ, prevsample.dir, prevsample.total-1, prevsample.scoreval + Δscoreval) 
end

function applyflip(samplers, i, nextτ)
    prevsample = samplers[i]
    return samplers[i] = Sample(prevsample.g, nextτ, -1*prevsample.dir, prevsample.total, prevsample.scoreval)
end

function applycopy(samplers, i, nextτ, j)
    copysample = samplers[j]
    return samplers[i] = Sample(copysample.g, nextτ, copysample.dir, copysample.total, copysample.scoreval)
end

function applykill(samplers, i, nextτ)
    prevsample = samplers[i]
    return samplers[i] = Sample(prevsample.g, nextτ, prevsample.dir, prevsample.total, prevsample.scoreval, false)
end

function applynothing(samplers, i, nextτ)
    @assert false
    sample = samplers[i]
    return samplers[i] = Sample(sample.g, nextτ, sample.dir, sample.total, sample.scoreval, sample.alive)
end

# for starters without turn move
const baseline_ = [0.0]

function sampleaction(samplers, i, M, balance, prior, score, σ, ρ, κ, coldness, Dcoldness, threshold, keep) 
    # preprocess 
    prevsample = samplers[i]
    prevsample.alive || return Action(i, Inf, applynothing, ())
 
    sup, sdown, Δscorevalup, Δscorevaldown, argsup, argsdown = count_moves_new(prevsample.g, κ, balance, prior, score, coldness(prevsample.τ), prevsample.total)
    global baseline_
    # propose moves
    λdir = prevsample.dir == 1 ? sup : sdown 
    λupdown = sup + sdown 
    λflip = max(prevsample.dir*(-sup + sdown), 0.0)
    if baseline_[] - prevsample.scoreval <= 0 # assert exp(score) < 1.0, -score > 0
        baseline_[] = prevsample.scoreval 
    end
    λterm = exp(ULogarithmic, 0.0)*Dcoldness(prevsample.τ) * clamp(baseline_[] - prevsample.scoreval, eps(), threshold) # TODO: prior
    Δτdir = randexp()/(ρ*λdir)
    Δτupdown = randexp()/(σ*λupdown)
    Δτflip = randexp()/(ρ*λflip)
    Δτterm = randexp()/abs(λterm) 
    Δτmin, a = findmin((Δτdir, Δτupdown, Δτflip, Δτterm))
    A = (:dir, :updown, :flip, :term)[a]
    @assert Δτmin >= 0 
    if :dir == A
        if prevsample.dir == 1
            return Action(i, prevsample.τ + Δτdir, applyup, (argsup..., Δscorevalup))
        else 
            return Action(i, prevsample.τ + Δτdir, applydown, (argsdown..., Δscorevaldown))
        end
    end

    if :updown == A
        λup = sup
        if rand() < λup/λupdown
            return Action(i, prevsample.τ + Δτupdown, applyup, (argsup..., Δscorevalup))
        else 
            return Action(i, prevsample.τ + Δτupdown, applydown, (argsdown..., Δscorevaldown))
        end
    end
    
    if :flip == A
        return Action(i, prevsample.τ + Δτflip, applyflip, ())
    end

    if :term == A
        if rand() < keep
            return Action(i, prevsample.τ + Δτterm, applycopy, (rand(1:M),)) 
        else
            return Action(i, prevsample.τ + Δτterm, applykill, ()) 
        end    
    end

    @assert false
end
# remark: chose κ = n-1 as default
function multisampler(n, G = (DiGraph(n), 0); M = 10, balance = metropolis_balance, prior = (_,_) -> 1.0, score=UniformScore(), σ = 0.0, ρ = 1.0, κ = n - 1, baseline = 0.0, iterations = min(3*n^2, 50000), schedule=(expcoldness, Dexpcoldness), threshold=Inf, keep=1.0) #, verbose = false, save = true)
    if κ >= n 
        κ = n - 1
        @warn "Truncate κ to $κ"
    end
    coldness, Dcoldness = schedule

    global baseline_
    baseline_[] = baseline
    # init M samplers
    samplers = [Sample(G[1], 0.0, 1, G[2], 0.0) for _ = 1:M] # pass correct initial score?!
    queue = PriorityQueue{Action, Float64}()
    
    for i = 1:M 
        action = sampleaction(samplers, i, M, balance, prior, score, σ, ρ, κ, coldness, Dcoldness, threshold, keep)
        enqueue!(queue, action, action.τ)
    end

    # todo: multiply iterations by M to keep passed iteration number indep of M?
    # could also stop if one sampler has more than iterations many samples
    # but then @showprogress does not work so nicely?! 
    iterations *= M
    bestgraph = DiGraph(n)
    bestscore = 0.0 # fix if correct initial score is given above
    count = 0
    particles = M
    t = 0.0
    pr = Progress(iterations)
    @showprogress for iter in 1:iterations 
        next!(pr; showvalues = [(:M,particles), (:temp, round(schedule[1](t), sigdigits=4))])
        action = dequeue!(queue)
        t = action.τ
        count += (action.apply! == applycopy) || (action.apply! == applykill)
        if action.apply! == applykill
            particles -= 1
        end

        nextsample = action.apply!(samplers, action.i, action.τ, action.args...)
        particles == 0 && break

        if nextsample.alive && nextsample.scoreval > bestscore 
            bestgraph = nextsample.g
            bestscore = nextsample.scoreval
        end
        action = sampleaction(samplers, action.i, M, balance, prior, score, σ, ρ, κ, coldness, Dcoldness, threshold, keep)
        enqueue!(queue, action, action.τ)
        # todo: applyflip shouldn't increase counter
    end
    killratio = count/iterations
    β = schedule[1](t)
    @show particles killratio t β

    return bestgraph, [sample for sample in samplers if sample.alive]
end
