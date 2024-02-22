@inline subset_range(n) = 0:2^n-1
@inline subset_size(n) = 2^n
@inline element_range(n) = 0:n-1

get_subset(j, n) = [k+1 for k = element_range(n) if ((1 << k) & j) != 0]
# for vertex i and subset j find the indexed parent set
get_parents(i, j, n) = [k+1 for k = element_range(n) if (k < i && 1 << k & j != 0) || (k > i && 1 << (k-1) & j != 0)]

function compute_bestparents(i, n, localscores) 
    bestscores = zeros(subset_size(n-1))
    bestparents = zeros(Int64, subset_size(n-1))
    for j = subset_range(n-1)
        bestscores[j+1] = localscores[i+1][j+1]
        bestparents[j+1] = j
        for k = element_range(n-1)
            subj = j & ~(1 << k)
            subj == j && continue
            if bestscores[subj+1] > bestscores[j+1]
                bestscores[j+1] = bestscores[subj+1]
                bestparents[j+1] = bestparents[subj+1]
            end
        end
    end
    return bestparents
end

function compute_bestsinks(n, bestparents, localscores)
    bestscores = zeros(subset_size(n))
    bestsinks = -1 * ones(Int64, subset_size(n))
    for j = subset_range(n)
        for k = element_range(n)
            1 << k & j == 0 && continue
            subj = j & ~(1 << k)
            before = (subj & (1 << k - 1))
            after = subj - before
            psubj = before + (after >> 1)
            subscore = bestscores[subj+1] + localscores[k+1][bestparents[k+1][psubj+1]+1]
            if bestsinks[j+1] == -1 || subscore > bestscores[j+1]
                bestscores[j+1] = subscore 
                bestsinks[j+1] = k
            end
        end
    end
    return bestsinks
end 

function compute_ordering(n, bestsinks)
    ordering = zeros(Int64, n)
    remaining = subset_size(n) - 1
    for i = reverse(element_range(n))
        ordering[i+1] = bestsinks[remaining+1]
        remaining = remaining - 1 << ordering[i+1]
    end
    return ordering
end

function compute_network(n, ordering, bestparents)
    g = SimpleDiGraph(n)
    predecs = 0
    for i in element_range(n)
        before = (predecs & (1 << ordering[i+1] - 1))
        after = predecs - before
        idx = before + (after >> 1)
        parents = get_parents(ordering[i+1], bestparents[ordering[i+1]+1][idx+1], n)
        for p in parents
            add_edge!(g, p, ordering[i+1]+1)
        end
        predecs += 1 << ordering[i+1]
    end
    return g
end

#########################################
# Main Entry point
# #######################################

""" 
    exactscorebased(X; method=:gaussian_bic, penalty=0.5, parallel=false, verbose=false)

Compute a CPDAG for the given observed data `X` (variables in columns) using the exact algorithm proposed by Silander and Myllymäki (2006) for optimizing the BIC score (or any decomposable score).  
The complexity is n*2^n and the algorithm should scale up to ~20-25 variables, afterwards memory becomes a problem.

* Silander, T., & Myllymäki, P. (2006). A simple approach for finding the globally optimal Bayesian network structure. In Uncertainty in Artificial Intelligence (pp. 445-452). 
"""
function exactscorebased(X::AbstractMatrix; method=:gaussian_bic, penalty=0.5, parallel=false, verbose=false)
    (N, n) = size(X)
    n > 25 && @warn "algorithm will take a long time to terminate and needs a lot of memory"
    n > 64 && @error "algorithm only works up to 64 variables (and will likely not be feasible for more than 25-30 variables)"
    
    if method == :gaussian_bic
        C = Symmetric(cov(X, dims = 1, corrected = false))
        S = GaussianScore(C, N, penalty)
        return exactscorebased(n, (p, v) -> local_score(S, p, v) ; parallel, verbose)
    elseif method == :gaussian_bic_raw
        S = GaussianScoreQR(X, penalty)
        return exactscorebased(n, (p, v) -> local_score(S, p, v); parallel, verbose)
    else 
        throw(ArgumentError("method=$method"))
    end
end
exactscorebased(X; method=:gaussian_bic, penalty=0.5, parallel=false, verbose=false) = exactscorebased(Tables.matrix(X); method, penalty, parallel, verbose)

function exactscorebased(n, local_score; parallel=false, verbose=false)
    # Step 1: calculate score for all n*2^(n-1) different (variable, variable set) pairs (variable and parents)
    localscores = [Vector{Float64}() for _ = 1:n]
    if parallel 
        Threads.@threads for i = element_range(n)
            localscores[i+1] = [local_score(get_parents(i, j, n), i+1) for j = subset_range(n-1)]
        end
    else 
        for i = element_range(n)
            localscores[i+1] = [local_score(get_parents(i, j, n), i+1) for j = subset_range(n-1)]
        end
    end

    # Step 2: find the best parents for all vertices and subsets (of candidate parents)
    bestparents = [Vector{Int64}() for _ = 1:n]
    if parallel 
        Threads.@threads for i = element_range(n)
            bestparents[i+1] = compute_bestparents(i, n, localscores)
        end
    else 
        for i = element_range(n)
            bestparents[i+1] = compute_bestparents(i, n, localscores)
        end
    end

    # Step 3: find the best sink for all subsets of variables
    bestsinks = compute_bestsinks(n, bestparents, localscores)

    # Step 4: find best ordering of vertices
    ordering = compute_ordering(n, bestsinks)
    verbose && println("ordering" * string(ordering))

    # Step 5: find best network
    return alt_cpdag(compute_network(n, ordering, bestparents))
end

