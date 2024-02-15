

# Simple structure to hold the current edge, a subset of neighbors, and a score change
Base.@kwdef struct Step{A,B}
    edge::Graphs.SimpleEdge{A} = Edge{A}(1,2)
    subset::Vector{A} = A[]
    Δscore::B = zero(B)
end

# Print method to display the Step
function show(io::IO, nextstep::Step{A,B}) where {A,B}
    print(io, "$(nextstep.edge), Subset $(nextstep.subset), Δscore $(round(nextstep.Δscore, sigdigits=5))")
end



####################################################################
# Main Entry point for the Algorithm
####################################################################

"""
    ges(X; method=:gaussian_bic, penalty=0.5, parallel=false, verbose=false)

Compute a causal graph for the given observed data `X` (variables in columns) using GES.
Returns the CPDAG, the score improvement relative to the empty graph and time measurements
of first and second phase.
"""
function ges(X::AbstractMatrix; method=:gaussian_bic, penalty=0.5, parallel=false, verbose=false)
    score = zero(eltype(X)) # initial score
    n, d = size(X)
    d ≥ n && @warn "High dimensional data (n ≤ p), ges might not terminate."
    if method == :gaussian_bic
        C = Symmetric(cov(X, dims = 1, corrected=false))
        return ges(d, GaussianScore(C, n, penalty); score, parallel, verbose)
    elseif method == :gaussian_bic_raw
        return ges(d, GaussianScoreQR(X, penalty); score, parallel, verbose)
    else
        throw(ArgumentError("method=$method"))
    end
end
ges(X; method=:gaussian_bic, penalty=0.5, parallel=false, verbose=false) = ges(Tables.matrix(X); method, penalty, parallel, verbose)


"""
    ges(n, local_score; score=0.0, parallel=false, verbose=false)

Internal method called by `ges`. 
"""
function ges(n, data; score=0.0, parallel=false, verbose=false)
    # Create an empty graph with one node for each feature
    g = DiGraph(n)
    parallel && Threads.nthreads() == 1 && @warn "Only one thread available"
    verbose && println("Start forward search")
    t1 = @elapsed g, score = ges_forward_search!(g, score, data, parallel, verbose)

    verbose && println("Start backward search")
    t2 = @elapsed g, score = ges_backward_search!(g, score, data, verbose)
   
    return g, score, (t1, t2)
end

####################################################################
# Insert and Delete Operators
####################################################################
"""
    Insert!(g, x, y, T)

Inserts x->y and directs previously undirected edges t->y, t ∈ T.
Here x and y are not adjacent and T are undirected-neighbors of y 
that are not adjacent to x.
"""
function Insert!(g, x, y, T)
    add_edge!(g, x → y)
    
    # Orient all edges in T incident into child node
    for t ∈ T
        orientedge!(g, t → y)
    end
    return g
end

"""
    Delete!(g, x, y, H)

Deletes x-y or x->y and directs previously undirected edges x->h and y->h
for h in H.
"""
function Delete!(g, x, y, H)

    # Remove directed or undirected edges (x→y and x-y)
    rem_edge!(g, x → y)
    rem_edge!(g, y → x)
    
    # Orient all vertices in H toward x and y
    for h ∈ H
        if has_both(g, x, h) # reading literally Definition 13, Chickering
            orientedge!(g, x → h) 
        end
        orientedge!(g, y → h)
    end

    return nothing
end

function Insert!(g, nextstep::Step) 
    x, y = Pair(nextstep.edge)
    T = nextstep.subset
    Insert!(g, x, y, T)
end
function Delete!(g, nextstep::Step) 
    x, y = Pair(nextstep.edge)
    H = nextstep.subset
    Delete!(g, x, y, H)
end

####################################################################
# Forward and Backward search
####################################################################

function ges_forward_search!(g, score, data, parallel, verbose)
    # Continually add edges to the graph until the score stops increasing
    while ne(g) < nv(g)*(nv(g)-1) # there are still missing edges
        # Get the new best step
        step = find_best_insert(score, data, g, parallel, verbose)
        
        # If the score did not improve...
        if step.Δscore ≤ 0
            #verbose && println(vpairs(g))
            break
        end
        score += step.Δscore
        verbose && println(step, " ", round(score, sigdigits=5))
        # Use the insert or delete operator update the graph
        Insert!(g, step)
        
        # Convert the PDAG to a complete PDAG
        # Undirect all edges unless they participate in a v-structure
        vskel!(g)
        # Apply the 3 Meek rules to orient some edges in the graph
        meek_rules!(g)
    end
    return g, score
end

function ges_backward_search!(g, score, data, verbose)
    # Continually remove edges to the graph until the score stops increasing
    while ne(g) > 0 # there are still edges
        
        # Get the new best step
        step = find_best_delete(score, data, g, verbose)
        
        # If the score did not improve...
        if step.Δscore ≤ 0
            #verbose && println(vpairs(g))
            break
        end
        score += step.Δscore
        verbose && println(step, " ", round(score, sigdigits=5))
        # Use the insert or delete operator update the graph
        Delete!(g, step)
        #verbose && println(vpairs(g))

        # Convert the PDAG to a complete PDAG
        # Undirect all edges unless they participate in a v-structure
        vskel!(g)
        # Apply the 3 Meek rules to orient some edges in the graph
        meek_rules!(g)
    end
    return g, score
end

best(a::Step, b::Step) = a.Δscore > b.Δscore ? a : b

function find_best_insert(score, data, g, parallel, verbose)
    # Loop through all possible node combinations, skip over diagonal and adjacent edges
    if parallel
        return ThreadsX.reduce(best, (score_edge_inserts(score, data, g, x, y) for x in vertices(g), 
                                                                               y in vertices(g) 
                                      if x != y && !isadjacent(g, x, y)),
                               init=Step(Edge(0,0), Int[], typemin(typeof(score))))
    else
        return reduce(best, (score_edge_inserts(score, data, g, x, y) for x in vertices(g), 
                                                                               y in vertices(g) 
                                      if x != y && !isadjacent(g, x, y)))
    end
end

function find_best_delete(score, data, g, verbose)

    # Loop through all possible edges
   return reduce(best, (score_edge_deletions(score, data, g, src(e), dst(e)) for e in edges(g) 
                        if !(dst(e) < src(e) && has_edge(g, reverse(e)))  # undirected only once ✓
                            && has_edge(g, e)))
end



# for x-y, get undirected neighbors of y connected to x
#calcNAyx(g, y::Integer, x::Integer) = intersect(inneighbors(g,y), outneighbors(g,y), all_neighbors(g,x))
#for x-y, undirected neighbors of y not connected to x
#calcT(g, y::Integer, x::Integer) = setdiff(neighbors_undirected(g,y), all_neighbors(g,x), x)

function tails_and_adj_neighbors(g, x, y)
    Nb = neighbors_undirected(g, y)
    a = Bool[isadjacent(g, t, x) for t in Nb]
    Nb[.~ a], Nb[a]
end 
function adj_neighbors(g, x, y)
#    a = intersect(inneighbors(g,y), outneighbors(g,y), all_neighbors(g,x))
    sorted_intersect_(neighbors_undirected(g,y), all_neighbors(g,x))
end 


function score_edge_inserts(score, dataParsed, g, x, y)
    isblocked(g, x, y, nodesRemoved) = !has_a_path(g, [x], y, nodesRemoved)

    Tyx, NAyx = tails_and_adj_neighbors(g, x, y)


    # Best found score and best subset of Tyx
    bestΔ = typemin(typeof(score))
    bestT = Vector{Int}()

    # Keep a list of invalid sets
    invalid = Vector{Vector{Int}}()
    
    # Loop through all possible subsets of Tyx
    for T in powerset(Tyx)
        if checkSupersets(T, invalid)
            NAyxT = NAyx ∪ T
            # Validity of insert operator
            if isclique(g, NAyxT) && isblocked(g, y, x, NAyxT)

                # Score the insert operator
                PAy = parents(g, y)
                newΔ = Δscoreinsert(dataParsed, NAyxT ∪ PAy, x, y, T)
                
                # Save the new score if it was better than any previous
                if newΔ > bestΔ
                    bestT = T
                    bestΔ = newΔ
                end
            end
        else
            # Record that the subset T is invalid
            push!(invalid, T)
        end
    end

    return Step(Edge(x,y), bestT, bestΔ)
end

# Check if the set T is a superset of any invalid set
function checkSupersets(T, invalid)
    for i ∈ invalid
        if i ⊆ T
            return false
        end
    end
    return true
end


function score_edge_deletions(score, dataParsed, g, x, y)
    # Calculate two (possibly empty) sets of nodes
    # NAxy: any nodes that are undirected neighbors of y and connected to x by any edge
    # Hyx: any subset of the undirected neighbors of y that are connected to x
    NAyx = adj_neighbors(g, x, y)
    Hyx = NAyx

    # Best found score difference
    # and best subset of Hyx
    bestΔ = zero(score)
    bestH = Vector{Int}()

    # Loop through all possible subsets of Hyx
    for H in powerset(Hyx)
        # Calculate NAyx \ {H}
        NAyx_H = setdiff(NAyx, H)

        # Check if the operator is valid
        if isclique(g, NAyx_H)

            # Score the operator 
            PAy = parents(g, y)
            PAy⁻ = setdiff(PAy, x)
            newΔ = Δscoredelete(dataParsed, NAyx_H ∪ PAy⁻, x, y, H)
            
            if newΔ > bestΔ
                bestH = H
                bestΔ = newΔ
            end
        end
    end

    return Step(Edge(x,y), bestH, bestΔ)
end


####################################################################
# Scoring function
####################################################################


Δscoreinsert(data, parents, x, v, _) = Δscore(data, parents, x, v)
Δscoredelete(data, parents, x, v, _) = -Δscore(data, parents, x, v)

Δscore(data, parents, x, v) = local_score(data, sort(push!(copy(parents), x)), v) - local_score(data, sort(parents), v)



export score_dag
function score_dag(g, data) # g dag
    s = 0.0
    for v in vertices(g)
        s += local_score(data, inneighbors(g, v), v)
    end
    s
end

struct GaussianScore{T, S<:AbstractMatrix{T}}
    C::S # correlation compares identically
    n::Int # hypothetical number of obs
    penalty::Float64
    hash::UInt
end
GaussianScore(C, n, penalty) = GaussianScore(C, n, penalty, hash((C, n, penalty)))
export GaussianScore
import Base.:(==), Base.hash
is_equal(a::T, b::T) where T <: GaussianScore = ((a.C === b.C) || is_equal(a.C, b.C)) && is_equal(a.n, b.n) &&  is_equal(a.penalty, b.penalty)
hash(a::GaussianScore, u::UInt) = hash(a.hash, u)

# compare https://github.com/py-why/causal-learn/blob/f51195473b316662b6f7dce68cd73d734766a6a3/causallearn/score/LocalScoreFunction.py
"""
    local_score(os::GaussianScore, p, v)

Local Gaussian BIC score. Memoized for `GaussianScore{Float64, Symmetric{Float64, Matrix{Float64}}}`.
"""
function local_score(os::GaussianScore, p, v)
    length(p) > 2 && return local_score_mem(os, p, v)
    local_score_(os, p, v)
end
@memoize LRU{Tuple{Tuple{GaussianScore{Float64, Symmetric{Float64, Matrix{Float64}}}, Vector{Int64}, Int64}, Tuple{}}, Float64}(maxsize=100_000) function local_score_mem(os::GaussianScore{Float64, Symmetric{Float64, Matrix{Float64}}}, p, v)
    local_score_(os, p, v)
end
function local_score_mem(os::GaussianScore, p, v)
    local_score_(os, p, v)
end

function local_score_(os::GaussianScore, p, v)
    k = length(p)
    C = os.C 
    penalty = os.penalty
    n = os.n
    if k == 0 
        Cp = C[v, v]
    elseif k == 1
        p_ = p[]
        c = C[p_, v]
        Cp = C[v, v] - c*(C[p_, p_]\c)
    else # compute conditional correlation
        c = @view C[p, v]
        Cp = C[v, v] - dot(c, (@view C[p, p])\c)
    end
    (-n*(1 + log(max(0,Cp))) - penalty*(1  + k)*log(n))/2
end



struct GaussianScoreQR{T, S<:AbstractMatrix{T}}
    X::S # Data matrix
    penalty::Float64
    hash::UInt
end
GaussianScoreQR(X, penalty) = (Xc = X .-  mean(X, dims=1); GaussianScoreQR(Xc, penalty, hash((Xc, penalty))))
export GaussianScoreQR
is_equal(a::T, b::T) where T <: GaussianScoreQR = ((a.X === b.X) || is_equal(a.X, b.X)) &&  is_equal(a.penalty, b.penalty)
hash(a::GaussianScoreQR, u::UInt) = hash(a.hash, u)


@memoize LRU(maxsize=1_000_000) function local_score(os::GaussianScoreQR, p, v)
    X = os.X 
    penalty = os.penalty
    n = size(X, 1)
    k = length(p) # dimension
    y = @view X[:, v]
    if k == 0 
        Cp = var(y; mean=0.0, corrected=false)
    else # compute conditional correlation
        x = @view X[:, p]
        Cp = var(y - x*(x\y);  mean=0.0, corrected=false) 
    end
    (-n*(1 + log(max(0,Cp))) - penalty*(1 + k)*log(n))/2
end





