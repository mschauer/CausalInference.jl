####################################################################
# Data Structures
####################################################################

# Collection of variables to pass along to different functions in the GES algorithm.  
struct ParseData{M, F<:AbstractFloat}
    data::M # original data
    normAugData::M # standardized by the mean and std of each column, appended with ones column at end
    numFeatures::Int # number of columns
    numObservations::Int # number of rows
    penalty::F

    function ParseData(data::M, penalty) where M 

        # Determine the data type of the data inputted
        baseType = eltype(data)

        # Get the dimensions of the input data
        numObservations, numFeatures = size(data)

        # Copy the data and standardize each column
        normAugData = copy(data)
        normAugData .-= mean(normAugData, dims=1) #subtract mean
        normAugData ./= std(normAugData, dims=1) #divide by standard deviation

        # Augment a column of ones on the end for linear regression
        normAugData = [normAugData ones(baseType, numObservations)]

        # Ensure that the penalty is the same type as the data (e.g. Float32)
        penalty = baseType(penalty)

        return new{M, baseType}(data, normAugData, numFeatures, numObservations, penalty)
    end
end  

# Simple structure to hold the current edge, a subset of neighbors, and a score change
Base.@kwdef struct Step{A,B}
    edge::Graphs.SimpleEdge{A} = Edge{A}(1,2)
    subset::Vector{A} = A[]
    Δscore::B = zero(B)
end

####################################################################
# Base Function Overloads
####################################################################

# Print method to display the Step
function show(io::IO, nextstep::Step{A,B}) where {A,B}
    print(io, "Edge: $(nextstep.edge), Subset: $(nextstep.subset), Δscore: $(nextstep.Δscore)")
end

# The @memoize macro has to check if ParseData is the same argument. 
# We only ever define one immutable ParseData. We check equality of data.
==(a::T, b::T) where T <: ParseData = a.data === b.data


####################################################################
# Main Entry point for the Algorithm
####################################################################

"""
    ges(data; verbose=false)

Compute a causal graph for the given observed data using GES.
"""
function ges(data; penalty=1.0, parallel=false, verbose=false)

    # Data type / precision
    score = zero(eltype(data))

    # Parse the inputted data
    dataParsed = ParseData(data, penalty)

  
    return ges(dataParsed.numFeatures, dataParsed; score, parallel, verbose)
end

"""
    ges(n, data; score=0.0, parallel=false, verbose=false)

Estimate a causal graph for the given observed data using GES.
"""
function ges(n, data; score=0.0, parallel=false, verbose=false)
    # Create an empty graph with one node for each feature
    g = DiGraph(n)

    verbose && println("Start forward search")
    g, score = ges_search_insert!(g, score, data, parallel, verbose)

    verbose && println("Start backward search")
    g, score = ges_search_delete!(g, score, data, verbose)
   
    return g, score
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

function ges_search_insert!(g, score, data, parallel, verbose)
    # Continually add edges to the graph until the score stops increasing
    while ne(g) < nv(g)*(nv(g)-1) # there are still missing edges
        # Get the new best step
        step = find_insert(score, data, g, parallel, verbose)
        
        # If the score did not improve...
        if step.Δscore ≤ 0
            verbose && println(vpairs(g))
            break
        end
        verbose && println(step)
        # Use the insert or delete operator update the graph
        Insert!(g, step)
        score += step.Δscore
        
        # Convert the PDAG to a complete PDAG
        # Undirect all edges unless they participate in a v-structure
        vskel!(g)
        g2 = copy(g)
        # Apply the 3 Meek rules to orient some edges in the graph
        meek_rules!(g)
    end
    return g, score
end

function ges_search_delete!(g, score, data, verbose)
    # Continually remove edges to the graph until the score stops increasing
    while ne(g) > 0 # there are still edges
        
        # Get the new best step
        step = find_delete(score, data, g, verbose)
        
        # If the score did not improve...
        if step.Δscore ≤ 0
            verbose && println(vpairs(g))
            break
        end
        verbose && println(step)
        # Use the insert or delete operator update the graph
        Delete!(g, step)
        score += step.Δscore
        verbose && println(vpairs(g))

        # Convert the PDAG to a complete PDAG
        # Undirect all edges unless they participate in a v-structure
        vskel!(g)
        g2 = copy(g)
        # Apply the 3 Meek rules to orient some edges in the graph
        meek_rules!(g)
    end
    return g, score
end

best(a::Step, b::Step) = a.Δscore > b.Δscore ? a : b

function find_insert(score, data, g, parallel, verbose)
    # Loop through all possible node combinations, skip over diagonal and adjacent edges
    if parallel
        return ThreadsX.reduce(best, (findBestInsert(score, data, g, x, y) for x in vertices(g), 
                                                                               y in vertices(g) 
                                      if x != y && !isadjacent(g, x, y)),
                               init=Step(Edge(0,0), Int[], typemin(typeof(score))))
    else
        return reduce(best, (findBestInsert(score, data, g, x, y) for x in vertices(g), 
                                                                               y in vertices(g) 
                                      if x != y && !isadjacent(g, x, y)))
    end
end

function find_delete(score, data, g, verbose)

    # Loop through all possible edges
   return reduce(best, (findBestDelete(score, data, g, src(e), dst(e)) for e in edges(g) 
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
    intersect(inneighbors(g,y), outneighbors(g,y), all_neighbors(g,x))
end 


function findBestInsert(score, dataParsed, g, x, y)
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


function findBestDelete(score, dataParsed, g, x, y)
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

# score equivalent (?) oracle score
# two dag with the same cpdag need to have the same score sum
function Δscoreinsert(cpdag::DiGraph, vparents, x, v, T)
    !dsep(cpdag, x, v, vparents)
end
function Δscoredelete(cpdag::DiGraph, vparents, x, v, H)
    dsep(cpdag, x, v, vparents)
end

Δscoreinsert(data, parents, x, v, _) = Δscore(data, parents, x, v)
Δscoredelete(data, parents, x, v, _) = -Δscore(data, parents, x, v)

Δscore(data, parents, x, v) = local_score(data, sort(push!(copy(parents), x)), v) - local_score(data, sort(parents), v)



export score_dag
function score_dag(g, data) # g dag
    s = 0.0
    for v in vertices(g)
        s += score(data, inneighbors(g, v), v)
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
    score(os::GaussianScore, p, v)

Local Gaussian BIC score.
"""
@memoize LRU(maxsize=1_000_000) function local_score(os::GaussianScore, p, v)
    C = os.C 
    penalty = os.penalty
    n = os.n
    k = length(p) # dimension
    if k == 0 
        Cp = C[v, v]
    else # compute conditional correlation
        Cp = C[v, v] - dot(C[v, p], C[p, p]\C[p, v])
    end
    - n*log(Cp) - penalty*k*log(n) 
end



@memoize LRU(maxsize=1_000_000) function local_score(dataParsed::ParseData{Matrix{A}}, nodeParents, node) where A

    # Unpack some variables from the dataParsed structure
    n = A(dataParsed.numObservations) #convert datatype
    data = dataParsed.normAugData
    p = dataParsed.penalty

    # The last column of dataParsed.normAugData is all ones which is required for a linear regression with an intercept. If there are no node parents, model the child node with just the intercept, else use the parents and the intercept
    if isempty(nodeParents)
        parentsAndIncept = [dataParsed.numFeatures+1]
    else
        parentsAndIncept = [nodeParents; dataParsed.numFeatures+1]
    end

    # To calculate the score we need a mean-squared 
    # error which we can get by regessing the the child node onto the parents
    # Create variables for a linear regression y = X*b

    # Use views to avoid creating copies of the data in memory
    # X is the design matrix, augmented with a column of ones at the end
    # X is also been standardized so mean(columns)=0 and std(column)=1
    # y is data from the child node being tested
    @views begin
        y = data[:,node]
        X = data[:,parentsAndIncept]
    end

    # Perform a linear regression
    b = X \ y

    # Get the estimation
    ŷ = X*b

    # Next we want to calculate the log likelihood that these are the parents of our node
    # score = log(P(data|Model)) ≈ -BIC/2
    # because we're only comparing log likelihoods we'll ignore the 1/2 factor
    # when P(⋅) is Gaussian, log(P(data|Model)) takes the form:
    # -k⋅log(n) - n⋅log(mse)
    # k is the number of free parameters and mse is mean squared error
    k = length(parentsAndIncept) #includes the intercept
    mse = sum(x->x^2, y-ŷ) / n

    # return the final score we want to maximize (which is technically -2BIC)
    return -n*log(mse) - p*k*log(n) 
end



