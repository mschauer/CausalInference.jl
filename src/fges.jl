####################################################################
# Data Structures
####################################################################

# Collection of variables to pass along to different functions in the FGES algorithm.  
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
Base.@kwdef mutable struct Step{A,B}
    edge::Edge = Edge(1,2)
    subset::Vector{A} = A[]
    Δscore::B = zero(B)
end


####################################################################
# Base Function Overloads
####################################################################

# Print method to display the Step
function show(io::IO, newStep::Step{A,B}) where {A,B}
    print(io, "Edge: $(newStep.edge), Subset: $(newStep.subset), Δscore: $(newStep.Δscore)")
end

# The @memoize macro has to check if ParseData is the same argument. 
# We only ever define one immutable ParseData. We check equality of data.
==(a::T, b::T) where T <: ParseData = a.data == b.data


####################################################################
# Main Entry point for the Algorithm
####################################################################

"""
    fges(data; penalty = 1.0, verbose=false)

Compute a causal graph for the given observed data using GES.
"""
function fges(data; penalty = 1.0, verbose=false)

    # data type / precision
    scoreT = eltype(data)

    # Parse the inputted data
    dataParsed = ParseData(data, penalty)

    # Create an empty graph with one node for each feature
    g = DiGraph(dataParsed.numFeatures)

    return fges_internal(g, scoreT, dataParsed; penalty, verbose)
end

function fges_internal(g, scoreT, data; penalty = 1.0, verbose=false)

    verbose && println("Start forward search")
    # Perform the forward search 
    
    ges_search!(g, Step{Int, scoreT}(), data, Insert!, verbose)

    verbose && println("Start backward search")
    # Perform the backward search 
    ges_search!(g, Step{Int, scoreT}(), data, Delete!, verbose)
    
    return g # Return the graph

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

    Deletes!(g, x, y, T)

Deletes x-y or x->y and directs previously undirected edges t->y, t ∈ T.
Here x and y are not adjacent and T are undirected-neighbors of y 
that are not adjacent to x.
"""
function Delete!(g, x, y, H)

    # Remove directed or undirected edges (x→y and x-y)
    rem_edge!(g, x → y)
    rem_edge!(g, y → x)
    
    # Orient all vertices in H toward x and y
    for h ∈ H
        orientedge!(g, x → h) 
        orientedge!(g, y → h)
    end

    return nothing
end

function Insert!(g, newStep::Step) 
    x, y = Pair(newStep.edge)
    T = newStep.subset
    Insert!(g, x, y, T)
end
function Delete!(g, newStep::Step) 
    x, y = Pair(newStep.edge)
    T = newStep.subset
    Delete!(g, x, y, T)
end

####################################################################
# Forward and Backward search
####################################################################

function ges_search!(g, newStep, data, operator, verbose)

    # Continually add/remove edges to the graph until the score stops increasing
    while true
        
        # Get the new best step (depends on if we're inserting or deleting)
        if operator == Insert!
            find_insert!(newStep, data, g, verbose)
        else
            find_delete!(newStep, data, g, verbose)
        end
        
        
        # If the score did not improve...
        if newStep.Δscore ≤ 0
            #...stop searching
            break
        end
        
        # Use the insert or delete operator update the graph
        operator(g, newStep)
        
        # Convert the PDAG to a complete PDAG
        # Undirect all edges unless they participate in a v-structure
        vskel!(g)
        # Apply the 4 Meek rules to orient some edges in the graph
        meek_rules!(g)

        # Reset the score
        newStep.Δscore = 0
        
    end

    return nothing
end


function find_insert!(newStep, data, g, verbose)
    # Loop through all possible node combinations
    for (x,y) in allpairs(vertices(g)) 
        # Skip over adjacent edges if we're trying to insert
        if !isadjacent(g,x,y)
            #For this pair of nodes, the following function will:
                #(1) check if a valid operator exists
                #(2) score all valid operators and return the one with the highest score
            (bestSubset, bestScore) = findBestInsert(newStep, data, g, x, y)

            # If the best valid operator was better than any previously found...
            if bestScore > newStep.Δscore 
                #...update newStep
                newStep.edge = Edge(x,y)
                newStep.subset = bestSubset
                newStep.Δscore = bestScore

                verbose && println(newStep)
            end
        end
    end
end

function find_delete!(newStep, data, g, verbose)
    # Loop through all possible node combinations
    for e in edges(g) # go through edges
        x, y = Pair(e)
        if y < x && has_edge(g, y, x) # undirected only once
            continue
        end
         # Check if there is an edge x->y between two vertices
        if has_edge(g, x, y)
            #For this pair of nodes, the following function will:
                #(1) check if a valid operator exists
                #(2) score all valid operators and return the one with the highest score
                #findBestOperation is either "findBestInsert" or "findBestDelete"
            (bestSubset, bestΔ) = findBestDelete(newStep, data, g, x, y)

            # If the best valid operator was better than any previously found...
            if bestΔ > newStep.Δscore 
                newStep.edge = Edge(x,y)
                newStep.subset = bestSubset
                newStep.Δscore = bestΔ

                verbose && println(newStep)
            end
        end
    end
end



# for x-y, get undirected neighbors of y connected to x
calcNAyx(g, y::Integer, x::Integer) = intersect(inneighbors(g,y), outneighbors(g,y), all_neighbors(g,x))
calcNAyx(g, edge) = calcNAyx(g, dst(edge), src(edge))

#for x-y, undirected neighbors of y not connected to x
calcT(g, y::Integer, x::Integer) = setdiff(neighbors_undirected(g,y), all_neighbors(g,x), x)
calcT(g, edge) = calcT(g, dst(edge), src(edge))

function tails_and_adj_neighbors(g, x, y)
    Nb = neighbors_undirected(g, y)
    a = Bool[isadjacent(g, t, x) for t in Nb]
    sort!(Nb[.~ a]), sort!(Nb[a])
end 
function adj_neighbors(g, x, y)
    Nb = neighbors_undirected(g, y)
    a = Bool[isadjacent(g, t, x) for t in Nb]
    sort!(Nb[a])
end 


function findBestInsert(step, dataParsed, g, x, y)
    isblocked(g, x, y, nodesRemoved) = !has_a_path(g, [x], y, nodesRemoved)

    Tyx = calcT(g, y, x)
    NAyx = calcNAyx(g, y, x)

    Ta, Nb = tails_and_adj_neighbors(g, x, y)
    @assert sort(Tyx) == Ta
    @assert sort(NAyx) == Nb

    # Create two containers to hold the best found score and best subset of Tyx
    bestΔ = zero(step.Δscore)
    bestT = Vector{Int}()

    # Keep a list of invalid sets
    invalid = Vector{Vector{Int}}()
    
    # Loop through all possible subsets of Tyx
    for T in powerset(Tyx)
        if checkSupersets(T, invalid)
            NAyxT = NAyx ∪ T
            # Validity of insert operator
            if isclique(g, NAyxT) && isblocked(g, y, x, NAyxT)

                # Score the Insert
                PAy = inneighbors(g, y)
                PAy⁺ = PAy ∪ x
                newΔ = score(dataParsed, NAyxT ∪ PAy⁺, y) - score(dataParsed, NAyxT ∪ PAy, y)
                
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

    return (bestT, bestΔ)
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


function findBestDelete(step, dataParsed, g, x, y)
    # Calculate two (possibly empty) sets of nodes
    # NAxy: any nodes that are undirected neighbors of y and connected to x by any edge
    # Hyx: any subset of the undirected neighbors of y that are connected to x
    NAyx = calcNAyx(g, y, x)
    Na = adj_neighbors(g, x, y)
    Hyx = NAyx
    @assert sort(NAyx) == Na

    # Best found score difference
    # and best subset of Hyx
    bestΔ = zero(step.Δscore)
    bestH = Vector{Int}()

    # Loop through all possible subsets of Hyx
    for H in powerset(Hyx)
        # Calculate NAyx \ {H}
        NAyx_H = setdiff(NAyx, H)

        # Check if the operator is valid
        if isclique(g, NAyx_H)

            # Score the valid operator 
            PAy = inneighbors(g, y)
            PAy⁻ = setdiff(PAy, x)
            newΔ = score(dataParsed, NAyx_H ∪ PAy⁻, y) - score(dataParsed, NAyx_H ∪ PAy, y)
            
            if newΔ > bestΔ
                bestH = H
                bestΔ = newΔ
            end
        end
    end

    return (bestH, bestΔ)
end


####################################################################
# Scoring function
####################################################################

# score equivalent (?) oracle score
# two dag with the same cpdag need to have the same score sum
function score(cpdag::DiGraph, vparents, v)
    uparents = neighbors_undirected(cpdag, v) 
    # possible parents are good (otherwise we could learn only the v structures)
    # impossible parents are bad
    iparents = setdiff(vparents, inneighbors(cpdag, v))
    # new v-structures are bad
    ns = vparents ∩ uparents
    protected = falses(length(ns)) # mark in-neighbours which are v structures
    for (i, j) in combinations(1:length(ns), 2) 
        if !isadjacent(cpdag, ns[i], ns[j]) 
            protected[i] = protected[j] = true
        end
    end

    return length(ns) - length(protected) -  length(iparents)
end
# missing necessary parents are bad? add - length(setdiff(parents_(cpdag, v), vparents))

@memoize LRU(maxsize=1_000_000) function score(dataParsed::ParseData{Matrix{A}}, nodeParents, node) where A

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

    # To calculate the score we need a mean-squared error which we can get by regessing the the child node onto the parents
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
    return -p*k*log(n) - n*log(mse)
end



