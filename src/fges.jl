####################################################################
# Data Structures
####################################################################

#Collection of variables to pass along to different functions in the FGES algorithm.  
struct ParseData{M, F<:AbstractFloat}
    data::M #orginal data
    normAugData::M #standardized by the mean and std of each column, appended with ones column at end
    numFeatures::Int #number of columns
    numObservations::Int #number of rows
    penalty::F

    function ParseData(data::M, penalty) where M 

        #Determine the data type of the data inputted
        baseType = eltype(data)

        #Get the dimensions of the input data
        numObservations, numFeatures = size(data)

        #Copy the data and standardize each column
        normAugData = copy(data)
        normAugData .-= mean(normAugData, dims=1) #subtract mean
        normAugData ./= std(normAugData, dims=1) #divide by standard deviation

        #Augment a column of ones on the end for linear regression
        normAugData = [normAugData ones(baseType, numObservations)]

        #Ensure that the penalty is the same type as the data (e.g. Float32)
        penalty = baseType(penalty)

        return new{M, baseType}(data, normAugData, numFeatures, numObservations, penalty)
    end
end  

#Simple structure to hold the current edge, a subset of neighbors, and a score change
Base.@kwdef mutable struct Step{A,B}
    edge::Edge = Edge(1,2)
    subset::Vector{A} = A[]
    Δscore::B = zero(B)
end


####################################################################
# Base Function Overloads
####################################################################

#Print method to display the Step
function show(io::IO, newStep::Step{A,B}) where {A,B}
    print(io, "Edge: $(newStep.edge), Subset: $(newStep.subset), Δscore: $(newStep.Δscore)")
end

#This could be dangerous. The @memoize macro has to check if ParseData is the same argument. We only ever define one immutable ParseData, so it will never change. Any equality check should be true.
==(a::T, b::T) where T <: ParseData = true


####################################################################
# Main Entry point for the Algorithm
####################################################################

"""
    fges(data; penalty = 1.0, verbose=false)
Compute a causal graph for the given observed data.
"""
function fges(data; penalty = 1.0, verbose=false)

    #Parse the inputted data
    dataParsed = ParseData(data, penalty)

    #Create an empty graph with one node for each feature
    g = DiGraph(dataParsed.numFeatures)

    verbose && println("Start forward search")
    #Perform the forward search 
    Search!(g, dataParsed, Insert!, verbose)

    verbose && println("Start backward search")
    #Perform the backward search 
    Search!(g, dataParsed, Delete!, verbose)
    
    #Return the graph
    return g
end

####################################################################
# Insert and Delete Operators
####################################################################

function Insert!(g, newStep::Step) 
    edge = newStep.edge
    T = newStep.subset

    #Add a directed edge x→y
    add_edge!(g,edge)
    
    #Orient all edges incident into child node
    y = dst(edge)
    for t ∈ T
        orientedge!(g,t,y) #t→y
    end

    return nothing
end


function Delete!(g, newStep::Step)
    edge = newStep.edge
    H = newStep.subset

    #remove directed and unidirected edges (x→y and x-y)
    rem_edge!(g, edge)
    rem_edge!(g, reverse(edge))
    
    #Orient all vertices in H toward x and y
    (x, y) = src(edge), dst(edge)
    for h ∈ H
        orientedge!(g,x,h) #x→h
        orientedge!(g,y,h) #y→h
    end

    return nothing
end

####################################################################
# Forward and Backward search
####################################################################

function Search!(g, dataParsed::ParseData{Matrix{A}}, operator, verbose) where A

    #Create a container to hold information about the next step
    newStep = Step{Int,A}()

    #Continually add/remove edges to the graph until the score stops increasing
    while true
        
        #Get the new best step (depends on if we're inserting or deleting)
        if operator == Insert!
            findNextEquivClass!(newStep, dataParsed, g, findBestInsert, verbose)
        else
            findNextEquivClass!(newStep, dataParsed, g, findBestDelete, verbose)
        end
        
        
        #If the score did not improve...
        if newStep.Δscore ≤ zero(A)
            #...stop searching
            break
        end
        
        #Use the insert or delete operator update the graph
        operator(g, newStep)
        
        #Covert the PDAG to a complete PDAG
        #undirect all edges unless they participate in a v-structure
        graphVStructure!(g)
        #Apply the 4 Meek rules to orient some edges in the graph
        meekRules!(g)

        #Reset the score
        newStep.Δscore = zero(A)

        
    end

    return nothing
end


function findNextEquivClass!(newStep, dataParsed, g, findBestOperation::Function, verbose)

    #Parallelization with synced threads (to avoid race conditions)
    @sync begin
        #Loop through all possible node combinations
        for (x,y) in allpairs(vertices(g)) 
            #Spawn a new thread
            Threads.@spawn begin
                #Check if there is an edge between two vertices
                hasEdge = isadjacent(g,x,y)

                #Skip over adjacent edges if we're trying to insert and vice versa
                if findBestOperation == findBestInsert ? !hasEdge : hasEdge
                    #For this pair of nodes, the following function will:
                        #(1) check if a valid operator exists
                        #(2) score all valid operators and return the one with the highest score
                        #findBestOperation is either "findBestInsert" or "findBestDelete"
                    (bestSubset, bestScore) = findBestOperation(dataParsed, g, x, y)

                    #If the best valid operator was better than any previously found...
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
    end
end


function findBestInsert(dataParsed::ParseData{Matrix{A}}, g, x, y) where A
    #Calculate two (possibly empty) sets of nodes
    # NAxy: any nodes that are undirected neighbors of y and connected to x by any edge
    # Txy: any subset of the undirected neighbors of y not connected to x
    Tyx = calcT(g, y, x)
    NAyx = calcNAyx(g, y, x)


    #Ceate two containers to hold the best found score and best subset of Tyx
    bestScore = zero(A)
    bestT = Vector{Int}()

    #Keep a list of invalid sets
    invalid = Vector{Vector{Int}}()
    
    #Loop through all possible subsets of Tyx
    for T in powerset(Tyx)
        if checkSupersets(T,invalid)
            NAyxT = NAyx ∪ T
            if isclique(g, NAyxT) && isblocked(g, y, x, NAyxT)

                #Score the valid Insert
                PAy = parents(g,y)
                PAy⁺ = PAy ∪ x
                newScore = score(dataParsed, NAyxT ∪ PAy⁺, y) - score(dataParsed, NAyxT ∪ PAy, y)
                
                #Save the new score if it was better than any previous
                if newScore > bestScore
                    bestT = T
                    bestScore = newScore
                end
            end
        else
            #Record that the subset T is invalid
            push!(invalid,T)
        end
    end

    return (bestT, bestScore)
end

#Check if the set T is contained in any invalid set
function checkSupersets(T,invalid)
    for i ∈ invalid
        if i ⊆ T
            return false
        end
    end
    return true
end


function findBestDelete(dataParsed::ParseData{Matrix{A}}, g, x, y) where A
    #Calculate two (possibly empty) sets of nodes
    # NAxy: any nodes that are undirected neighbors of y and connected to x by any edge
    # Hyx: any subset of the undirected neighbors of y that are connected to x
    NAyx = calcNAyx(g, y, x)
    Hyx = NAyx

    #Ceate two containers to hold the best found score and best subset of Tyx
    bestScore = zero(A)
    bestH = Vector{Int}()

    #Loop through all possible subsets of Hyx
    for H in powerset(Hyx)
        #Calculate NAyx \ {H}
        NAyx_H = setdiff(NAyx,H)

        #Check if the operator is valid
        if isclique(g, NAyx_H)

            #Score the valid operator 
            PAy = parents(g,y)
            PAy⁻ = setdiff(PAy, x)
            newScore = score(dataParsed, NAyx_H ∪ PAy⁻, y) - score(dataParsed, NAyx_H ∪ PAy, y)
            
            if newScore > bestScore
                bestH = H
                bestScore = newScore
            end
        end
    end

    return (bestH, bestScore)
end


####################################################################
# Scoring function
####################################################################

@memoize LRU(maxsize=1_000_000) function score(dataParsed::ParseData{Matrix{A}}, nodeParents, node) where A

    #Unpack some variables from the dataParsed structure
    n = A(dataParsed.numObservations) #convert datatype
    data = dataParsed.normAugData
    p = dataParsed.penalty

    #The last column of dataParsed.normAugData is all ones which is required for a linear regression with an intercept. If there are no node parents, model the child node with just the intercept, else use the parents and the intercept
    if isempty(nodeParents)
        parentsAndIncept = [dataParsed.numFeatures+1]
    else
        parentsAndIncept = [nodeParents; dataParsed.numFeatures+1]
    end

    #To calculate the score we need a mean-squared error which we can get by regessing the the child node onto the parents
    #Create variables for a linear regression y = X*b

    #Use views to avoid creating copies of the data in memory
        # X is the design matrix, augmented with a column of ones at the end
        # X is also been standardized so mean(columns)=0 and std(column)=1
        # y is data from the child node being tested
    @views begin
        y = data[:,node]
        X = data[:,parentsAndIncept]
    end

    #Perform a linear regression
    b = X \ y

    #Get the estimation
    ŷ = X*b

    #Next we want to calculate the log likelihood that these are the parents of our node
    # score = log(P(data|Model)) ≈ -BIC/2
    # because we're only comparing log likelihoods we'll ignore the 1/2 factor
    # when P(⋅) is Guassian, log(P(data|Model)) takes the form:
    # -k⋅log(n) - n⋅log(mse)
    # k is the number of free parameters and mse is mean squared error
    k = length(parentsAndIncept) #includes the intercept
    mse = sum(x->x^2, y-ŷ) / n

    #return the final score we want to maximize (which is technically -2BIC)
    return -p*k*log(n) - n*log(mse)
end



####################################################################
# Using Meek's Rules to Update PDAG
####################################################################

#Revert a graph to undirected edges and unshielded colliders (i.e. parents not adjacent)
function graphVStructure!(g)
    
    #loop through all vertices
    for x in vertices(g)
        #get the parents of current vertex
        parentsX = parents(g,x)
        #if all the parents are adjacent, undirect the edges
        #(if there is only 1 parent, then it is still a clique)
        if isclique(g, parentsX)
            for p in parentsX
                add_edge!(g, x, p)
                add_edge!(g, p, x)
            end
        end
    end
end


#Thoughts on improvement
    #Would it be cleaner to find all unique triples where one edge is undirected?
    #Can this update be more local? Most edges will remain unchanged.

#This is set up so we only loop through all the vertices and edges once
function meekRules!(g)
    
    #Loop through all the edges in the graph (x-y)
    for edge in edges(g)
        #We only need to update undirected edges
        if !isoriented(g, edge)

            #For clarity extract the edge vertices
            (x, y) = src(edge), dst(edge)

            #Label the neighbors of the nodes that comprise the current edge
            #Neighbors can be labelled "parent", "child", or "undirected"
            #e.g. v₁→x-v₂, v₁ is a parent and v₂ is undirected
            xNeighbors, yNeighbors =  categorizeNeighbors(g, x, y)

            #Rule 1: If x or y have a unique parent, direct edge away from it
            #Check if x has a unique parent 
            if R1(xNeighbors, yNeighbors)
                #Add x→y
                orientedge!(g, x, y)
            #Check if y has a unique parent
            elseif R1(yNeighbors, xNeighbors)
                #Add y→x
                orientedge!(g, y, x)
            end

            #Rule 2: Direct the edge away from a potential cycle
            #Check if x has a child that is a parent of y
            if R2(xNeighbors, yNeighbors)
                #Add x→y
                orientedge!(g, x, y)
            #Check if y has a child that is a parent of x
            elseif R2(yNeighbors, xNeighbors)
                #Add y→x
                orientedge!(g, y, x)
            end

            #Rule 3: Double Triangle, diagonal
            if R3(xNeighbors, yNeighbors)
                #Add x→y
                orientedge!(g, x, y)
            elseif R3(yNeighbors, xNeighbors)
                #Add y→x
                orientedge!(g, y, x)
            end

            #Rule 4: Double Triangle, side
            #Requires an additional test that can't be check with the neighbor dictionaries 
            #So, we need to pass the graph through as well
            if R4(xNeighbors, yNeighbors, g)
                #Add x→y
                orientedge!(g, x, y)
            elseif R4(yNeighbors, xNeighbors, g)
                #Add y→x
                orientedge!(g, y, x)
            end
            
        end
    end

    return nothing
end


function R1(neighborSet1, neighborSet2)
    #given x-y, look for patterns that match v₁→x and not(v₁→y)
    for (v₁, category) in neighborSet1
        if category == :parent && !haskey(neighborSet2,v₁)
            return true
        end
    end
    return false
end

function R2(neighborSet1, neighborSet2)
    #given x-y, look for patterns that match x→v₁→y
    for (v₁, category) in neighborSet1
        if category == :child && haskey(neighborSet2,v₁) && neighborSet2[v₁]==:parent
            return true
        end
    end
    return false
end

function R3(neighborSet1, neighborSet2)
        
    #given x-y, count the number of patterns that match x-v₁→y
    counter = 0
    for (v₁, category) in neighborSet1
        if category == :undirected && haskey(neighborSet2,v₁) && neighborSet2[v₁]==:parent
            counter += 1
            #If we find two such paths, the Rule 3 pattern is matched
            if counter == 2
                return true
            end
        end
    end

    return false
end

function R4(neighborSet1, neighborSet2, g)
        
    #given x-y, look for patterns that match x-v₁→v₂→y and x-v₂
    for (v₂, category) in neighborSet1
        #first look for x-v₂→y
        if category == :undirected && haskey(neighborSet2,v₂) && neighborSet2[v₂]==:parent
            #check for x-v₁ and v₁→v₂
            for (v₁, category) in neighborSet1
                if category == :undirected && isparent(g, v₁, v₂)
                    return true
                end
            end
        end 
    end

    return false
end


#Categorize neighboring edges as "parent", "child", or "undirected"
function categorizeNeighbors(g, x, y)
    
    #Create a dictionary of neighbors (e.g. vertex 5 => :parent)
    xNeighbors = Dict{Int, Symbol}()
    yNeighbors = Dict{Int, Symbol}()

    #Loop through all the vertices in the graph
    for vᵢ in vertices(g)

        #If vᵢ is adjacent, give it a category
        if isadjacent(g, x, vᵢ)
            xNeighbors[vᵢ] = setCategory(g, vᵢ, x)
        end

        #Same procedure to the other vertex in the current edge
        if isadjacent(g, y, vᵢ)
            yNeighbors[vᵢ] = setCategory(g, vᵢ, y)
        end
    end

    return (xNeighbors, yNeighbors)
end


function setCategory(g, v₁, v₂)

    #Test if v₁ → v₂
    if isparent(g,v₁,v₂)
        return :parent
    #Test if v₁ ← v₂
    elseif ischild(g,v₁,v₂)
        return :child
    #If not then we must have an undirected edge
    else
        return :undirected
    end
end