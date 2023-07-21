####################################################################
# Neighborhood functions 
####################################################################

# Should these be non-allocating iterators? 
neighbors_general(f, g, x) = Iterators.filter(vᵢ -> f(g, vᵢ, x), vertices(g))

# Note: inneighbors and outneighbors are already defined in Graphs.jl

descendent(g, x) = neighbors_general(isdescendent, g, x)


####################################################################
# Counting Functions
####################################################################

#Very similar to the neighborsGeneral function except we only keep a count of the valid nodes
function countGeneral(isFunction::Function, g, x)
    counter=0
    for vᵢ in vertices(g)
        if isFunction(g,vᵢ,x)
            counter += 1
        end
    end
    return counter
end

#All vertices connected to x
countNeighbors(g, x) = countGeneral(isadjacent, g, x)

#All vertices pointing to x
countNeighbors_in(g, x) = countGeneral(isparent, g, x)

#All vertices pointing away from x
countNeighbors_out(g, x) = countGeneral(ischild, g, x)

#All vertices connected to x by an undirected edge
countNeighbors_undirected(g, x) = countGeneral(isundirected, g, x)

#these are just aliases for the functions above
countParents(g, x) = countNeighbors_in(g, x)
countChildren(g, x) = countNeighbors_out(g, x)

####################################################################
# all* functions
####################################################################

#Generate an iterator equivalent to combinations(x,2) but produces a vector of tuples
allpairs(v) = Iterators.filter(i -> isless(i...), Iterators.product(v,v))

#Get every undirected edge in the entire graph g
allundirected(g) = [Edge(nodePair) for nodePair in allpairs(vertices(g)) if isundirected(g,nodePair)]

#Get every directed edge in the entire graph g
alldirected(g) = [Edge(nodePair) for nodePair in allpairs(vertices(g)) if isoriented(g,nodePair)]

#combine and sort all directed and undirected edges
alledges(g) = [Edge(nodePair) for nodePair in allpairs(vertices(g)) if isadjacent(g,nodePair)]




####################################################################
# Misc
####################################################################

# for x-y, get undirected neighbors of y and any neighbor of x
#calcNAyx(g, y::Integer, x::Integer) = setdiff(neighbors_undirect(g,y) ∩ neighbors(g,x), x)
function calcNAyx(g, y, x)
    neighborList = eltype(g)[]

    #loop through all vertices except for x
    for vᵢ in Iterators.filter(!isequal(x),vertices(g))
        #look for neighbors of y and adjacencies to x
        if isundirected(g,vᵢ,y) && isadjacent(g,vᵢ,x)
            push!(neighborList,vᵢ)
        end
    end

    return neighborList
end
calcNAyx(g, edge) = calcNAyx(g, dst(edge), src(edge))

#for x-y, undirected neighbors of y not connected to x
#calcT(g, y::Integer, x::Integer) = setdiff(neighbors_undirect(g,y), neighbors(g,x), x)
function calcT(g, y, x)
    neighborList = Int64[]

    #loop through all vertices except for x
    for vᵢ in Iterators.filter(!isequal(x),vertices(g))
        #look for neighbors of y and adjacencies to x
        if isundirected(g,vᵢ,y) && !isadjacent(g,vᵢ,x)
            push!(neighborList,vᵢ)
        end
    end

    return neighborList
end
calcT(g, edge) = calcT(g, dst(edge), src(edge))