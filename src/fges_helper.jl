####################################################################
# is* functions
####################################################################

"""
    isundirected(g, edge::Edge)
    isundirected(g, x, y)
Test if `x` and `y` are connected by a undirected edge in the graph `g`.
Can also perform the same test given an `edge`.
"""
isundirected(g, x, y) = has_edge(g,x,y) && has_edge(g,y,x)
isundirected(g, edge) = has_edge(g,edge) && has_edge(g,reverse(edge))

"""
    isadjacent(g, edge::Edge)
    isadjacent(g, x, y)
Test if `x` and `y` are connected by a any edge in the graph `g`.
Can also perform the same test given an `edge`.
"""
isadjacent(g, edge) = has_edge(g,edge) || has_edge(g,reverse(edge))
#isadjacent(g, x, y) = has_edge(g,x,y) || has_edge(g,y,x) #defined in pc.jl

"""
    isparent(g, edge::Edge)
    isparent(g, x, y)
Test if `x` is a parent of `y` in the graph `g`, meaning x→y.
Can also perform the same test given an `edge`.
"""
isparent(g, edge) = has_edge(g,edge) && !has_edge(g,reverse(edge))
isparent(g, x, y) = has_edge(g,x,y) && !has_edge(g,y,x)


"""
    ischild(g, edge::Edge)
    ischild(g, x, y)
Test if `x` is a child of `y` in the graph `g`, meaning x←y.
Can also perform the same test given an `edge`.
"""
ischild(g, edge) = !has_edge(g,edge) && has_edge(g,reverse(edge))
ischild(g, x, y) = !has_edge(g,x,y) && has_edge(g,y,x)


"""
    isdescendent(g, edge::Edge)
    isdescendent(g, x, y)
Return `true` if `x`←`y` OR `x`-`y` in the graph `g`.
"""
isdescendent(g, edge) = has_edge(g,reverse(edge))
isdescendent(g, x, y) = has_edge(g,y,x)

"""
    isoriented(g, edge::Edge)
    isoriented(g, x, y)
Test if `x` and `y` are connected by a directed edge in the graph `g`, either x←y OR x→y.
Can also perform the same test given an `edge`.
"""
isoriented(g, edge) = has_edge(g,edge) ⊻ has_edge(g,reverse(edge)) # xor
isoriented(g, x, y) = has_edge(g,x,y) ⊻ has_edge(g,y,x)


"""
    isclique(g, nodes)
Return `true` if all vertices in `nodes` are connected to each other in the graph `g`.
"""
function isclique(g, nodes)

    for nodePair in allpairs(nodes)
        if !isadjacent(g, nodePair)
            return false
        end
    end

    return true
end


#Do the nodesRemoved block all semi-directed paths between src and dest?
"""
    isblocked(g, src, dest, nodesRemoved)
Return `true` if there is no semi-directed path between `src` and `dest` in the graph `g`.
A set of vertices (`nodesRemoved`) can be removed from the graph before searching for a semi-directed path.
A semi-directed path between `src` and `dest` is a list of edges in `g` where every edge is either undirected or points toward `dest`. 
    src → x₁ - x₂ → dest ✓
    src → x₁ ← x₂ - dest ✖
"""
function isblocked(g, x, y, nodesRemoved)

    #Keep track of all the nodes visited
    visited = zeros(Bool, nv(g))

    # mark excluded vertices as visited
    for vᵢ in nodesRemoved 
        visited[vᵢ] = true
    end

    #If src or dest were in nodesRemoved, the path is blocked
    (visited[x] || visited[y]) && return true

    #if the src and dest are the same, path is itself
    x == y && return false

    #check if scr or dest have no neighbors
    if countNeighbors(g, x)==0 || countNeighbors(g, y)==0
        return true
    end

    queue = [x]
    visited[x] = true
    while !isempty(queue)
        currentNode = popfirst!(queue) # get new element from queue
        for vᵢ in descendent(g,currentNode)
            vᵢ == y && return false
            if !visited[vᵢ]
                push!(queue, vᵢ) # push onto queue
                visited[vᵢ] = true
            end
        end
    end
    return true
end
# function isblocked(g, x, y, nodesRemoved)
    
#     #If x or y are in nodesRemoved, the path is blocked
#     (x ∈ nodesRemoved || y ∈ nodesRemoved) && return true

#     #if the src and dest are the same, path is itself
#     x == y && return false

#     #Check out algorithm 3: https://leetcode.com/problems/find-if-path-exists-in-graph/solutions/2715942/official-solution/
#     edgeSet = IntDisjointSets(nv(g))

#     for edge in edges(g)
#         (source, destination) = src(edge), dst(edge)

#         if (source ∉ nodesRemoved) && (destination ∉ nodesRemoved)
#             union!(edgeSet, source, destination)
#         end
#     end

#     return !in_same_set(edgeSet, x, y)
# end


####################################################################
# Modify Edges
####################################################################

"""
    orientedge!(g, x, y)
Update the edge `x`-`y` to `x`→`y` in the graph `g`. 
"""
orientedge!(g, edge) = rem_edge!(g, reverse(edge))
orientedge!(g, x, y) = rem_edge!(g, y, x)



####################################################################
# Neighborhood functions 
####################################################################

#Should these be non-allocating iterators? 
neighborsGeneral(isFunction::Function, g, x) = Iterators.filter(vᵢ -> isFunction(g,vᵢ,x), vertices(g))

#Note: inneighbors and outneighbors are already defined in Graphs.jl

#All vertices connected to x by an undirected edge
neighbors_undirect(g, x) = neighborsGeneral(isundirected, g, x)

#these are just aliases for the functions above
parents(g, x) = inneighbors(g, x)
children(g, x) = outneighbors(g, x)

descendent(g, x) = neighborsGeneral(isdescendent, g, x)


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
countNeighbors_undirect(g, x) = countGeneral(isundirected, g, x)

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
    neighborList = Int64[]

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