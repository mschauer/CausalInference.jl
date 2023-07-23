


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