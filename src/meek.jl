
# This is set up so we only loop through all the vertices and edges once
function meekRules!(g)
    
    # Loop through all the edges in the graph (x-y)
    for edge in edges(g)
        # We only need to update undirected edges
        if !isoriented(g, edge)

            # For clarity extract the edge vertices
            (x, y) = src(edge), dst(edge)

            # Label the neighbors of the nodes that comprise the current edge
            # Neighbors can be labelled "parent", "child", or "undirected"
            # e.g. v₁→x-v₂, v₁ is a parent and v₂ is undirected
            xNeighbors, yNeighbors =  categorizeNeighbors(g, x, y)

            # Rule 1: If x or y have a unique parent, direct edge away from it
            # Check if x has a unique parent 
            if meek_R1(xNeighbors, yNeighbors)
                #Add x→y
                orientedge!(g, x, y)
            # Check if y has a unique parent
            elseif meek_R1(yNeighbors, xNeighbors)
                # Add y→x
                orientedge!(g, y, x)
            end

            # Rule 2: Direct the edge away from a potential cycle
            # Check if x has a child that is a parent of y
            if meek_R2(xNeighbors, yNeighbors)
                #Add x→y
                orientedge!(g, x, y)
            # Check if y has a child that is a parent of x
            elseif meek_R2(yNeighbors, xNeighbors)
                # Add y→x
                orientedge!(g, y, x)
            end

            # Rule 3: Double Triangle, diagonal
            if meek_R3(xNeighbors, yNeighbors)
                # Add x→y
                orientedge!(g, x, y)
            elseif meek_R3(yNeighbors, xNeighbors)
                # Add y→x
                orientedge!(g, y, x)
            end

            # Rule 4: Double Triangle, side
            # Requires an additional test that can't be check with the neighbor dictionaries 
            # So, we need to pass the graph through as well
            if meek_R4(xNeighbors, yNeighbors, g)
                # Add x→y
                orientedge!(g, x, y)
            elseif meek_R4(yNeighbors, xNeighbors, g)
                # Add y→x
                orientedge!(g, y, x)
            end
            
        end
    end

    return nothing
end


function meek_R1(neighborSet1, neighborSet2)
    # given x-y, look for patterns that match v₁→x and not(v₁→y)
    for (v₁, category) in neighborSet1
        if category == :parent && !haskey(neighborSet2,v₁)
            return true
        end
    end
    return false
end

function meek_R2(neighborSet1, neighborSet2)
    # given x-y, look for patterns that match x→v₁→y
    for (v₁, category) in neighborSet1
        if category == :child && haskey(neighborSet2,v₁) && neighborSet2[v₁]==:parent
            return true
        end
    end
    return false
end

function meek_R3(neighborSet1, neighborSet2)
        
    # given x-y, count the number of patterns that match x-v₁→y
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

function meek_R4(neighborSet1, neighborSet2, g)
        
    # given x-y, look for patterns that match x-v₁→v₂→y and x-v₂
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

