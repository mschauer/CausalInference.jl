using Graphs

function gensearch(G, S, rules)
    visited = falses(2*nv(G))
    result = falses(nv(G))

    function genvisit(G, v, prevedge)
        prevedge != -1 && (visited[2*v-prevedge] = true)
        for nextedge in [0, 1]
            nextedge == 0 && (neighbors = outneighbors(G, v))
            nextedge == 1 && (neighbors = inneighbors(G, v))
            for w in neighbors
                (cont, yld) = rules(prevedge, nextedge, v, w)
                yld && (result[w] = true)
                cont && !visited[2*w-nextedge] && genvisit(G, w, nextedge)
            end
        end
    end

    for v in S
        genvisit(G, v, -1) 
    end
    return result
end


