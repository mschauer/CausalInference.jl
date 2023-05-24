using Graphs

# helper function
function inunion(w, sets...)
    for S in sets
        w in S && return true
    end
    return false
end

function gensearch(g, S, rules)
    visited = falses(2*nv(g))
    result = falses(nv(g))

    function genvisit(g, v, prevedge)
        prevedge != -1 && (visited[2*v-prevedge] = true)
        for nextedge in [0, 1]
            nextedge == 0 && (neighbors = inneighbors(g, v))
            nextedge == 1 && (neighbors = outneighbors(g, v))
            for w in neighbors
                (cont, yld) = rules(prevedge, nextedge, v, w)
                yld && (result[w] = true)
                cont && !visited[2*w-nextedge] && genvisit(g, w, nextedge)
            end
        end
    end

    for v in S
        genvisit(g, v, -1) 
    end
    return Set{Integer}(findall(any, result))
end

function ancestors(g, X, inrem = Set{Integer}(), outrem = Set{Integer}())
    function rules(prevedge, nextedge, v, w)
        nextedge == 0 && !(v in inrem) && !(w in outrem) && return true, true
        return false, false
    end
    return gensearch(g, X, rules)
end

function descendants(g, X, inrem = Set{Integer}(), outrem = Set{Integer}())
    function rules(prevedge, nextedge, v, w)
        nextedge == 1 && !(v in outrem) && !(w in inrem) && return true, true
        return false, false
    end
    return gensearch(g, X, rules)
end

function alt_test_dsep(g, X, Y, S, inrem = Set{Integer}(), outrem = Set{Integer}())
    function rules(prevedge, nextedge, v, w)
       (v in inrem || w in outrem) && nextedge == 1 && return false, false 
       (v in outrem || w in inrem) && nextedge == 0 && return false, false
       if prevedge == -1 || (v in S && prevedge == 1 && nextedge == 0) || (!(v in S) && !(prevedge == 1 && nextedge ==0))
           return true, true
       end
       return false, false
    end
    R = gensearch(g, X, rules)
    return length(intersect(R, Y)) == 0 
end 

function alt_test_backdoor(g, X, Y, S)
    function rules(prevedge, nextedge, v, w)
        if (prevedge == -1 && nextedge == 0) || (prevedge != -1 && nextedge == 1 && !(inunion, X, S)) || (prevedge == 0 && nextedge == 0 && !(v in S)) || (prevedge == 1 && nextedge == 0 && v in S)
            return true, true
        end
        return false, false
    end
    R = gensearch(g, X, rules)
    return length(intersect(R, Y)) == 0 
end

function find_dsep(g, X, Y, I, R, inrem = Set{Integer}(), outrem = Set{Integer}())
    
end

function find_min_dsep(g, X, Y, I, R, inrem = Set{Integer}(), outrem = Set{Integer}())

end

function pcp(g, X, Y)
    DeX = descendants(g, X, X, Set{Integer}())
    AnY = ancestors(g, Y, Set{Integer}(), X)
    return intersect(setdiff(DeX, X), AnY)
end

function find_covariate_adjustment(g, X, Y, I, R)
    PCPXY = pcp(g, X, Y)
    DpcpXY = descendants(g, PCPXY)
    Z = setdiff(intersect(ancestors(g, union(X, Y, I)), R), union(X, Y, DpcpXY))
    if alt_dsep(g, X, Y, Z, PCPXY, Set{Integer}())
        return Z
    else
        return false
    end
end

function find_backdoor_adjustment(g, X, Y, I, R)
    Z = find_covariate_adjustment(g, X, Y, I, R)
    DeX = descendants(g, X, X, Set{Integer}())
    bdZ = setdiff(Z, DeX)
    if alt_dsep(g, X, Y, bdZ, Set{Integer}(), X)
        return Z
    else
        return false
    end
end

function find_frontdoor_adjustment(g, X, Y, I, R)

end

function find_min_frontdoor_adjustment(g, X, Y, I, R)

end

# TODO: write iterators with state as call stack
function list_covariate_adjustment(g, X, Y, I, R)

end

function list_backdoor_adjustment(g, X, Y, I, R)

end

function list_frontdoor_adjustment(g, X, Y, I, R)

end
