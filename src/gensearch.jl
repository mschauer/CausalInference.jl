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

function bayesball(g, X, S, inrem = Set{Integer}(), outrem = Set{Integer}())
    function rules(prevedge, nextedge, v, w)
       (v in inrem || w in outrem) && nextedge == 1 && return false, false 
       (v in outrem || w in inrem) && nextedge == 0 && return false, false
       if prevedge == -1 || (v in S && prevedge == 1 && nextedge == 0) || (!(v in S) && !(prevedge == 1 && nextedge ==0))
           return true, true
       end
       return false, false
    end
    return gensearch(g, X, rules)
end

function alt_test_dsep(g, X, Y, S, inrem = Set{Integer}(), outrem = Set{Integer}())
    R = bayesball(g, X, S, inrem, outrem)
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
    Z = intersect(R, setdiff(ancestors(g, union(X, Y, I), inrem, outrem), union(X, Y)))
    if alt_test_dsep(g, X, Y, Z, inrem, outrem)
        return Z
    else
        return false
    end
end

function closure(g, X, A, Z, inrem, outrem)
    function rules(prevedge, nextedge, v, w)
        !(w in A) && return false, false
       (v in inrem || w in outrem) && nextedge == 1 && return false, false 
       (v in outrem || w in inrem) && nextedge == 0 && return false, false
        if ((prevedge == 1 && nextedge == 1) || prevedge == 0) && v in Z
            return false, false
        end
        return true, true
    end
    return gensearch(g, X, rules)
end

function find_min_dsep(g, X, Y, I, R, inrem = Set{Integer}(), outrem = Set{Integer}())
    A = ancestors(g, union(X, Y, I), inrem, outrem)
    Z = find_dsep(g, X, Y, I, R, inrem, outrem)
    Z == false && return false
    Xstar = closure(g, X, A, Z, inrem, outrem)
    ZX = union(intersect(Z, Xstar), I)
    Ystar = closure(g, Y, A, ZX, inrem, outrem)
    return union(intersect(ZX, Ystar), I)
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
    if alt_test_dsep(g, X, Y, Z, PCPXY, Set{Integer}())
        return Z
    else
        return false
    end
end

function find_backdoor_adjustment(g, X, Y, I, R)
    Z = find_covariate_adjustment(g, X, Y, I, R)
    DeX = descendants(g, X, X, Set{Integer}())
    bdZ = setdiff(Z, DeX)
    if alt_test_dsep(g, X, Y, bdZ, Set{Integer}(), X)
        return Z
    else
        return false
    end
end

# this is of course also minimal bd set!
function find_min_covariate_adjustment(g, X, Y, I, R)
    PCPXY = pcp(g, X, Y)
    DePCP = descendants(g, PCPXY)
    Rd = setdiff(R, DePCP) 
    return find_min_dsep(g, X, Y, I, Rd, PCPXY, Set{Integer}())
end

function find_frontdoor_adjustment(g, X, Y, I, R)
    Za = bayesball(g, X, Set{Integer}(), Set{Integer}(), X)
    A = ancestors(g, Y)
    # could check for za \in Za not in I already here
    # paper is wrong here -> double check and correct
    function Zabrules(prevedge, nextedge, v, w)
        nextedge == 1 && return true, true
        prevedge == 1 && v in A && !(w in Za) && return true, true 
        (prevedge == -1 || prevedge == 0) && nextedge == 0 && w in Za && return true, false
        return false, false
    end
    Zab = setdiff(Za, gensearch(g, Za, Zabrules))
    if issubset(I, Zab) && alt_test_dsep(g, X, Y, Zab, X, Zab)
        return Zab
    else
        return false
    end
end

# double check these rules
function find_min_frontdoor_adjustment(g, X, Y, I, R)
    Zii = find_frontdoor_adjustment(g, X, Y, I, R)
    Zii == false && return false 
    function Za_rules(prevedge, nextedge, v, w)
        cont, yld = false, false
        if prevedge in [-1, 0] && nextedge == 0
            !inunion(w, X, Y, Zii) && (cont = true)
            w in Zii && (yld = true)
        end
        return cont, yld
    end
    Za = gensearch(g, Y, Za_rules) 
    function Zxy_rules(prevedge, nextedge, v, w)
        cont, yld = false, false
        if prevedge in [-1, 1] && nextedge == 1
            !inunion(w, X, Y, I, Za) && (cont = true)
            w in Za && (yld = true)
        end
        return cont, yld
    end
    Zxy = gensearch(g, X, Zxy_rules)
    function Zzy_rules(prevedge, nextedge, v, w)
        cont, yld = false, false
        if prevedge in [-1, 0] && nextedge == 0
            !inunion(w, X, I, Zxy) && (cont = true)
            w in Za && (yld = true)
        end

        if prevedge in [0, 1] && nextedge == 1
            !(w in X) && !inunion(v, I, Za) && (cont = true)
            w in Za && !inunion(v, I, Za) && (yld = true)
        end

        if prevedge == 1 && nextedge == 0
            inunion(v, I, Za) && !inunion(w, X, I, Zxy) && (cont = true)
            inunion(v, I, Za) && w in Za && (yld = true)
        end
        return cont, yld
    end
    Zzy = gensearch(g, union(I, Zxy), Zzy_rules) 
    return union(I, Zxy, Zzy)
end

# TODO: write iterators with state as call stack
function list_covariate_adjustment(g, X, Y, I, R)

end

function list_backdoor_adjustment(g, X, Y, I, R)

end

function list_frontdoor_adjustment(g, X, Y, I, R)

end
