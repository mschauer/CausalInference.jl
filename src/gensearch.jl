using Graphs

# helper function
function inunion(w, sets...)
    for S in sets
        w in S && return true
    end
    return false
end

function gensearch(g, S, pass)
    visited = falses(3*nv(g))
    function genvisit(g, v, prevedge)
        visited[3*v - 1 + prevedge] = true
        for nextedge in [0, 1]
            nextedge == 0 && (neighbors = inneighbors(g, v))
            nextedge == 1 && (neighbors = outneighbors(g, v))
            for w in neighbors
                pass(prevedge, nextedge, v, w) && !visited[3*w - 1 + nextedge] && genvisit(g, w, nextedge)
            end
        end
    end
    foreach(s -> genvisit(g, s, -1), S)
    return Set{Integer}([x for x in 1:nv(g) if visited[3*x-2] || visited[3*x-1] || visited[3*x]])
end

function ancestors(g, X, inrem = Set{Integer}(), outrem = Set{Integer}())
    function pass(prevedge, nextedge, v, w)
        nextedge == 0 && !(v in inrem) && !(w in outrem) && return true
        return false
    end
    return gensearch(g, X, pass)
end

function descendants(g, X, inrem = Set{Integer}(), outrem = Set{Integer}())
    function rules(prevedge, nextedge, v, w)
        nextedge == 1 && !(v in outrem) && !(w in inrem) && return true
        return false
    end
    return gensearch(g, X, rules)
end

function bayesball(g, X, S, inrem = Set{Integer}(), outrem = Set{Integer}())
    function pass(prevedge, nextedge, v, w)
       (v in inrem || w in outrem) && nextedge == 0 && return false
       (v in outrem || w in inrem) && nextedge == 1 && return false
       if prevedge == -1 || (v in S && prevedge == 1 && nextedge == 0) || (!(v in S) && !(prevedge == 1 && nextedge ==0))
           return true
       end
       return false
    end
    return gensearch(g, X, pass)
end

function alt_test_dsep(g, X, Y, S, inrem = Set{Integer}(), outrem = Set{Integer}())
    return length(intersect(bayesball(g, X, S, inrem, outrem), Y)) == 0 
end 

function alt_test_backdoor(g, X, Y, S)
    return length(intersect(bayesball(g, X, S, Set{Integer}(), X), Y)) == 0 
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
    function pass(prevedge, nextedge, v, w)
        !(w in A) && return false
       (v in inrem || w in outrem) && nextedge == 0 && return false
       (v in outrem || w in inrem) && nextedge == 1 && return false
        if ((prevedge == 1 && nextedge == 1) || prevedge == 0) && v in Z
            return false
        end
        return true
    end
    return gensearch(g, X, pass)
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
    if issubset(I, Z) && alt_test_dsep(g, X, Y, Z, PCPXY, Set{Integer}())
        return Z
    else
        return false
    end
end

function find_backdoor_adjustment(g, X, Y, I, R)
    Z = find_covariate_adjustment(g, X, Y, I, R)
    DeX = descendants(g, X, X, Set{Integer}())
    bdZ = setdiff(Z, DeX)
    if issubset(I, bdZ) && alt_test_dsep(g, X, Y, bdZ, Set{Integer}(), X)
        return bdZ
    else
        return false
    end
end

# this is of course also minimal bd set!
function find_min_covariate_adjustment(g, X, Y, I, R)
    PCPXY = pcp(g, X, Y)
    DePCP = descendants(g, PCPXY)
    Rd = setdiff(R, DePCP) 
    # TODO: implement vetos!!!
    function vetos(prevedge, nextedge, v, w)
        v in X && w in PCPXY && return true
        return false
    end
    return find_min_dsep(g, X, Y, I, Rd, PCPXY, Set{Integer}())
end

function find_frontdoor_adjustment(g, X, Y, I, R)
    Za = bayesball(g, X, Set{Integer}(), Set{Integer}(), X)
    A = ancestors(g, Y)
    # could check for za \in Za not in I already here
    # paper is wrong here -> double check and correct
    function Zab_pass(prevedge, nextedge, v, w)
        nextedge == 1 && return true
        prevedge == 1 && v in A && !(w in Za) && return true
        (prevedge == -1 || prevedge == 0) && nextedge == 0 && w in Za && return true
        return false
    end
    Zab = setdiff(Za, gensearch(g, Za, Zab_pass))
    if issubset(I, Zab) && alt_test_dsep(g, X, Y, Zab, X, Zab)
        return Zab
    else
        return false
    end
end

# write explanation for these rules
function find_min_frontdoor_adjustment(g, X, Y, I, R)
    Zii = find_frontdoor_adjustment(g, X, Y, I, R)
    Zii == false && return false 
    function Za_rules(prevedge, nextedge, v, w)
        prevedge in [-1, 0] && nextedge == 0 && !(v in Zii) && !inunion(w, X, Y) && return true
        return false
    end
    Za = intersect(gensearch(g, Y, Za_rules), Zii) 
    function Zxy_rules(prevedge, nextedge, v, w)
        prevedge in [-1, 1] && nextedge == 1 && !(v in Za) && !inunion(w, X, Y, I) && return true
        return false
    end
    Zxy = intersect(gensearch(g, X, Zxy_rules), Za)
    function Zzy_rules(prevedge, nextedge, v, w)
        prevedge in [-1, 0] && nextedge == 0 && !inunion(w, X, I, Zxy) && return true
        prevedge in [0, 1] && nextedge == 1 && !(w in X) && !inunion(v, I, Za) && return true
        prevedge == 1 && nextedge == 0 && inunion(v, I, Za) && !inunion(w, X, I, Zxy) && return true
        return false
    end
    Zzy = intersect(gensearch(g, union(I, Zxy), Zzy_rules), Za) 
    return union(I, Zxy, Zzy)
end

struct AdjustmentIterator{T<:Integer, F<:Function}
    g::SimpleDiGraph{T}
    X::Set{T}
    Y::Set{T}
    I::Set{T}
    R::Set{T}
    find_adjustment::F
end

function downwards(state, I, R)
    v = first(setdiff(R, I))
    push!(state, ("up", "I", v, I, R))
    push!(state, ("down", "I", v, I, R))
    push!(state, ("up", "R", v, I, R))
    push!(state, ("down", "R", v, I, R))
end

# TODO: do this more elegantly
function Base.iterate(A::AdjustmentIterator)
    R = deepcopy(A.R)
    I = deepcopy(A.I)
    state = Vector{Tuple{String, String, Int64, Set{Int64}, Set{Int64}}}()
    !A.find_adjustment(A.g, A.X, A.Y, I, R) && return nothing
    issetequal(I, R) && return I, state
    downwards(state, I, R)
    Base.iterate(A, state)
end

function Base.iterate(A::AdjustmentIterator, state)
    while !isempty(state)
        dir, set, v, I, R = pop!(state)
        dir == "down" && set == "I" && push!(I, v)
        dir == "up" && set == "I" && delete!(I, v)
        dir == "down" && set == "R" && delete!(R, v)
        dir == "up" && set == "R" && push!(R, v)
        !A.find_adjustment(A.g,  A.X, A.Y, I, R) && continue
        issetequal(I, R) && return deepcopy(I), state
        if dir == "down"
            downwards(state, I, R)
        end
    end
    return nothing
end

function list_covariate_adjustment(g, X, Y, I, R)
    return AdjustmentIterator(g, X, Y, I, R, find_covariate_adjustment)
end

function list_backdoor_adjustment(g, X, Y, I, R)
    return AdjustmentIterator(g, X, Y, I, R, find_backdoor_adjustment)
end

function list_frontdoor_adjustment(g, X, Y, I, R)
    return AdjustmentIterator(g, X, Y, I, R, find_frontdoor_adjustment)
end
