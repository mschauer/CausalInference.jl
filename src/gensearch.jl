using Graphs

const INIT = 0
const LEFT = 1
const RIGHT = 2

# helper functions
function inunion(w, sets...)
    for S in sets
        w in S && return true
    end
    return false
end

function no_veto(pe, ne, v, w)
    return false
end

function no_outgoing(S, pe, ne, v, w)
   return v in S && ne == RIGHT 
end

function no_incoming(S, pe, ne, v, w)
    return v in S && ne == LEFT
end

function gensearch(g, S, pass)
    visited = falses(3*nv(g))
    function genvisit(g, v, pe)
        visited[3*v - pe] = true
        for ne in [LEFT, RIGHT]
            ne == LEFT && (neighbors = inneighbors(g, v))
            ne == RIGHT && (neighbors = outneighbors(g, v))
            for w in neighbors
                pass(pe, ne, v, w) && !visited[3*w - ne] && genvisit(g, w, ne)
            end
        end
    end
    foreach(s -> genvisit(g, s, -1), S)
    return Set{Integer}([x for x in 1:nv(g) if visited[3*x-2] || visited[3*x-1] || visited[3*x]])
end

function ancestors(g, X, veto = no_veto)
    return gensearch(g, X, (pe, ne, v, w) -> ne == LEFT && !veto(pe, ne, v, w))
end

function descendants(g, X, veto = no_veto)
    return gensearch(g, X, (pe, ne, v, w) -> ne == RIGHT && !veto(pe, ne, v, w))
end

function bayesball(g, X, S, veto = no_veto)
    return gensearch(g, X, (pe, ne, v, w) -> !veto(pe, ne, v, w) && pe == INIT || (v in S && pe == RIGHT && ne == LEFT) || (!(v in S) && !(pe == RIGHT && ne == LEFT)))
end

function alt_test_dsep(g, X, Y, S, veto = no_veto)
    return length(intersect(bayesball(g, X, S, veto), Y)) == 0 
end 

function alt_test_backdoor(g, X, Y, S)
    return length(intersect(bayesball(g, X, S, (pe, ne, v,w) -> no_outgoing(X, pe, ne, v, w)), Y)) == 0
end

function find_dsep(g, X, Y, I, R, veto)
    Z = intersect(R, setdiff(ancestors(g, union(X, Y, I), veto), union(X, Y)))
    if alt_test_dsep(g, X, Y, Z, veto)
        return Z
    else
        return false
    end
end

function closure(g, X, A, Z, veto)
    return gensearch(g, X, (pe, ne, v, w) -> (w in A) && !veto(pe, ne, v, w) && !(((pe == 1 && ne == 1) || pe == 0) && v in Z))
end

function find_min_dsep(g, X, Y, I, R, veto)
    A = ancestors(g, union(X, Y, I), veto)
    Z = find_dsep(g, X, Y, I, R, veto)
    Z == false && return false
    Xstar = closure(g, X, A, Z, veto)
    ZX = union(intersect(Z, Xstar), I)
    Ystar = closure(g, Y, A, ZX, veto)
    return union(intersect(ZX, Ystar), I)
end

function pcp(g, X, Y)
    DeX = descendants(g, X, (pe, ne, v, w) -> no_outgoing(X, pe, ne, v, w))
    AnY = ancestors(g, Y, (pe, ne, v, w) -> no_incoming(X, pe, ne, v, w))
    return intersect(setdiff(DeX, X), AnY)
end

function find_covariate_adjustment(g, X, Y, I, R)
    PCPXY = pcp(g, X, Y)
    DpcpXY = descendants(g, PCPXY)
    Z = setdiff(intersect(ancestors(g, union(X, Y, I)), R), union(X, Y, DpcpXY))
    if issubset(I, Z) && alt_test_dsep(g, X, Y, Z, (pe, ne, v, w) -> v in X && w in PCPXY && ne == RIGHT)
        return Z
    else
        return false
    end
end

function find_backdoor_adjustment(g, X, Y, I, R)
    Z = find_covariate_adjustment(g, X, Y, I, R)
    DeX = descendants(g, X, (pe, ne, v, w) -> no_incoming(X, pe, ne, v, w))
    bdZ = setdiff(Z, DeX)
    if issubset(I, bdZ) && alt_test_dsep(g, X, Y, bdZ, (pe, ne, v, w) -> no_outgoing(X, pe, ne, v, w))
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
    return find_min_dsep(g, X, Y, I, Rd, (pe, ne, v, w) -> v in X && w in PCPXY && ne == RIGHT)
end

function find_frontdoor_adjustment(g, X, Y, I, R)
    Za = bayesball(g, X, Set{Integer}(), (pe, ne, v, w) -> no_outgoing(X, pe, ne, v, w))
    A = ancestors(g, Y)
    function Zab_pass(pe, ne, v, w)
        ne == 1 && return true
        pe == 1 && v in A && !(w in Za) && return true
        (pe == -1 || pe == 0) && ne == 0 && w in Za && return true
        return false
    end
    Zab = setdiff(Za, gensearch(g, Za, Zab_pass))
    if issubset(I, Zab) && alt_test_dsep(g, X, Y, Zab, (pe, ne, v, w) -> no_incoming(X, pe, ne, v, w) || no_outgoing(Zab, pe, ne, v, w))
        return Zab
    else
        return false
    end
end

# write explanation for these rules
function find_min_frontdoor_adjustment(g, X, Y, I, R)
    Zii = find_frontdoor_adjustment(g, X, Y, I, R)
    Zii == false && return false 
    function Za_rules(pe, ne, v, w)
        pe in [-1, 0] && ne == 0 && !(v in Zii) && !inunion(w, X, Y) && return true
        return false
    end
    Za = intersect(gensearch(g, Y, Za_rules), Zii) 
    function Zxy_rules(pe, ne, v, w)
        pe in [-1, 1] && ne == 1 && !(v in Za) && !inunion(w, X, Y, I) && return true
        return false
    end
    Zxy = intersect(gensearch(g, X, Zxy_rules), Za)
    function Zzy_rules(pe, ne, v, w)
        pe in [-1, 0] && ne == 0 && !inunion(w, X, I, Zxy) && return true
        pe in [0, 1] && ne == 1 && !(w in X) && !inunion(v, I, Za) && return true
        pe == 1 && ne == 0 && inunion(v, I, Za) && !inunion(w, X, I, Zxy) && return true
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
