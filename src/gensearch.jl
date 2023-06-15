using Graphs

## missing features:
# - make functions callable with non-sets (e.g. X = 1, Y = 4) -> what is the best way to do this?

const INIT = 0
const LEFT = 1
const RIGHT = 2

## helper functions

function inunion(w, sets...)
    for S in sets
        w in S && return true
    end
    return false
end

# maybe replace unused arguments by type like so:
# function no_veto(::T, ::T, ::T, ::T)
function no_veto(pe, ne, v, w)
    return false
end

function no_outgoing(S)
    return (pe, ne, v, w) -> v in S && ne == RIGHT 
end

function no_incoming(S)
    return (pe, ne, v, w) -> v in S && ne == LEFT
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
    foreach(s -> genvisit(g, s, INIT), S)
    return Set([x for x in 1:nv(g) if visited[3*x-2] || visited[3*x-1] || visited[3*x]])
end

"""
    ancestors(g, X, veto = no_veto)

Return the set of ancestors of the set of vertices `X` in graph `g`. 

That is, a vertex is in the output set if it is an ancestor of at least one of the vertices in `X`.
"""
function ancestors(g, X, veto = no_veto)
    return gensearch(g, X, (pe, ne, v, w) -> ne == LEFT && !veto(pe, ne, v, w))
end

"""
    descendants(g, X, veto = no_veto)

Return the set of descendants of the set of vertices `X` in graph `g`. 

That is, a vertex is in the output set if it is a descendant of at least one of the vertices in `X`.
"""
function descendants(g, X, veto = no_veto)
    return gensearch(g, X, (pe, ne, v, w) -> ne == RIGHT && !veto(pe, ne, v, w))
end

function bayesball(g, X, S = Set{eltype(g)}(), veto = no_veto)
    return gensearch(g, X, (pe, ne, v, w) -> !veto(pe, ne, v, w) && (pe == INIT || (v in S && pe == RIGHT && ne == LEFT) || (!(v in S) && !(pe == RIGHT && ne == LEFT))))
end

"""
    alt_test_dsep(g, X, Y, S, veto = no_veto)

Check if sets of vertices `X` and `Y` are d-separated in `g` given `S`. 

An alternative to the `test_dsep` function, which uses gensearch under the hood. Might be (a bit) slower. 
"""
function alt_test_dsep(g, X, Y, S, veto = no_veto)
    return length(intersect(bayesball(g, X, S, veto), Y)) == 0 
end 

"""
    test_covariate_adjustment(g, X, Y, S)

Check if `S` is a covariate adjustment set relative to `(X, Y)` in graph `g`. 

Done by checking the graphical criterion given in https://arxiv.org/abs/1203.3515 using the algorithmic approach proposed in https://arxiv.org/abs/1803.00116. Output is a boolean.
"""
function test_covariate_adjustment(g, X, Y, S) 
    PCPXY = pcp(g, X, Y)
    length(intersect(S, descendants(g, PCPXY))) != 0 && return false
    return alt_test_dsep(g, X, Y, S, (pe, ne, v, w) -> v in X && w in PCPXY && ne == RIGHT)
end

"""
    alt_test_backdoor(g, X, Y, S)

Check if `S` satisfies the backdoor criterion relative to `(X, Y)` in graph `g`. 

Follows the algorithmic approach proposed in https://arxiv.org/abs/1803.00116.
"""
function alt_test_backdoor(g, X, Y, S)
    length(intersect(S, descendants(g, X))) != 0 && return false
    return alt_test_dsep(g, X, Y, S, no_outgoing(X))
end

"""
    find_dsep(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y), veto = no_veto)

Find a d-separator `Z` with `I subseteq Z subseteq R` for sets of vertices `X` and `Y` in `g`, else return `false`. 

Follows the algorithmic approach proposed in https://arxiv.org/abs/1803.00116. 
"""
function find_dsep(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y), veto = no_veto)
    Z = intersect(R, setdiff(ancestors(g, union(X, Y, I), veto), X, Y))
    if alt_test_dsep(g, X, Y, Z, veto)
        return Z
    else
        return false
    end
end

function closure(g, X, A, Z, veto)
    return gensearch(g, X, (pe, ne, v, w) -> (w in A) && !veto(pe, ne, v, w) && !(((pe == RIGHT && ne == RIGHT) || pe == LEFT) && v in Z))
end

"""
    find_min_dsep(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y), veto = no_veto)

Find an inclusion minimal d-separator `Z` with `I subseteq Z subseteq R` for sets of vertices `X` and `Y` in `g`, i.e., one for which no subset is a d-separator, else return `false`.

Follows the algorithmic approach proposed in https://arxiv.org/abs/1803.00116. 
"""
function find_min_dsep(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y), veto = no_veto)
    A = ancestors(g, union(X, Y, I), veto)
    Z = find_dsep(g, X, Y, I, R, veto)
    Z == false && return false
    ZX = union(intersect(Z, closure(g, X, A, Z, veto)), I)
    return union(intersect(ZX, closure(g, Y, A, ZX, veto)), I)
end

function pcp(g, X, Y)
    return intersect(setdiff(descendants(g, X, no_incoming(X)), X), ancestors(g, Y, no_outgoing(X))
)
end

"""
    find_covariate_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))

Find a covariate adjustment set `Z` with `I subseteq Z subseteq R` for sets of vertices `X` and `Y` in `g`, else return `false`. 

Follows the algorithmic approach proposed in https://arxiv.org/abs/1803.00116. 
"""
function find_covariate_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    PCPXY = pcp(g, X, Y)
    Z = setdiff(intersect(ancestors(g, union(X, Y, I)), R), X, Y, descendants(g, PCPXY))
    if issubset(I, Z) && alt_test_dsep(g, X, Y, Z, (pe, ne, v, w) -> v in X && w in PCPXY && ne == RIGHT)
        return Z
    else
        return false
    end
end

"""
    find_backdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))

Find a backdoor adjustment set `Z` with `I subseteq Z subseteq R` for sets of vertices `X` and `Y` in `g`, else return `false`. 

Follows the algorithmic approach proposed in https://arxiv.org/abs/1803.00116. 
"""
function find_backdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    bdZ = setdiff(find_covariate_adjustment(g, X, Y, I, R), descendants(g, X, no_incoming(X)))
    if issubset(I, bdZ) && alt_test_dsep(g, X, Y, bdZ, no_outgoing(X))
        return bdZ
    else
        return false
    end
end

"""
    find_covariate_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))

Find an inclusion minimal covariate adjustment set `Z` with `I subseteq Z subseteq R` for sets of vertices `X` and `Y` in `g`, else return `false`. 

Follows the algorithmic approach proposed in https://arxiv.org/abs/1803.00116.  
"""
function find_min_covariate_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    PCPXY = pcp(g, X, Y)    
    # this is also minimal bd set!
    return find_min_dsep(g, X, Y, I, setdiff(R, descendants(g, PCPXY)), (pe, ne, v, w) -> v in X && w in PCPXY && ne == RIGHT)
end

"""
    find_frontdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))

Find a frontdoor adjustment set `Z` with `I subseteq Z subseteq R` for sets of vertices `X` and `Y` in `g`, else return `false`. 

Follows the algorithm given in  https://arxiv.org/abs/2211.16468. 
"""
function find_frontdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    Za = setdiff(R, bayesball(g, X, Set{Integer}(), no_outgoing(X)))
    A = ancestors(g, Y)
    function Zab_pass(pe, ne, v, w)
        v in X && return false
        ne == RIGHT && return true
        pe == RIGHT && v in A && !(w in Za) && return true
        (pe == INIT || pe == LEFT) && !(w in Za) && return true
        return false
    end
    Zab = setdiff(Za, gensearch(g, Y, Zab_pass))
    if issubset(I, Zab) && length(intersect(descendants(g, X, no_outgoing(Zab)), Y)) == 0
        return Zab
    else
        return false
    end
end

"""
    find_min_frontdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))

Find an inclusion minimal frontdoor adjustment set `Z` with `I subseteq Z subseteq R` for sets of vertices `X` and `Y` in `g`, else returns `false`. 

Follows the algorithm given in  https://arxiv.org/abs/2211.16468. 
"""
# write explanation for the pass rules
    function find_min_frontdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    Zii = find_frontdoor_adjustment(g, X, Y, I, R)
    Zii == false && return false 
    function Za_pass(pe, ne, v, w)
        pe in [INIT, LEFT] && ne == LEFT && !(v in Zii) && !inunion(w, X, Y) && return true
        return false
    end
    Za = intersect(gensearch(g, Y, Za_pass), Zii) 
    function Zxy_pass(pe, ne, v, w)
        pe in [INIT, RIGHT] && ne == RIGHT && !(v in Za) && !inunion(w, X, Y, I) && return true
        return false
    end
    Zxy = intersect(gensearch(g, X, Zxy_pass), Za)
    function Zzy_pass(pe, ne, v, w)
        pe in [INIT, LEFT] && ne == LEFT && !inunion(w, X, I, Zxy) && return true
        pe in [LEFT, RIGHT] && ne == RIGHT && !(w in X) && !inunion(v, I, Za) && return true
        pe == RIGHT && ne == LEFT && inunion(v, I, Za) && !inunion(w, X, I, Zxy) && return true
        return false
    end
    Zzy = intersect(gensearch(g, union(I, Zxy), Zzy_pass), Za) 
    return union(I, Zxy, Zzy)
end

struct ConstraintIterator{T<:Integer, F<:Function}
    g::SimpleDiGraph{T}
    X::Set{T}
    Y::Set{T}
    I::Set{T}
    R::Set{T}
    find::F
end

function downwards(state, I, R)
    v = first(setdiff(R, I))
    push!(state, (:up, :I, v, I, R))
    push!(state, (:down, :I, v, I, R))
    push!(state, (:up, :R, v, I, R))
    push!(state, (:down, :R, v, I, R))
end

# TODO: do this more elegantly
function Base.iterate(A::ConstraintIterator)
    R = deepcopy(A.R)
    I = deepcopy(A.I)
    state = Vector{Tuple{String, String, Int64, Set{Int64}, Set{Int64}}}()
    A.find(A.g, A.X, A.Y, I, R) == false && return nothing
    issetequal(I, R) && return I, state
    downwards(state, I, R)
    Base.iterate(A, state)
end

function Base.iterate(A::ConstraintIterator, state)
    while !isempty(state)
        dir, set, v, I, R = pop!(state)
        dir == :down && set == :I && push!(I, v)
        dir == :up && set == :I && delete!(I, v)
        dir == :down && set == :R && delete!(R, v)
        dir == :up && set == :R && push!(R, v)
        A.find(A.g,  A.X, A.Y, I, R) == false && continue
        issetequal(I, R) && return deepcopy(I), state
        if dir == :down
            downwards(state, I, R)
        end
    end
    return nothing
end

Base.IteratorSize(::ConstraintIterator) = Base.SizeUnknown()

"""
    list_dseps(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))

List all d-separators `Z` with `I subseteq Z subseteq R` for sets of vertices `X` and `Y` in `g`. 
"""
function list_dseps(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    return ConstraintIterator(g, X, Y, I, R, find_dsep)
end

"""
    list_covariate_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))

List all covariate adjustment sets `Z` with `I subseteq Z subseteq R` for sets of vertices `X` and `Y` in `g`. 
"""
function list_covariate_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    return ConstraintIterator(g, X, Y, I, R, find_covariate_adjustment)
end

"""
    list_backdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))

List all back-door adjustment sets `Z` with `I subseteq Z subseteq R` for sets of vertices `X` and `Y` in `g`. 
"""
function list_backdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    return ConstraintIterator(g, X, Y, I, R, find_backdoor_adjustment)
end

"""
    list_frontdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))

List all front-door adjustment sets `Z` with `I subseteq Z subseteq R` for sets of vertices `X` and `Y` in `g`. 
"""
function list_frontdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    return ConstraintIterator(g, X, Y, I, R, find_frontdoor_adjustment)
end

# note that this approach does *not* work for *minimal* dseps and adjustment sets
