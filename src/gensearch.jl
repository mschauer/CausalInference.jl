using Graphs

const INIT = 1
const LEFT = 2
const RIGHT = 3

## helper functions

function inunion(w, sets...)
    for S in sets
        w in S && return true
    end
    return false
end

function toset(object)
    if !isa(object, Set) 
        return Set(object)
    else 
        return object
    end
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
    visited = falses(3, nv(g))
    function genvisit(g, v, pe)
        visited[pe, v] = true
        for ne in (LEFT, RIGHT)
            ne == LEFT && (neighbors = inneighbors(g, v))
            ne == RIGHT && (neighbors = outneighbors(g, v))
            for w in neighbors
                pass(pe, ne, v, w) && !visited[ne, w] && genvisit(g, w, ne)
            end
        end
    end
    foreach(s -> genvisit(g, s, INIT), S)
    return Set(getindex.(findall(any(visited, dims=1)),2))
end

"""
    ancestors(g, X, veto = no_veto)

Return the set of ancestors of the set of vertices `X` in graph `g`. 

Every vertex is an ancestor of itself.
"""
function ancestors(g, X, veto = no_veto)
    X = toset(X)
    return gensearch(g, X, (pe, ne, v, w) -> ne == LEFT && !veto(pe, ne, v, w))
end

"""
    descendants(g, X, veto = no_veto)

Return the set of descendants of the set of vertices `X` in graph `g`. 

Every vertex is a descendant of itself.
"""
function descendants(g, X, veto = no_veto)
    X = toset(X)
    return gensearch(g, X, (pe, ne, v, w) -> ne == RIGHT && !veto(pe, ne, v, w))
end

"""
    bayesball(g, X, S = Set{eltype(g)}())

Return the set of vertices d-connected to the set of vertices X given set of vertices S in dag g.
"""
function bayesball(g, X, S = Set{eltype(g)}(), veto = no_veto)
    return gensearch(g, X, (pe, ne, v, w) -> !veto(pe, ne, v, w) && (pe == INIT || (v in S && pe == RIGHT && ne == LEFT) || (!(v in S) && !(pe == RIGHT && ne == LEFT))))
end


"""
    bayesball_graph(g, X, S = Set{eltype(g)}(); back=false)

Return an mixed graph `b` containing edges for possible moves of the
Bayes ball. Vertex `x` of `g` is vertex  "`x` forward" at `2x-1` of `b`
if entered forward and "`x` backward" at `2x` if entered backward.
`y` is d-connected to `x` given `S` if and only if there is a 
semi-directed path in `b` from "`x` backward" to "`y` forward"
 or "`y` backward"). `back=true` allows path through `X`
"""
function bayesball_graph(g, X, S = Set{eltype(g)}(); back=false)
    ι(e, i) = e == RIGHT ? 2i-1 : 2i 
    edges = Pair{Int,Int}[]
    CausalInference.gensearch(g, X, (pe, ne, v, w) ->  (pe == INIT || (v in S && pe == RIGHT && ne == LEFT) || (!(v in S) && !(pe == RIGHT && ne == LEFT))) && (back || w ∉ X) && (push!(edges, ι(pe, v)=>ι(ne, w)); true))
    CausalInference.digraph(edges, 2*nv(g))
end

"""
    alt_test_dsep(g, X, Y, S, veto = no_veto)

Check if sets of vertices `X` and `Y` are d-separated in `g` given `S`. 

An alternative to the `test_dsep` function, which uses gensearch under the hood. Might be (a bit) slower. 
"""
function alt_test_dsep(g, X, Y, S, veto = no_veto)
    X, Y, S = toset.((X, Y, S))
    return length(intersect(bayesball(g, X, S, veto), Y)) == 0 
end 

"""
    test_covariate_adjustment(g, X, Y, S)

Check if `S` is a covariate adjustment set relative to `(X, Y)` in graph `g`. 

Based on the sound and complete graphical criterion for covariate adjustment given in https://arxiv.org/abs/1203.3515 using the algorithmic approach proposed in https://arxiv.org/abs/1803.00116. Output is a boolean.
"""
function test_covariate_adjustment(g, X, Y, S)
    X, Y, S = toset.((X, Y, S))
    PCPXY = pcp(g, X, Y)
    length(intersect(S, descendants(g, PCPXY))) != 0 && return false
    return alt_test_dsep(g, X, Y, S, (pe, ne, v, w) -> v in X && w in PCPXY && ne == RIGHT)
end

"""
    alt_test_backdoor(g, X, Y, S)

Check if `S` satisfies the backdoor criterion relative to `(X, Y)` in graph `g`. 

The generalization to sets X and Y differs from, e.g., Pearl (2009). See the Example section (TODO: ref).
"""
function alt_test_backdoor(g, X, Y, S)
    X, Y, S = toset.((X, Y, S))
    length(intersect(S, descendants(g, X))) != 0 && return false
    return alt_test_dsep(g, X, Y, S, no_outgoing(X))
end

"""
    find_dsep(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y), veto = no_veto)

Find a d-separator `Z` with ``I ⊆ Z ⊆ R`` for sets of vertices `X` and `Y` in `g`, else return `false`. 

Based on the algorithmic approach proposed in https://arxiv.org/abs/1803.00116. 
"""
function find_dsep(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y), veto = no_veto)
    X, Y, I, R = toset.((X, Y, I, R))
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

Find an inclusion minimal d-separator `Z` with ``I ⊆ Z ⊆ R`` for sets of vertices `X` and `Y` in `g`, i.e., one for which no subset is a d-separator, else return `false`.

Based on the algorithmic approach proposed in http://auai.org/uai2019/proceedings/papers/222.pdf. 
"""
function find_min_dsep(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y), veto = no_veto)
    X, Y, I, R = toset.((X, Y, I, R))
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

Find a covariate adjustment set `Z` with ``I ⊆ Z ⊆ R`` for sets of vertices `X` and `Y` in `g`, else return `false`. 

Based on the algorithmic approach proposed in https://arxiv.org/abs/1803.00116. 
"""
function find_covariate_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    X, Y, I, R = toset.((X, Y, I, R))
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

Find a backdoor adjustment set `Z` with ``I ⊆ Z ⊆ R`` for sets of vertices `X` and `Y` in `g`, else return `false`. 

The generalization to sets X and Y differs from, e.g., Pearl (2009). See the Example section (TODO: ref).
"""
function find_backdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    (X, Y, I, R) = toset.((X, Y, I, R))
    bdZ = setdiff(find_covariate_adjustment(g, X, Y, I, R), descendants(g, X, no_incoming(X)))
    # we generalize bd to sets by testing (X indep Y given Z) in G with outgoing edges from X removed
    # instead of checking the bd criterion for pairs of vertices (x in X and y in Y) as originally proposed by Pearl which is more restrictive
    if issubset(I, bdZ) && alt_test_dsep(g, X, Y, bdZ, no_outgoing(X))
        return bdZ
    else
        return false
    end
end

"""
    find_min_covariate_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))

Find an inclusion minimal covariate adjustment set `Z` with ``I ⊆ Z ⊆ R`` for sets of vertices `X` and `Y` in `g`, else return `false`. 

Based on the algorithmic approach proposed in https://arxiv.org/abs/1803.00116.  
"""
function find_min_covariate_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    X, Y, I, R = toset.((X, Y, I, R))
    PCPXY = pcp(g, X, Y)    
    return find_min_dsep(g, X, Y, I, setdiff(R, descendants(g, PCPXY)), (pe, ne, v, w) -> v in X && w in PCPXY && ne == RIGHT)
end

"""
    find_min_backdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))

Find an inclusion minimal backdoor adjustment set `Z` with ``I ⊆ Z ⊆ R`` for sets of vertices `X` and `Y` in `g`, else return `false`. 

The generalization to sets X and Y differs from, e.g., Pearl (2009). See the Example section (TODO: ref).
"""
function find_min_backdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    X, Y, I, R = toset.((X, Y, I, R))
    return find_min_dsep(g, X, Y, I, setdiff(R, descendants(g, X)), (pe, ne, v, w) -> v in X && ne == RIGHT)
end

"""
    find_frontdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))

Find a frontdoor adjustment set `Z` with ``I ⊆ Z ⊆ R`` for sets of vertices `X` and `Y` in `g`, else return `false`. 

Based on the algorithm given in  https://arxiv.org/abs/2211.16468. 
"""
function find_frontdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    X, Y, I, R = toset.((X, Y, I, R))
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

Find an inclusion minimal frontdoor adjustment set `Z` with ``I ⊆ Z ⊆ R`` for sets of vertices `X` and `Y` in `g`, else returns `false`. 

Based on the algorithm given in https://arxiv.org/abs/2211.16468. 
"""
function find_min_frontdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    X, Y, I, R = toset.((X, Y, I, R))
    # write explanation for the pass rules
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

struct ConstraintIterator{T<:Integer, S, U<:AbstractSet{T}, F<:Function}
    g::SimpleDiGraph{T}
    X::S
    Y::S
    I::U
    R::U
    find::F
end

function downwards(state, I, R)
    v = first(setdiff(R, I))
    push!(state, (:up, :I, v, I, R))
    push!(state, (:down, :I, v, I, R))
    push!(state, (:up, :R, v, I, R))
    push!(state, (:down, :R, v, I, R))
end

function Base.iterate(A::ConstraintIterator)
    R = copy(A.R)
    I = copy(A.I)
    state = Vector{Tuple{Symbol, Symbol, eltype(A.g), Set{eltype(A.g)}, Set{eltype(A.g)}}}()
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
        issetequal(I, R) && return copy(I), state
        if dir == :down
            downwards(state, I, R)
        end
    end
    return nothing
end

Base.IteratorSize(::ConstraintIterator) = Base.SizeUnknown()

"""
    list_dseps(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))

List all d-separators `Z` with ``I ⊆ Z ⊆ R`` for sets of vertices `X` and `Y` in `g`. 
"""
function list_dseps(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    X, Y, I, R = toset.((X, Y, I, R))
    return ConstraintIterator(g, X, Y, I, R, find_dsep)
end

"""
    list_covariate_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))

List all covariate adjustment sets `Z` with ``I ⊆ Z ⊆ R`` for sets of vertices `X` and `Y` in `g`. 
"""
function list_covariate_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    X, Y, I, R = toset.((X, Y, I, R))
    return ConstraintIterator(g, X, Y, I, R, find_covariate_adjustment)
end

"""
    list_backdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))

List all back-door adjustment sets `Z` with ``I ⊆ Z ⊆ R`` for sets of vertices `X` and `Y` in `g`. 
"""
function list_backdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    X, Y, I, R = toset.((X, Y, I, R))
    return ConstraintIterator(g, X, Y, I, R, find_backdoor_adjustment)
end

# --- Bool-array frontdoor find (Wienöbst et al.; cf. mwien/frontdoor-adjustment) ---
# Used by `list_frontdoor_adjustment` for lower constant factors than the
# gensearch-based `find_frontdoor_adjustment` while preserving the same I/R
# ConstraintIterator enumeration (O(n(n+m)) delay).

function _frontdoor_bdbb!(G, Z, X, v, edgetype, visited, reachable)
    visited[2 * v + edgetype] = true
    reachable[v] = true
    if !(v in Z) && !(v in X)
        for u in outneighbors(G, v)
            if !visited[2 * u]
                _frontdoor_bdbb!(G, Z, X, u, 0, visited, reachable)
            end
        end
    end
    if (edgetype == 0 && (v in Z)) || (edgetype == 1 && !(v in Z))
        for u in inneighbors(G, v)
            if !visited[2 * u + 1]
                _frontdoor_bdbb!(G, Z, X, u, 1, visited, reachable)
            end
        end
    end
end

function _frontdoor_reachbd(G, X, Z)
    visited = falses(2 * nv(G) + 1)
    reachable = falses(nv(G))
    for x in X
        for edgetype in (0, 1)
            if !visited[2 * x + edgetype]
                _frontdoor_bdbb!(G, Z, X, x, edgetype, visited, reachable)
            end
        end
    end
    return Set{eltype(G)}(findall(identity, reachable))
end

function _frontdoor_modbb!(G, X, v, edgetype, visited, continuelater, forbidden)
    visited[2 * v + edgetype] = true
    forbidden[v] = true
    if !(v in X)
        for u in outneighbors(G, v)
            if !visited[2 * u]
                _frontdoor_modbb!(G, X, u, 0, visited, continuelater, forbidden)
            end
        end
    end
    if (edgetype == 0 && v in X) || (edgetype == 1 && !(v in X))
        for u in inneighbors(G, v)
            if forbidden[u] && !visited[2 * u + 1]
                _frontdoor_modbb!(G, X, u, 1, visited, continuelater, forbidden)
            elseif !visited[2 * u + 1]
                continuelater[u] = true
            end
        end
    end
    if continuelater[v] && !visited[2 * v + 1]
        _frontdoor_modbb!(G, X, v, 1, visited, continuelater, forbidden)
    end
end

function _frontdoor_dfs!(G, Z, u, visited)
    visited[u] = true
    for v in outneighbors(G, u)
        if !(v in Z) && !visited[v]
            _frontdoor_dfs!(G, Z, v, visited)
        end
    end
end

function _frontdoor_checkfirst(G, X, Y, Z)
    visited = falses(nv(G))
    for x in X
        if !visited[x]
            _frontdoor_dfs!(G, Z, x, visited)
        end
    end
    for y in Y
        visited[y] && return false
    end
    return true
end

"""
Find a frontdoor adjustment set with ``I ⊆ Z ⊆ R`` using the bool-array
formulation from Wienöbst et al. (2024) / mwien/frontdoor-adjustment.
Same contract as [`find_frontdoor_adjustment`](@ref): return `Z` or `false`.
"""
function _find_frontdoor_adjustment_bool(g, X, Y, I, R)
    Za = setdiff(R, _frontdoor_reachbd(g, X, Set{eltype(g)}()))
    visited = falses(2 * nv(g) + 1)
    continuelater = falses(nv(g))
    save = setdiff(Set(vertices(g)), Za)
    forbidden = [x in save for x in 1:nv(g)]
    for y in Y
        for edgetype in (0, 1)
            _frontdoor_modbb!(g, X, y, edgetype, visited, continuelater, forbidden)
        end
    end
    Zab = setdiff(Za, Set{eltype(g)}(findall(identity, forbidden)))
    if issubset(I, Zab) && _frontdoor_checkfirst(g, X, Y, Zab)
        return Zab
    else
        return false
    end
end

"""
    list_frontdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))

List all front-door adjustment sets `Z` with ``I ⊆ Z ⊆ R`` for sets of vertices `X` and `Y` in `g`.

Uses the same `ConstraintIterator` I/R enumeration as the other `list_*`
methods (``O(n(n + m))`` delay with linear-time find; Wienöbst et al. 2024,
Corollary 1). The find oracle is a bool-array implementation of the paper's
Algorithm 2 (cf. [mwien/frontdoor-adjustment](https://github.com/mwien/frontdoor-adjustment)),
which has lower constant factors than [`find_frontdoor_adjustment`](@ref) under
heavy enumeration while returning the same family of sets.
"""
function list_frontdoor_adjustment(g, X, Y, I = Set{eltype(g)}(), R = setdiff(Set(vertices(g)), X, Y))
    X, Y, I, R = toset.((X, Y, I, R))
    return ConstraintIterator(g, X, Y, I, R, _find_frontdoor_adjustment_bool)
end

# note that this enumeration approach does *not* work for *minimal* dseps and adjustment sets
