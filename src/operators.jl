using CliqueTrees: mcs
using Graphs


"""
    count_mcs(G)

Perform a Maximum Cardinality Search on graph G. 
"""
function count_mcs(G)
    invmcsorder, size = mcs(G)
    
    for i in vertices(G)
        invmcsorder[i] = 1 + nv(G) - invmcsorder[i]
        size[i] = -size[i]
    end

    size, invmcsorder
end

"""
    operator_mcs(G, K)

Perform a Maximum Cardinality Search on graph G. The elements of clique K are of prioritized and chosen first.
"""
function operator_mcs(G, K)
    invmcsorder = first(mcs(G, reverse(K)))
    
    for i in vertices(G)
        invmcsorder[i] = 1 + nv(G) - invmcsorder[i]
    end

    invmcsorder
end

# EXPORTS:
"""
    next_CPDAG(g, op, x, y, S)

`next_CDPAG` applies an operator and computes the resulting CPDAG in linear-time 
"""
function next_CPDAG(g, op, x, y, S)
    n = nv(g)
    c = copy(g)
    c.ne = 0
    for i = 1:n
        filter!(j->has_edge(g, j, i), c.fadjlist[i])
        filter!(j->has_edge(g, i, j), c.badjlist[i])
        c.ne += length(c.fadjlist[i])
    end
    K = Vector{Int64}()
    if op == :up
        append!(K, S)
        append!(K, adj_neighbors(g, x, y))
        push!(K, y)
    end
    if op == :down
        append!(K, setdiff(adj_neighbors(g, x, y), S))
        if isundirected(g, x, y)
            push!(K, x)
        end
        push!(K, y)
    end
    invmcsorder = operator_mcs(c, K)
    
    # maybe this is not necessary
    comp = Vector{Int64}(undef, n)
    components = connected_components(c)
    for i = 1:length(components)
        for x in components[i]
            comp[x] = i
        end
    end
    
    d = copy(g)
    d.ne = 0
    for i = 1:n
        filter!(j->comp[i] != comp[j] || invmcsorder[j] > invmcsorder[i], d.fadjlist[i])
        filter!(j->comp[i] != comp[j] || invmcsorder[j] < invmcsorder[i], d.badjlist[i])
        d.ne += length(d.fadjlist[i])
    end

    # now apply operation
    if op == :up
        add_edge!(d, x, y)
    end
    if op == :down
        rem_edge!(d, x, y)
        rem_edge!(d, y, x)
    end

    return alt_cpdag(d)
end

## use for tests
#g = SimpleDiGraph(Edge.([(1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2), (2, 4), (4, 2), (3, 4), (4, 3)]))
#vmap = [1, 2, 3, 4]
#h = SimpleDiGraph(Edge.([(1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2), (2, 4), (4, 2), (3, 4), (4, 3), (1, 5), (4, 5), (4, 6), (6, 4), (4, 7), (7, 4), (6, 7), (7, 6), (2, 6), (6, 2), (3, 6), (6, 3)]))
#x = 1
#y = 4

# assumes that g is a chordal graph
struct CliqueIterator{T<:Integer}
    g::SimpleDiGraph{T}
    vmap::Vector{T}
end

function Base.iterate(C::CliqueIterator)
    g = C.g
    n = nv(g)
    _, invmcsorder = count_mcs(g)
    P = [Vector{Int64}() for _ = 1:n]
    for i = 1:n
        for j in neighbors(g, i)
            invmcsorder[j] < invmcsorder[i] && push!(P[i], j)
        end
    end
    # TODO: replace int by bitvector or something
    state = (1, 0, P)
    return Vector{Int64}(), state
end

function Base.iterate(C::CliqueIterator, state)
    v = state[1]
    if v > nv(C.g)
        return nothing
    end
    idx = state[2]
    P = state[3][v]
    if idx >= 2^length(P)
        return Base.iterate(C, (v+1, 0, state[3])) 
    else
        clique = P[digits(Bool, idx, base=2, pad=length(P))]
        push!(clique, v)
    end
    return map(v -> C.vmap[v], clique), (v, idx+1, state[3])
end

"""
    precompute_semidirected(g, y)

Computes for vertex y all vertices reachable via semidirected path from any undirected neighbor and y 
itself with all vertices in this same set blocked.
"""
function precompute_semidirected(g, y)
    n = nv(g)
    blocked = falses(n)
    nu = neighbors_undirected(g, y)
    push!(nu, y)
    for v in nu
        blocked[v] = true
    end
    semidirected = Vector{BitVector}()
    for z in nu
        vis = falses(n)
        q = Vector{Int64}()
        push!(q, z)
        vis[z] = true
        while !isempty(q)
            w = popfirst!(q)
            for v in outneighbors(g, w)
                vis[v] && continue
                blocked[v] && continue
                push!(q, v)
                vis[v] = true
            end
        end
        push!(semidirected, vis)
    end 
    return semidirected
end

# needs semidirected precomputed with precompute_semidirected(g, y) 
# for efficiency have y in outer loop!
"""
    InsertIterator{T<:Integer}

Lists all insert operators for pair x,y (needs semidirected from prev point)
"""
struct InsertIterator{T<:Integer}
    g::SimpleDiGraph{T}
    x::T
    y::T
    semidirected::Vector{BitVector}
end

Base.IteratorSize(::InsertIterator) = Base.SizeUnknown()

function Base.iterate(O::InsertIterator)
    g = O.g
    x = O.x
    y = O.y
    semidirected = O.semidirected
    if x == y || has_edge(g, x, y) || last(semidirected)[x]
        return nothing
    end
    nu = neighbors_undirected(g, y)
    if length(nu) == 0
        (empty, it) = Iterators.peel([Vector{Int64}()])
        return empty, (it, Vector{Int64}())
    end
    push!(nu, y)
    an = adj_neighbors(g, x, y)
    musttake = Set(an)
    for i = 1:length(nu)
        if semidirected[i][x]
            push!(musttake, nu[i])
        end
    end
    !isclique(g, musttake) && return nothing
    cantake = Vector{Int64}()
    for v in nu
        v == y && continue
        v in musttake && continue
        fullyconnected = true
        for w in musttake
            if !isadjacent(g, v, w)
                fullyconnected = false
                break
            end
        end
        fullyconnected && !isadjacent(g, x, v) && push!(cantake, v)
    end
    cliqueit = CliqueIterator(induced_subgraph(g, cantake)...)
    state = (cliqueit, setdiff(musttake, an))
    return Base.iterate(O, state)
end

function Base.iterate(O::InsertIterator, state)
    cliqueit = state[1]
    musttake = state[2]
    res = Iterators.peel(cliqueit)
    res === nothing && return nothing 
    clique = res[1]
    cliqueit = res[2]
    append!(clique, musttake)
    return clique, (cliqueit, musttake)
end

"""
    DeleteIterator{T<:Integer}

Lists all delete operators for pair x,y.
"""
struct DeleteIterator{T<:Integer}
    g::SimpleDiGraph{T}
    x::T
    y::T
end

Base.IteratorSize(::DeleteIterator) = Base.SizeUnknown()

function Base.iterate(O::DeleteIterator)
    g = O.g
    x = O.x
    y = O.y
    if !has_edge(g, x, y)
        return nothing
    end
    an = adj_neighbors(g, x, y)
    if length(an) == 0
        (empty, it) = Iterators.peel([Vector{Int64}()])
        return empty, (it, an)
    end
    if length(an) == 1
        (empty, it) = Iterators.peel([Vector{Int64}(), Vector{Int64}(an)])
        return setdiff(an, empty), (it, an)
    end
    cliqueit = CliqueIterator(induced_subgraph(g, an)...)
    state = (cliqueit, an)
    return Base.iterate(O, state)
end

function Base.iterate(O::DeleteIterator, state)
    cliqueit = state[1]
    an = state[2]
    res = Iterators.peel(cliqueit)
    res === nothing && return nothing
    clique = res[1]
    cliqueit = res[2]
    return setdiff(an, clique), (cliqueit, an)
end

## TODO: refactor ##
function countcliques(g)
    n = nv(g)
    preceding, _ = count_mcs(g)
    # maybe use BigInt at some point
    cnt = 1 # don't forget "empty" clique
    for i = 1:n
        cnt += 2^preceding[i]
    end
    return cnt
end

function sampleclique(g, r)
    n = nv(g)
    preceding, invmcsorder = count_mcs(g)
    cnt = 1 # don't forget "empty" clique
    r <= cnt && return Vector{Int64}()
    for i = 1:n
        cnt += 2^preceding[i]
        if r <= cnt
            p = Vector{Int64}()
            for j in neighbors(g, i)
                invmcsorder[j] < invmcsorder[i] && push!(p, j) 
            end
             ret = randsubseq(p, 0.5)
        #    subset = rand(0:2^preceding[i]-1)
        #    ret = Vector{Int64}()
         #   for d = 0:length(p)
         #       if ((1 << (d-1)) & subset) > 0
         #           push!(ret, p[d])
         #       end
         #   end
            push!(ret, i)
            return ret
        end
    end
end

function exactup(g, κ)
    n = nv(g)
    ttcnt = 0
    cnts = Vector{Int64}()
    cntids = Vector{Tuple{Int64, Int64}}()
    need = Vector{Set{Int64}}()
    pot = Vector{Vector{Int64}}()
    for y in vertices(g)
        length(neighbors_adjacent(g, y)) == κ && continue
        blocked = falses(n)
        nu = neighbors_undirected(g, y)
        push!(nu, y)
        for v in nu
            blocked[v] = true
        end
        visfrom = Vector{BitVector}()
        for z in nu
            # => find all vertices reachable by semidirected path (with other neighbors of y and y itself blocked) from z and save somewhere
            vis = falses(n)
            q = Vector{Int64}()
            push!(q, z)
            vis[z] = true
            while !isempty(q)
                w = popfirst!(q)
                for v in outneighbors(g, w)
                    vis[v] && continue
                    blocked[v] && continue
                    push!(q, v)
                    vis[v] = true
                end
            end
            push!(visfrom, vis)
        end
        for x in vertices(g)
            x == y && continue
            isadjacent(g, x, y) && continue
            last(visfrom)[x] && continue
            length(neighbors_adjacent(g, x)) == κ && continue
            needtotake = Set(adj_neighbors(g, x, y))
            # => get all other undirected neighbors we have to take to close semidirected paths
            for i = 1:length(nu)
                if visfrom[i][x]
                    push!(needtotake, nu[i])
                end
            end
            !isclique(g, needtotake) && continue
            # => keep all other undirected neighbors which are connected to all the must have ones
            potential = Vector{Int64}()
            for v in nu
                v == y && continue
                v in needtotake && continue
                fullyconnected = true
                for w in needtotake
                    if !isadjacent(g, v, w)
                        fullyconnected = false
                        break
                    end
                end
                fullyconnected && push!(potential, v)
            end
            # => count number of cliques in chordal graph
            sg, _ = induced_subgraph(g, potential)
            cnt = countcliques(sg)
            ttcnt += cnt
            push!(cnts, ttcnt)
            push!(cntids, (x,y))
            push!(need, needtotake)
            push!(pot, potential)
            # println(x, " ", y, " ", ttcnt) 
        end
    end
    # store counts for each pair x,y and sample one of them
    # then sample a clique randomly in chordal graph -> by same procedure find# highest index node and then take random subset of prev visited neighbors
    ttcnt == 0 && return ttcnt, (0, 0, Int[]) 
    r = rand(1:ttcnt)
    for i = 1:length(cnts)
        if r <= cnts[i]
            (x, y) = cntids[i]
            cnt = cnts[i]
            i > 1 && (cnt -= cnts[i-1])
            r2 = rand(1:cnt)
            sg, vmap = induced_subgraph(g, pot[i])
            cl = map(v -> vmap[v], sampleclique(sg, r2))
            append!(cl, need[i])
            return ttcnt, (x, y, setdiff(cl, adj_neighbors(g, x, y)))
        end
    end
end

function exactdown(g)
    ttcnt = 0
    cnts = Vector{Int64}()
    cntids = Vector{Tuple{Int64, Int64}}()
    for x in vertices(g)
        for y in [neighbors_undirected(g, x); children(g, x)]
            NAyx = adj_neighbors(g, x, y)
            if length(NAyx) == 0
                ttcnt += 1
            elseif length(NAyx) == 1
                ttcnt += 2
            else
                sg, _ = induced_subgraph(g, NAyx)
                cnt = countcliques(sg)
                ttcnt += cnt
            end
            push!(cnts, ttcnt)
            push!(cntids, (x,y))
        end
    end 
    ttcnt == 0 && return 0, (0, 0, Int[])
    # sample clique
    r = rand(1:ttcnt)
    for i = 1:length(cnts)
        if r <= cnts[i]
            (x, y) = cntids[i]
            cnt = cnts[i]
            i > 1 && (cnt -= cnts[i-1])
            r2 = rand(1:cnt)
            NAyx = adj_neighbors(g, x, y)
            sg, vmap = induced_subgraph(g, NAyx)
            cl = map(v -> vmap[v], sampleclique(sg, r2))
            return ttcnt, (x, y, setdiff(NAyx, cl)) 
        end
    end
end

"""
    count_moves_uniform(g, κ=nv(g) - 1) = s1, s2, (x1, y1, T1), (x2, y2, H2)

Counts and samples operator in polynomial-time by avoiding full enumeration (only works for uniform score.) 
Count the number `s1` of Insert and `s2` of Delete operators for CPDAG `g` with 
degree bound `κ` and return a uniformly selected `Insert(x1, y1, T1)`` and a uniform selected 
`Delete(x2, y2, H2)` operator. 
"""
function count_moves_uniform(g, κ=nv(g) - 1)
    s1, (x1, y1, T1) = exactup(g, κ)
    s2, (x2, y2, H2) = exactdown(g)
    return s1, s2, (x1, y1, T1), (x2, y2, H2)
end
