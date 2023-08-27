using Graphs
using LinkedLists


"""
    vispush!(l::LinkedList, pointers, x, vis)

"""
@inline function vispush!(l::LinkedLists.LinkedList, pointers, x, vis)
    if vis
        pointers[x] = push!(l,x)
    else
        pointers[x] = pushfirst!(l,x)
    end
end

"""
    countmcs(G)

TODO.

"""

function countmcs(G)
    n = nv(G)
    sets = [LinkedLists.LinkedList{Int64}() for i = 1:n+1]
    pointers = Vector{ListNode{Int64}}(undef,n)
    size = ones(Int64, n)
    visited = falses(n) 
    invmcsorder = zeros(Int64, n)

    # init
    for v in vertices(G)
        vispush!(sets[1], pointers, v, visited[v])
    end
    maxcard = 1

    for i = 1:n
        v = first(sets[maxcard])
        size[v] = -size[v] + 1
        invmcsorder[v] = i

        deleteat!(sets[maxcard], pointers[v])

        # update the neighbors
        for w in inneighbors(G, v)
            if size[w] >= 1
                deleteat!(sets[size[w]], pointers[w])
                size[w] += 1
                vispush!(sets[size[w]], pointers, w, visited[w])
            end
        end
        maxcard += 1
        while maxcard >= 1 && isempty(sets[maxcard])
            maxcard -= 1
        end
    end

    return -size, invmcsorder
    
end

"""
    mcs(G, K)

Perform a Maximum Cardinality Search on graph G. The elements of clique K are of prioritized and chosen first.If K is empty a normal MCS is performed. Return the visit order.

"""
function mcs(G, K)
    n = nv(G)
    copy_K = copy(K)

    sets = [LinkedLists.LinkedList{Int64}() for i = 1:n+1]
    pointers = Vector{ListNode{Int64}}(undef,n)
    size = ones(Int64, n)
    visited = falses(n) 
    invmcsorder = zeros(Int64, n)

    # init
    for v in vertices(G)
        vispush!(sets[1], pointers, v, visited[v])
    end
    maxcard = 1

    for i = 1:n
        # first, the vertices in K are chosen
        # they are always in the set of maximum cardinality vertices
        if !isempty(copy_K)
            v = popfirst!(copy_K)
        # afterwards, the algorithm chooses any vertex from maxcard
        else
            v = first(sets[maxcard])
        end
        invmcsorder[v] = i
        size[v] = -1

        deleteat!(sets[maxcard], pointers[v])

        # update the neighbors
        for w in inneighbors(G, v)
            if size[w] >= 1
                deleteat!(sets[size[w]], pointers[w])
                size[w] += 1
                vispush!(sets[size[w]], pointers, w, visited[w])
            end
        end
        maxcard += 1
        while maxcard >= 1 && isempty(sets[maxcard])
            maxcard -= 1
        end
    end

    return invmcsorder
end

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
    invmcsorder = mcs(c, K)
    
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

    return cpdag(d)
end
