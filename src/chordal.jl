import LinkedLists

"""
    ischordal(g)

Return true if the given graph is chordal
"""
function ischordal(G)
    mcsorder, invmcsorder, _ = mcs(G, Set())
    
    n = length(mcsorder)
    
    f = zeros(Int, n)
    index = zeros(Int, n)
    for i=n:-1:1
        w = mcsorder[i]
        f[w] = w
        index[w] = i
        for v in neighbors(G, w)
            if invmcsorder[v] > i
                index[v] = i
                if f[v] == v
                    f[v] = w
                end
            end
        end
        for v in neighbors(G, w)
            if invmcsorder[v] > i
                if index[f[v]] > i
                    return false
                end
            end
        end
    end
    return true
end

function cliquetreefrommcs(G, mcsorder, invmcsorder)
    n = nv(G)
    # data structures for the algorithm
    K = Vector{Set}()
    push!(K, Set())
    s = 1
    edgelist = Set{Edge}()
    visited = falses(n)
    clique = zeros(Int, n)

    for i = 1:n
        x = mcsorder[i]
        S = Set{Int}()
        for w in inneighbors(G, x)
            if visited[w]
                push!(S, w)
            end
        end
        
        # if necessary create new maximal clique
        if K[s] != S
            s += 1
            push!(K, S)
            k, _ = findmax(map(x -> invmcsorder[x], collect(S)))
            p = clique[mcsorder[k]]
            push!(edgelist, Edge(p, s))
        end
        
        union!(K[s], x)
        clique[x] = s 
        visited[x] = true;
    end

    T = SimpleGraphFromIterator(edgelist)
    # ensure graph is not empty
    nv(T) == 0 && add_vertices!(T,1)
    return K, T
end

@inline function vispush!(l::LinkedList, pointers, x, vis)
    if vis
        pointers[x] = push!(l,x)
    else
        pointers[x] = pushfirst!(l,x)
    end
end

# TODO: separate mcs and mcs plus cgk??
# Returns the visit order of the vertices, its inverse and the subgraphs C_G(K) (see Def. 1 in [1,2]). If K is empty a normal MCS is performed.
function mcs(G, K)
    n = nv(G)
    copy_K = copy(K)
    
    # data structures for MCS
    sets = [LinkedList{Int}() for _ = 1:n+1]
    pointers = Vector(undef,n)
    size = Vector{Int}(undef, n)
    visited = falses(n)
    
    # output data structures
    mcsorder = Vector{Int}(undef, n)
    invmcsorder = Vector{Int}(undef, n)
    subgraphs = Array[]

    # init
    visited[collect(copy_K)] .= true
    for v in vertices(G)
        size[v] = 1
        vispush!(sets[1], pointers, v, visited[v])
    end
    maxcard = 1

    for i = 1:n
        # first, the vertices in K are chosen
        # they are always in the set of maximum cardinality vertices
        if !isempty(copy_K)
            v = pop!(copy_K)
        # afterwards, the algorithm chooses any vertex from maxcard
        else
            v = first(sets[maxcard])
        end
        # v is the ith vertex in the mcsorder
        mcsorder[i] = v
        invmcsorder[v] = i
        size[v] = -1

        # immediately append possible subproblems to the output
        if !visited[v]
            vertexset = Vector{Int}()
            for x in sets[maxcard]
                visited[x] && break
                visited[x] = true
                push!(vertexset, x)
            end
            sg = induced_subgraph(G, vertexset)
            subgraphs = vcat(subgraphs, (map(x -> sg[2][x], connected_components(sg[1]))))
        end

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

    return mcsorder, invmcsorder, subgraphs
end

"""
    cliquetree(G)

Computes a clique tree of a graph G. A vector K of maximal cliques and a tree T on 1,2,...,|K| is returned.

"""
function cliquetree(G)
    mcsorder, invmcsorder, _ = mcs(G, Set())
    K, T = cliquetreefrommcs(G, mcsorder, invmcsorder)
    return K, T
end
