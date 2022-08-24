function na(y, x, A)
    """Return all neighbors of y which are adjacent to x in A."""
    neighbors(g, X) & adj(x, A)
end

function neighbors(g, s)
    """The neighbors of i in A, i.e. all nodes connected to i by an
    undirected edge."""
    neighbors(g, s)
end

function adj(g, dg, s)
    """The adjacent nodes of i in A, i.e. all nodes connected by a
    directed or undirected edge."""
    union!(neighbors(g, s), neighbors(dg, s))
end

function pa(dg, s)
    """The parents of i in A."""
    inneighbors(dg, s)
end

function ch(dg, s)
    """The children of i in A."""
    outneighbors(dg, s)
end

function is_clique(g, s)
    """Check if the subgraph of A induced by nodes S is a clique."""
    s in maximal_cliques(g)
end

function is_dag(g, dg)
    """Checks wether the given adjacency matrix corresponds to a DAG."""
    """maybe not needed"""
    ne(g) == 0 && ne(dg) > 0
end


function topological_ordering(g)
    """Fill up later"""
end


function convert_undirected_to_directed(g, dg)
    dg_copy = deepcopy(dg)
    for s in node(g)
        for d in node(g)
            add_edge!(dg_copy,s,d) #forward
            add_edge!(dg_copy,d,s) #backward
        end
    end
end

function seperates(S, A, B, g)
    """Returns true if the set S separates A from B in G, i.e. if all
    paths in G from nodes in A to nodes in B contain a node in
    S. Exception is raised if S,A and B are not pairwise disjoint."""
    if length(intersect!(a,b)) || length(intersect!(a,s)) || length(intersect!(b,s))
        @error("The sets are not pairwise disjoint")
    end
    for a in A
        for b in B
            for path in collect!(all_simple_paths(a,b,g))
                if length(intersect!(Set(path),S)) == 0
                    return false
                end
            end
        end
    end
    return true
end

function chain_component(i, g)
    """Return all nodes in the connected component of node i after
    dropping all directed edges in G."""
    A = g
    visited = Set()
    to_visit = {i}
    while length(to_visit) > 0
        for j in to_visit
            push!(visited, j)
            to_visit = setdiff(union!(to_visit, neighbors(A, j)) , visited)
        end
    end
    return visited
end
 
function induced_subgraph(S,G)
    """Remove all edges which are not between nodes in S."""
    for s in nodes(S)
        for d in nodes(S)
            if not has_edges(S, s, d)
                rem_edge!(G, s, d)
            end
        end
    end
end

function make_double_edges(S,G)
    """Remove all edges which are not between nodes in S."""
    for s in nodes(S)
        for d in nodes(S)
            if has_edges(S, s, d)
                add_edge!(G, s, d)
                add_edge!(G, s, d)
            end
        end
    end
end

function vstructures(dg)
    """
    Return the v-structures of a DAG or PDAG, given its adjacency matrix.
    """
    colliders = Dict()
    for node in nodes(dg)
        all_parents = inneighbors(node)
        colliders[node] = collect(powerset(all_parents,2,2))
    end
    for node in nodes(dg)
        for collider in colliders[node]
            if has_edge(dg, collider[1], collider[2]) || has_edge(g, collider[1], collider[2]) 
                pop!(collider, colliders[node])
            end
        end
    end
end

function skeleton(g,dg)
    """Return the skeleton of a given graph."""
    g_copy = deepcopy(g)
    for s in nodes(dg)
        for d in nodes(dg)
            if has_edge(dg, s, d) || has_edge(dg,d,s)
                add_edge!(g_copy, s,d)
            end
        end
    end
    g_copy
end

function is_consistent_extension(dg_final, dg, g)
    """Returns True if the DAG G is a consistent extension of the PDAG
    P. Will raise a ValueError exception if the graph G is not a DAG
    (i.e. cycles or undirected edges)."""
    #figure out later how to check for undirected edges
    if simplecycles_hawick_james(dg_final) 
        throw("dg is not a dag")
    end

    same_vstructures = vstructures(dg) == vstructures(dg_final)
    same_skeleton = (skeleton(dg, g) == skeleton(dg_final,g))
    same_orientation = issubset(collect(edges(dg)), collect(edges(dg_final))) 
    return same_orientation && same_skeleton && same_vstructures
end

function pdag_to_cpdag(pdg)
    """
    Transform a PDAG into its corresponding CPDAG. Returns a ValueError
    exception if the given PDAG does not admit a consistent extension.
    """
    dag = pdag_to_dag(pdag)
    return dag_to_cpdag(dag)
end

function label_edges(ordered)
    """Given a DAG with edges labelled according to a total ordering,
    label each edge as being compelled or reverisble."""
    if !(is_dag(ordered))
        Throw("The given Graph is not a DAG")
    end 
end

"""
Find a consistent extension of the given PDAG. Return a ValueError
exception if the PDAG does not admit a consistent extension.
Parameters
""" 
function pdag_to_dag(dg_copy, g_copy)   
    G = deepcopy(dg_copy)
    g = deepcopy(g_copy)
    dg= deepcopy(dg_copy)
    indexes = collect(1:nv(dg))
    while nv(g) > 0
        @show nv(g)
        found = false
        i = 1
        while (!found) && i < nv(dg)
            #condition 1
            sink = length(outneighbors(dg, i)) == 0
            #condition 2
            n_i = neighbors(g, i)
            @show n_i
            adj_i = adj(g, dg, i)
            @show adj_i
            adj_neighbors = all([(setdiff(adj_i, y) <= adj(g,dg,y)) for y in n_i])
            found = sink && adj_neighbors
            @show found
            if found
                real_i = indexes[i]
                @show real_i
                real_neighbors = [indexes[j] for j in n_i]
                for j in real_neighbors
                    add_edge!(G, j, real_i)
                end
                rem_vertices!(g,[i])
                rem_vertices!(dg,[i])
                filter!(e->e≠real_i,indexes)
            else
                i+=1
                @show i
            end
        end
        if !found
            throw("PDAG does not admit consistent extension")
        end
    end
    G
end


function cpdag1(skel)
    g = copy(skel)
    ys = topological_sort_by_dfs(g)
    yp = sortperm(ys)
    n = nv(g)
    m = ne(g)
    label = Dict(zip(map(Pair, edges(g)), repeated(unknown)))

    i = 0 # src index
    j = 0 # dst index
    xs = Int[]
    while true
        while i == 0
            j = j + 1
            if j > n
                @goto ende 
            end
            xs = sort(inneighbors(g, ys[j]), by=x->yp[x])
            i = length(xs)
        end

        x, y = xs[i], ys[j]

        if label[x => y] == unknown
            for w in inneighbors(g, x)
                label[w => x] == compelled || continue
                if !has_edge(g, w => y)
                    for z in inneighbors(g, y)
                        label[z => y] = compelled
                    end
                    @goto next
                else
                    label[w => y] = compelled
                end
            end

            lbl = reversible
            for z in inneighbors(g, y)
                if x == z || has_edge(g, z => x) 
                    continue
                end
                lbl = compelled
            end
            for z in inneighbors(g, y)
                if label[z => y] == unknown
                    label[z => y] = lbl
                end
            end
        end
        @label next

        i = i - 1
    end
    @label ende
    for e in collect(edges(g))
        x, y = Tuple(e)
        if label[x=>y] == reversible
            add_edge!(g, y => x)
        end
    end
    g
end