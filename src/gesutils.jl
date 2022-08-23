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


function semi_directed_paths(dg,u,v, trajs=[[]], current = u,visited = zeros(Bool, nv(dg)))
    visited[u] = true
    while true
        for n in outneighbors(dg,current)
            

@memoize function semi_directed_paths(s, d, g, dg, dict)
    """Return all paths from i to j in A. Note: a path is a sequence
    (a_1,...,a_n) of non-repeating nodes where either a_i -> a_i+1 or
    a_i - a_i+1 are edges in the PDAG A."""
    all_trajs::Vector{Vector{Any}} = []
    possible_trajs = union(inneighbors(dg, d), neighbors(g, d))
    if s in possible_trajs
        return [[s]]
    elseif length(possible_trajs) == 0
        return [[]]
    end
    for end_node in possible_trajs
        @show end_node
        trajs = semi_directed_paths(s, end_node, g, dg, dict)
        
        for traj in trajs
            @show(traj)
            if !(end_node in traj)
                push!(traj, end_node)
                push!(all_trajs, traj)
            end
        end
    end
    @show all_trajs
    return all_trajs
end



