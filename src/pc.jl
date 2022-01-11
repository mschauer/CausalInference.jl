using Graphs, Tables, Statistics, Distributions
using TikzGraphs

"""
    insorted(a, x)

Check if `x` is in the sorted collection `a`
"""
insorted(a, x) = !isempty(searchsorted(a, x))

"""
    disjoint_sorted(u, v)

Check if the intersection of sorted collections is empty. The intersection
of empty collectios is empty.
"""
function disjoint_sorted(u, v)
    xa = iterate(u)
    xa === nothing && return true
    yb = iterate(v)
    yb === nothing && return true

    x, a = xa
    y, b = yb

    while true
        x == y && return false
        if x > y
            yb = iterate(v, b)
            yb === nothing && return true
            y, b = yb
        else
            xa = iterate(u, a)
            xa === nothing && return true
            x, a = xa
        end
    end
end
#=
Let e=(v, w).
Note that each unshielded triple is of type `∨` or `╎` or `∧` of shape
     v           v       v
    / \          |        \   z
   w   \         w         \ /
        z        |          w
                 z
where higher nodes have smaller vertex number.
=# 
"""
    orientable_unshielded(g, S)

Find the orientable unshielded triples in the skeleton. Triples are connected vertices v-w-z where
z is not a neighbour of v. Uses that `edges` iterates in lexicographical order.
""" 
function orientable_unshielded(g, S)
    Z = Tuple{Int64,Int64,Int64}[]
    for e in edges(g)
        v, w = Tuple(e)
        @assert(v < w)
        for z in neighbors(g, w) # case `∨` or `╎`
            z <= v && continue   # longer arm of `∨` is visited first
            insorted(neighbors(g, z), v) && continue
            w in S[Edge(v,z)] || push!(Z, (v, w, z))
        end
        for z in neighbors(g, v) # case `∧` 
            (z <= w) && continue # shorter arm is visited first
            insorted(neighbors(g, z), w) && continue
            v in S[Edge(minmax(z, w)...)] || push!(Z, (z, v, w))
        end
    end
    Z
end

"""
    unshielded(g)

Find the unshielded triples in the cyclefree skeleton. Triples are connected vertices v-w-z where
z is not a neighbour of v. Uses that `edges` iterates in lexicographical order.
""" 
function unshielded(g)
    Z = Tuple{Int64,Int64,Int64}[]
    for e in edges(g)
        v, w = Tuple(e)
        @assert(v < w)
        for z in neighbors(g, w) # case `∨` or `╎`
            z <= v && continue   # longer arm of `∨` is visited first
            insorted(neighbors(g, z), v) && continue
            push!(Z, (v, w, z))
        end
        for z in neighbors(g, v) # case `∧` 
            (z <= w) && continue # shorter arm is visited first
            insorted(neighbors(g, z), w) && continue
            push!(Z, (z, v, w))
        end
    end
    Z
end


isadjacent(dg, v, w) = has_edge(dg, v, w) || has_edge(dg, w, v)
has_both(dg, v, w) = has_edge(dg, v, w) && has_edge(dg, w, v)

remove!(dg::DiGraph, e::Pair) = rem_edge!(dg, Edge(e))
remove!(dg::Graph, e::Tuple) = rem_edge!(dg, Edge(e))

"""
    vskel(g)

Skeleton and `v`-structures. (Currently from the first step of the pc-Alg.)
"""
vskel(g) = _vskel(nv(g), dseporacle, g)

function _vskel(n::V, I, par...) where {V}
    # Step 1
    g, S = skeleton(n, I, par...)

    # Step 2: Apply Rule 0 once
    Z = orientable_unshielded(g, S)
    dg = DiGraph(g) # use g to keep track of unoriented edges

    for (u, v, w) in Z
        if has_edge(g, (u, v))
            remove!(dg, v => u)
            remove!(g, (v, u))
        end
        if has_edge(g, (v, w))
            remove!(dg, v => w)
            remove!(g, (v, w))
        end
    end
    dg
end


"""
    pcalg(n::V, I, par...)
    pcalg(g, I, par...)

Perform the PC algorithm for a set of 1:n variables using the tests

    I(u, v, [s1, ..., sn], par...)

Returns the CPDAG as DiGraph.   
"""
function pcalg(n, I, par...; kwargs...) 
    g = complete_graph(n)
    VERBOSE = false
    # Step 1
    g, S = skeleton(g, I, par...; kwargs...)
    dg = DiGraph(g) # use g to keep track of unoriented edges
    g, dg = orient_unshielded(g, dg, S) 
    apply_pc_rules(g, dg; kwargs...)
end


"""
    orient_unshielded(g, dg, S)

Orient unshielded triples using the seperating sets.
`g` is an undirected graph containing edges of unknown direction,
`dg` is an directed graph containing edges of known direction and 
both `v=>w` and `w=>v `if the direction of `Edge(v,w)`` is unknown.
`S` are the separating sets of edges.

Returns `g, dg`.
"""
function orient_unshielded(g, dg, S)
    VERBOSE = false
    V = eltype(g)
    n = nv(g)

    # Step 2: Apply Rule 0 once
    Z = orientable_unshielded(g, S)


    for (u, v, w) in Z
        if has_edge(g, (u, v))
            remove!(dg, v => u)
            remove!(g, (v, u))
        end
        if has_edge(g, (v, w))
            remove!(dg, v => w)
            remove!(g, (v, w))
        end
    end
    g, dg
end

"""
    apply_pc_rules(g, dg)


`g` is an undirected graph containing edges of unknown direction,
`dg` is an directed graph containing edges of known direction and 
both `v=>w` and `w=>v `if the direction of `Edge(v,w)`` is unknown.

Returns the CPDAG as DiGraph. 
"""  
function apply_pc_rules(g, dg; VERBOSE = false)
    # Step 3: Apply Rule 1-3 consecutively
    removed = Tuple{Int64,Int64}[]
    while true
        for e in edges(g)
            for e_ in (e, reverse(e))
                v, w = Tuple(e_)
                # Rule 1: Orient v-w into v->w whenever there is u->v
                # such that u and w are not adjacent
                for u in inneighbors(dg, v)
                    has_edge(dg, v => u) && continue # not directed
                    isadjacent(dg, u, w) && continue
                    VERBOSE && println("rule 1: ", v => w) 
                    remove!(dg, w => v)
                    push!(removed, (w, v))
                    @goto ende
                end
            

                # Rule 2: Orient v-w into v->w whenever there is a chain
                # v->k->w.
                outs = Int[]
                for k in outneighbors(dg, v)
                    !has_edge(dg, k => v) && push!(outs, k)
                end
                ins = Int[]
                for k in inneighbors(dg, w)
                    !has_edge(dg, w => k) && push!(ins, k)
                end
                
                if !disjoint_sorted(ins, outs)
                    VERBOSE && println("rule 2: ", v => w) 
                    remove!(dg, w => v)
                    push!(removed, (w, v))
                    @goto ende
                end

                # Rule 3: Orient v-w into v->w whenever there are two chains
                # v-k->w and v-l->w such that k and l are nonadjacent
                fulls = [] # Find nodes k where v-k
                for k in outneighbors(dg, v)
                    has_edge(dg, k => v) && push!(fulls, k)
                end
                for (k, l) in combinations(fulls, 2) # FIXME: 
                    isadjacent(dg, k, l) && continue
                    
                    # Skip if not k->w or if not l->w
                    if has_edge(dg, w => k) || !has_edge(dg, k => w)
                        continue
                    end
                    if has_edge(dg, w => l) || !has_edge(dg, l => w)
                        continue
                    end
                    VERBOSE && println("rule 3: ", v => w) 
                    remove!(dg, w => v)
                    push!(removed, (w, v))
                    @goto ende
                end
            end
        end
        
        @label ende
        for e in removed
            remove!(g, e)
        end
        isempty(removed) && break 
        empty!(removed)
    end
    dg
end

"""
    pcalg(t, p::Float64, test::typeof(gausscitest); kwargs...)

run PC algorithm for tabular input data t using a p-value p to test for 
conditional independeces using Fisher's z-transformation.
"""
function pcalg(t, p::Float64, test::typeof(gausscitest); kwargs...)
    @assert Tables.istable(t)

    c = Tables.columns(t)
    sch = Tables.schema(t)
    n = length(sch.names)
  
    X = reduce(hcat, map(c->Tables.getcolumn(Tables.columns(t), c), Tables.columnnames(t)))
    N = size(X,1)
    C = Statistics.cor(X)
    return pcalg(n, gausscitest, (C,N), quantile(Normal(), 1-p/2); kwargs...)
end


"""
    pcalg(t::T, p::Float64; cmitest::typeof(cmitest); kwargs...) where{T}

run PC algorithm for tabular input data t using a p-value p to detect 
conditional independeces using a conditional mutual information permutation test.
"""
function pcalg(t, p::Float64, test::typeof(cmitest); kwargs...)
    @assert Tables.istable(t)
    @assert all(t->t==Float64, Tables.schema(t).types)

    c = Tables.columns(t)
    sch = Tables.schema(t)
    n = length(sch.names)

    return pcalg(n, cmitest, c, p; kwargs...)
end


"""
    plot_pc_graph(g, df)

plot the output of the PC algorithm.
"""
function plot_pc_graph(g, node_labels::Array=[])
    plot_g = DiGraph(nv(g))

    if length(node_labels) != nv(g) 
        node_labels = map(string, 1:nv(g))
    end
        
    node_style = "draw, rounded corners, fill=blue!10"
    options = "scale=2"
    
    styles_dict = Dict()
        
    for e in edges(g)
        if has_edge(g, e.dst, e.src)
            if e.src < e.dst # if both, plot once
                add_edge!(plot_g, e.src, e.dst)
                push!(styles_dict, (e.src, e.dst)=>"--")
            end
        else # this is the only edge
            add_edge!(plot_g, e.src, e.dst)
            push!(styles_dict, (e.src, e.dst)=>"->")
        end
    end
    
    TikzGraphs.plot(plot_g, node_labels, edge_styles=styles_dict,
                    node_style=node_style, options=options)
end
