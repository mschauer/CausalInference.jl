using Graphs, Tables, Statistics, Distributions

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
    is_directed(typeof(g)) && throw(ArgumentError("Argument is directed."))
    Z = Tuple{Int64, Int64, Int64}[]
    for e in edges(g)
        v, w = Tuple(e)
        @assert(v<w)
        for z in neighbors(g, w) # case `∨` or `╎`
            z <= v && continue   # longer arm of `∨` is visited first
            has_edge(g, v,  z) && continue
            w in S[Edge(v, z)] || push!(Z, (v, w, z))
        end
        for z in neighbors(g, v) # case `∧` 
            (z <= w) && continue # shorter arm is visited first
            has_edge(g, w,  z) && continue
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
    Z = Tuple{Int64, Int64, Int64}[]
    for e in edges(g)
        v, w = Tuple(e)
        @assert(v<w)
        for z in neighbors(g, w) # case `∨` or `╎`
            z <= v && continue   # longer arm of `∨` is visited first
            has_edge(g, v, z) && continue
            push!(Z, (v, w, z))
        end
        for z in neighbors(g, v) # case `∧` 
            (z <= w) && continue # shorter arm is visited first
            has_edge(g, w, z) && continue
            push!(Z, (z, v, w))
        end
    end
    Z
end



alt_vskel(g) = _vskel(nv(g), dseporacle, g)

function _vskel(n::V, I, par...) where {V}
    # Step 1
    g, S = skeleton(n, I, par...)

    # Step 2: Apply Rule 0 once
    Z = orientable_unshielded(g, S)
    dg = DiGraph(g) # use g to keep track of unoriented edges

    for (u, v, w) in Z
        if has_edge(g, (u, v))
            remove!(dg, v → u)
            remove!(g, (v, u))
        end
        if has_edge(g, (v, w))
            remove!(dg, v → w)
            remove!(g, (v, w))
        end
    end
    dg
end

"""
    pcalg(g, I; stable=true)

Perform the PC algorithm for a set of 1:n variables using the tests

    I(u, v, [s1, ..., sn])

Use `IClosure(I, args)` to wrap a function f with signature

    f(u, v, [s1, ..., sn], par...)

Returns the CPDAG as DiGraph. By default uses a stable and threaded versions
of the skeleton algorithm. (This is the most recent interface)
"""
function pcalg(g::Graph, I; stable=true)
    g, S = stable ? skeleton_stable(g, I) : skeleton(g, I)
    dg = DiGraph(g) # use g to keep track of unoriented edges
    _, dg = orient_unshielded(g, dg, S)
    meek_rules!(dg)
    dg
end

"""
    pcalg(n::V, I, par...; stable=true)

Perform the PC algorithm for a set of 1:n variables using the tests

    I(u, v, [s1, ..., sn], par...)

Returns the CPDAG as DiGraph. By default uses a stable and threaded versions
of the skeleton algorithm.
"""
function pcalg(n::Integer, I, par...; stable=true)
    # Step 1
    g, S = skeleton(n, I, par...; stable)
    dg = DiGraph(g) # use g to keep track of unoriented edges
    g, dg = orient_unshielded(g, dg, S)
    #apply_pc_rules(g, dg)
    meek_rules!(dg)
end

"""
    orient_unshielded(g, dg, S)

Orient unshielded triples using the separating sets.
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
            remove!(dg, v → u)
            remove!(g, (v, u))
        end
        if has_edge(g, (v, w))
            remove!(dg, v → w)
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
    removed = Tuple{Int64, Int64}[]
    while true
        for e in edges(g)
            for e_ in (e, reverse(e))
                v, w = Tuple(e_)

                # Rule 1: Orient v-w into v->w whenever there is u->v
                # such that u and w are not adjacent
                if meek_rule1(dg, v, w)
                    VERBOSE && println("rule 1: ", v => w)
                    remove!(dg, w → v)
                    push!(removed, (w, v))
                    @goto ende
                end

                # Rule 2: Orient v-w into v->w whenever there is a chain
                # v->k->w.
                if meek_rule2(dg, v, w)
                    VERBOSE && println("rule 2: ", v => w)
                    remove!(dg, w → v)
                    push!(removed, (w, v))
                    @goto ende
                end

                # Rule 3: Orient v-w into v->w whenever there are two chains
                # v-k->w and v-l->w such that k and l are nonadjacent 
                if meek_rule3(dg, v, w)
                    VERBOSE && println("rule 3: ", v => w)
                    remove!(dg, w → v)
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

Run PC algorithm for tabular input data t using a p-value p to test for 
conditional independeces using Fisher's z-transformation.
"""
function pcalg(t, p::Float64, test::typeof(gausscitest); kwargs...)
    Tables.istable(t) || throw(ArgumentError("Argument does not support Tables.jl"))
    X = Tables.matrix(t)
    N, n = size(X)
    C = Statistics.cor(X, dims = 1)
    return pcalg(n, gausscitest, (C, N), quantile(Normal(), 1 - p / 2); kwargs...)
end

"""
    pcalg(t::T, p::Float64; cmitest::typeof(cmitest); kwargs...) where{T}

Run PC algorithm for tabular input data t using a p-value p to detect 
conditional independeces using a conditional mutual information permutation test.
"""
function pcalg(t, p::Float64, test::typeof(cmitest); stable=true, kwargs...)
    @assert Tables.istable(t)
    @assert all(t -> t == Float64, Tables.schema(t).types)

    c = Tables.columns(t)
    sch = Tables.schema(t)
    n = length(sch.names)
    cl = IClosure(cmitest, (c, p), kwargs)

    if stable
        g, S = skeleton_stable(complete_graph(n), cl)
    else 
        g, S = skeleton(complete_graph(n), cl)
    end
    dg = DiGraph(g) # use g to keep track of unoriented edges
    g, dg = orient_unshielded(g, dg, S)

    meek_rules!(dg)
end

"""
    prepare_pc_graph(g::AbstractGraph, node_labels::AbstractVector{<:AbstractString}=String[])

Prepare resulting graph for plotting with various backends.
"""
function prepare_pc_graph(g::AbstractGraph,
                          node_labels::AbstractVector{<:AbstractString} = String[])
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
                push!(styles_dict, (e.src, e.dst) => "--")
            end
        else # this is the only edge
            add_edge!(plot_g, e.src, e.dst)
            push!(styles_dict, (e.src, e.dst) => "->")
        end
    end

    (; plot_g, node_labels, edge_styles = styles_dict,
     node_style = node_style, options = options)
end

"""
    kwargs_pdag_graphmakie(g; ilabels=1:nv(g), arrowsize=25, ilabels_fontsize=25)

Generates the keywords for `GraphMakie.graphplot` to plot causal graphs and (C)PDAGs
represented as `SimpleDiGraph` as partially directed graphs.

Usage:
```
graphplot(g; kwargs_pdag_graphmakie(g)...)
```
"""
function kwargs_pdag_graphmakie(g; ilabels=1:nv(g), arrowsize=25, ilabels_fontsize=25)
    es = edges(g)
    arrow_size = Int[]
    edge_width = Dict{Graphs.SimpleEdge{Int}, Int}()
    arrow_shift = :end
    for e in es     
        if reverse(e) ∉ es 
            push!(arrow_size, arrowsize)
            continue
        end
        push!(arrow_size, 0) # remove arrows for double lines
        if e.dst < e.src
            edge_width[e] = 0 # set draw width of one of the double edges to 0
        end
    end
    
    (;curve_distance_usage=false, arrow_shift, arrow_size, edge_width, ilabels, ilabels_fontsize)
end

"""
    plot_pc_graph_text(g::AbstractGraph, node_labels::AbstractVector{<:AbstractString}=String[])

Plot the output of the PC algorithm (text-based output).

See also: `plot_pc_graph` and `plot_pc_graph_tikz` (for TikzGraphs.jl-based plotting), `plot_pc_graph_recipes` (for GraphRecipes.jl-based plotting)
"""
function plot_pc_graph_text(g::AbstractGraph,
                            node_labels::AbstractVector{<:AbstractString} = String[])
    objs = prepare_pc_graph(g, node_labels)
    graph_to_text(objs.plot_g, objs.node_labels, edge_styles = objs.edge_styles)
end

# methods to extend conditionally
"""
    plot_pc_graph_recipes(g, node_labels::AbstractVector{<:AbstractString}=String[])

Plot the output of the PC algorithm (GraphRecipes backend).

Requires GraphRecipes and Plots to be imported
"""
function plot_pc_graph_recipes end

"""
    CausalInference.plot_pc_graph_tikz(g, node_labels::AbstractVector{<:AbstractString}=String[])

Plot the output of the PC algorithm (TikzGraphs backend).

Requires TikzGraphs to be imported
"""
function plot_pc_graph_tikz end

# for backward compatibility, default to TikzGraphs when available
plot_pc_graph = plot_pc_graph_tikz
