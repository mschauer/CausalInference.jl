using Graphs, MetaGraphs
using Combinatorics: combinations, powerset

macro arrow_str(str)
    symbols = ['*', 'o', '>', '<', '-']
    symDict = Dict('o' => :circle,
                   '>' => :arrow,
                   '<' => :arrow,
                   '-' => :tail,
                   '*' => :star)
    
    @assert length(str)==3
    @assert str[1] ∈ symbols
    @assert str[2] == '-'
    @assert str[3] ∈ symbols

    return (symDict[str[1]], symDict[str[3]])
end

"""
    has_marks(dg, v1, v2, t::Tuple{Symbol, Symbol}

test if the edge between node v1 and v2 has the edge markers given by the tuple t (use
the arrow macro to simplify use)

Example:
has_marks(dg, 1, 2, arrow"o->")
"""
function has_marks(dg, v1, v2, t::Tuple{Symbol, Symbol})    
    if t[2]!=:star
        if t[1]!=:star
            result = get_prop(dg, v2, v1, :mark)==t[1] && get_prop(dg, v1, v2, :mark)==t[2]
        else
            result = get_prop(dg, v1, v2, :mark)==t[2]
        end
    else
        result = get_prop(dg, v2, v1, :mark)==t[1]
    end
    
    return result
end

"""
    set_marks!(dg, v1, v2, t::Tuple{Symbol, Symbol})

set edge marks between node v1 and v2.

Example:
set_marks!(dg, 1, 2, arrow"*->")
"""
function set_marks!(dg, v1, v2, t::Tuple{Symbol, Symbol})
    if t[1]!=:star
        set_prop!(dg, v2, v1, :mark, t[1])
    end

    if t[2]!=:star
        set_prop!(dg, v1, v2, :mark, t[2])
    end
end

"""
    is_collider(dg, v1, v2, v3)

check if egde v1, v2 and v3 form a collider
"""
function is_collider(dg, v1, v2, v3)
    return has_marks(dg, v1, v2, arrow"*->") && has_marks(dg, v2, v3, arrow"<-*")
end

"""
    is_parent(dg, v1, v2)

check if v1 is a parent of v2
"""
function is_parent(dg, v1, v2)
    return has_edge(dg, v1, v2) && has_marks(dg, v1, v2, arrow"-->")
end

"""
    is_triangle(dg, v1, v2, v3)

check if v1, v2 and v3 form a triangle
"""
function is_triangle(dg, v1, v2, v3)
    return isadjacent(dg, v1, v2) && isadjacent(dg, v2, v3) && isadjacent(dg, v3, v1)
end

"""
    is_discriminating_path(dg, path)

check if `path` is a discriminating path
"""
function is_discriminating_path(dg, path)
    # a discriminating path consists of at least four edges
    if length(path)<4
        return false
    end
    
    triples = collect(zip(path[1:end-2], path[2:end-1], path[3:end]))[1:end-1]
    colliders = map(t->is_collider(dg, t...), triples)
    parents = map(t->is_parent(dg, t[2], path[end]), triples)
    if all(colliders) && all(parents)
        return true
    else
        return false
    end
end

"""
    is_uncovered_circle_path(dg, path)

check if `path` is an uncovered circle path
"""
function is_uncovered_circle_path(dg, path)
    if length(path)<4 
        return false
    end
    
    edges = collect(zip(path[1:end-1], path[2:end]))
    triples = collect(zip(path[1:end-2], path[3:end]))
    
    unshielded = map(t->!isadjacent(dg, t[1], t[2]), triples)
    circles = map(e->has_marks(dg, e[1], e[2], arrow"o-o"), edges)
    
    return all(unshielded) && all(circles)
end

"""
    is_uncovered_PD_path(dg, path)

check if `path` is an uncovered potentially directed path
"""
function is_uncovered_PD_path(dg, path)
    if length(path)<4
        return false
    end
    
    edges = collect(zip(path[1:end-1], path[2:end]))
    triples = collect(zip(path[1:end-2], path[3:end]))
    unshielded = map(t->!isadjacent(dg, t[1], t[2]), triples)
    directions = map(e->(!has_marks(dg, e[1], e[2], arrow"<-*") &&
                         !has_marks(dg, e[1], e[2], arrow"*--")), edges)

    return all(unshielded) && all(directions)
end

"""
    fcialg(n::V, I, par...; augmented=true, verbose=false, kwargs...)

Perform the FCI algorithm for a set of `n` variables using the test

    I(u, v, [s1, ..., sn], par...)

Returns the PAG as a MetaDiGraph
"""
function fcialg(n::V, I, par...; augmented=true, verbose=false, kwargs...) where {V<:Integer}

    # Step F1 and F2
    g, S = skeleton(n, I, par...; kwargs...)

    # Apply R0 once
    Z = orientable_unshielded(g, S)
    dg = MetaDiGraph(g) # use g to keep track of unoriented edges

    # construct initial PAG
    for e in edges(dg)
        set_prop!(dg, e, :mark, :circle)
    end

    for (u, v, w) in Z
        # is this check actually needed?
        if has_edge(dg, (u, v))
            set_marks!(dg, u, v, arrow"*->")
        end
        if has_edge(dg, (v, w))
            set_marks!(dg, v, w, arrow"<-*")
        end
    end
    
    # find possible d-separation sets (MAGs are trickier than DAGs...)
    pdsep = Dict()
    
    for v in vertices(g)
        pdsep[v] = Set{Int64}()
        for w in vertices(g)
            if w == v
                continue
            end
            paths = yen_k_shortest_paths(g, v, w, Graphs.weights(g), 100).paths
            for path in paths
                if length(path)==2
                    push!(pdsep[v], w)
                else
                    triples = zip(path[1:end-2], path[2:end-1], path[3:end])
                    tconds = map(t->is_triangle(dg, t...) || is_collider(dg, t...), triples)
                    if all(tconds)
                        push!(pdsep[v], w)
                    end
                end
            end
        end
    end

    # need to collect edges here since graph could change while looping
    for e in collect(edges(g))
        v = src(e)
        w = dst(e)
        if !(has_edge(g, v, w) || has_edge(g, v, w))
            # edge has already been removed
            continue
        end
        seps1 = Set(powerset([d for d in pdsep[v] if d!=w]))
        seps2 = Set(powerset([d for d in pdsep[w] if d!=v]))
        seps = union(seps1, seps2)

        for s in seps
            if I(src(e), dst(e), s, par...; kwargs...)
                @debug "Found hidden dsep set: $(src(e)) - $(dst(e)) given $(s)"
                rem_edge!(g, e)
                if !(e in keys(S))
                    S[e] = s
                end
                break
            end
        end
    end

    #  step F3
    Z = orientable_unshielded(g, S)
    dg = MetaDiGraph(g)

    for e in edges(dg)
        set_marks!(dg, src(e), dst(e), arrow"o-o")
    end
    
    for (u, v, w) in Z
        verbose && println("R0 with ($(u), $(v), $(w))")
        if has_edge(dg, (u, v))
            set_marks!(dg, u, v, arrow"*->")
        end
        if has_edge(dg, (v, w))
            set_marks!(dg, v, w, arrow"<-*")
        end
    end

    # main loop for rules R1 to R4
    # loop is repeated until none of them apply
    loop = true
    while loop
        loop = false
        for e in edges(dg)
            (α, β) = Tuple(e)
            
            for γ in inneighbors(dg, β)
                if γ==α
                    continue
                end

                # R1
                if (!has_edge(dg, α, γ) &&
                    has_marks(dg, α, β, arrow"*->") &&
                    has_marks(dg, β, γ, arrow"o-*"))

                    set_marks!(dg, β, γ, arrow"-->")
                    loop = true
                    verbose && println("R1 with $(α)-$(β)-$(γ)")
                end
                
                # R2
                if (has_edge(dg, α, γ) &&
                    has_marks(dg, α, γ, arrow"*-o") &&
                    ((has_marks(dg, α, β, arrow"-->") && has_marks(dg, β, γ, arrow"o-*")) ||
                     (has_marks(dg, α, β, arrow"o->") && has_marks(dg, β, γ, arrow"-->"))))
                    
                    set_marks!(dg, α, γ, arrow"*->")
                    loop = true
                    verbose && println("R1 with $(α)-$(β)-$(γ)")
                end
                
                #R3
                if !isadjacent(dg, α, γ)
                    for θ in inneighbors(dg, γ)
                        if (θ ∈ inneighbors(dg, α) &&
                            θ ∈ inneighbors(dg, β) &&
                            has_marks(dg, α, β, arrow"*->") &&
                            has_marks(dg, β, γ, arrow"<-*") &&
                            has_marks(dg, α, θ, arrow"*-o") &&
                            has_marks(dg, γ, θ, arrow"o-*") &&
                            has_marks(dg, θ, β, arrow"*-o"))
                            set_marks!(dg, θ, β, arrow"*->")
                            loop = true
                            verbose && println("R3 with $(α)-$(β)-$(γ)")
                        end
                    end
                end

            end

            # R4
            for x in vertices(dg)
                paths = yen_k_shortest_paths(g, x, α, Graphs.weights(g), 100).paths
                for path in paths
                    if (is_discriminating_path(dg, path) &&
                        has_marks(dg, path[end-1], path[end], arrow"o-*"))
                        if (haskey(S, Edge(path[1], path[end])) &&
                            path[end-1] ∈ S[Edge(path[1], path[end])])
                            set_marks!(dg, path[end-1], path[end], arrow"-->")
                        else
                            set_marks!(dg, path[end-1], path[end], arrow"<->")
                            set_marks!(dg, path[end-2], path[end-1], arrow"<->")
                        end
                        loop=true
                        verbose && println("R4 with $(α) and $(path)")
                    end
                end
            end
        end

    end

    if !augmented
        return dg
    end
    
    # rules R5 to R10

    loop = true
    while loop
        loop = false
        for e in edges(dg)
            (α, β) = Tuple(e)

            # R5
            if has_marks(dg, α, β, arrow"o-o")
                paths = yen_k_shortest_paths(g, α, β, Graphs.weights(g), 100).paths
                
                for path in paths
                    if (is_uncovered_circle_path(dg, path) &&
                        !isadjacent(dg, path[1], path[end-1]) &&
                        !isadjacent(dg, path[2], path[end]))
                        verbose && println("R5: $(α)-$(β) with $(path)")
                        set_marks!(dg, α, β, arrow"---")
                        for (e1, e2) in zip(path[1:end-1], path[2:end])
                            set_marks!(dg, e1, e2, arrow"---")
                        end
                        loop = true
                        break
                    end
                end
            end

            for γ in inneighbors(dg, β)

                if γ == α
                    continue
                end
                
                #R6
                if has_marks(dg, α, β, arrow"---") && has_marks(dg, β, γ, arrow"o-*")
                    set_marks!(dg, β, γ, arrow"--*")
                    loop = true
                    verbose && println("R6: $(α)-$(β)-$(γ)")
                end

                # R7
                if(!isadjacent(dg, α, γ) &&
                   has_marks(dg, α, β, arrow"--o") &&
                   has_marks(dg, β, γ, arrow"o-*"))
                    set_marks!(dg, β, γ, arrow"--*")
                    loop = true
                    verbose && println("R7: $(α)-$(β)-$(γ)")
                end

                # R8
                if(has_edge(dg, α, γ) &&
                   has_marks(dg, α, γ, arrow"o->") &&
                   (has_marks(dg, α, β, arrow"-->") || has_marks(dg, α, β, arrow"--o")) &&
                   has_marks(dg, β, γ, arrow"-->"))
                    set_marks!(dg, α, γ, arrow"-->")
                    loop = true
                    verbose && println("R8: $(α)-$(β)-$(γ)")
                end
            end
            
            # R9
            if has_marks(dg, α, β, arrow"o->")
                paths = yen_k_shortest_paths(g, α, β, Graphs.weights(g), 100).paths
                for path in paths
                    if (length(path)>3 &&
                        is_uncovered_PD_path(dg, path) &&
                        !isadjacent(dg, path[2], path[end]))
                        set_marks!(dg, α, β, arrow"-->")
                        loop = true
                        verbose && println("R9: $(α)-$(β) with $(path)")
                        break
                    end
                end
            end
            
            #R10
            if has_marks(dg, α, β, arrow"o->")
                for (γ,θ) in combinations(inneighbors(dg, β), 2)
                    if(θ==α || γ==α)
                        continue
                    end
                    if(has_marks(dg, γ, β, arrow"-->") &&
                       has_marks(dg, β, θ, arrow"<--"))

                        p1 = yen_k_shortest_paths(g, α, γ, Graphs.weights(g), 100).paths
                        p2 = yen_k_shortest_paths(g, α, θ, Graphs.weights(g), 100).paths

                        for path1 in p1
                            if is_uncovered_PD_path(dg, path1)
                                for path2 in p2
                                    if is_uncovered_PD_path(dg, path2)
                                        μ = path1[2]
                                        ω = path2[2]
                                        if(μ != ω && !isadjacent(dg, μ, ω))
                                            set_marks!(dg, α, β, arrow"-->")
                                            loop=true
                                            verbose && println("R10: $(α)-$(β) with $(p1) and $(p2)")
                                            break
                                        end
                                    end
                                end
                            end
                        end
                    end                            
                end
            end
        end
        
    end
    dg
end


function fcialg(t, p::Float64, test::typeof(gausscitest); kwargs...)
    @assert Tables.istable(t)

    c = Tables.columns(t)
    sch = Tables.schema(t)
    n = length(sch.names)
  
    X = reduce(hcat, map(c->Tables.getcolumn(Tables.columns(t), c), Tables.columnnames(t)))
    N = size(X,1)
    C = Statistics.cor(X)
    return fcialg(n, gausscitest, (C,N), quantile(Normal(), 1-p/2); kwargs...)
end

function fcialg(t, p::Float64, test::typeof(cmitest); kwargs...)
    @assert Tables.istable(t)
    @assert all(t->t==Float64, Tables.schema(t).types)

    c = Tables.columns(t)
    sch = Tables.schema(t)
    n = length(sch.names)

    return fcialg(n, cmitest, c, p; kwargs...)
end

"""
    plot_fci_graph(g, node_labels)

plot the output of the FCI algorithm.
"""
function plot_fci_graph(g, node_labels::Array=[])
    plot_g = DiGraph(nv(g))

    if length(node_labels) != nv(g) 
        node_labels = map(string, 1:nv(g))
    end
        
    node_style = "draw, rounded corners, fill=blue!10"
    options = "scale=2"
    
    styles_dict = Dict()
        
    for e in edges(g)
       if e.src < e.dst
           add_edge!(plot_g, e.src, e.dst)

           if has_marks(g, e.src, e.dst, arrow"o-o")
               push!(styles_dict, (e.src, e.dst)=>"o-o")
           elseif has_marks(g, e.src, e.dst, arrow"o->")
               push!(styles_dict, (e.src, e.dst)=>"o->")
           elseif has_marks(g, e.src, e.dst, arrow"<-o")
               push!(styles_dict, (e.src, e.dst)=>"<-o")
           elseif has_marks(g, e.src, e.dst, arrow"-->")
               push!(styles_dict, (e.src, e.dst)=>"->")
           elseif has_marks(g, e.src, e.dst, arrow"<--")
               push!(styles_dict, (e.src, e.dst)=>"<-")
           elseif has_marks(g, e.src, e.dst, arrow"---")
               push!(styles_dict, (e.src, e.dst)=>"--")
           end
            
        end
    end
    
    TikzGraphs.plot(plot_g, node_labels, edge_styles=styles_dict,
                    node_style=node_style, options=options)
end
