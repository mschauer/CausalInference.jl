using LightGraphs, MetaGraphs
using Combinatorics: combinations, powerset

function has_marks(dg, v1, v2, s::String)
    symbols = ['*', 'o', '>', '<', '-']
    symDict = Dict('o' => :circle,
                   '>' => :arrow,
                   '<' => :arrow,
                   '-' => :tail)
    
    @assert length(s)==3
    @assert s[1] ∈ symbols
    @assert s[2] == '-'
    @assert s[3] ∈ symbols

    if s[1]=='*'
        check1 = false 
    else
        check1 = true
        mark1 = symDict[s[1]]
    end

    if s[3]=='*'
        check2 = false
    else
        check2 = true
        mark2 = symDict[s[3]]
    end
    
    if check2
        if check1
            result = get_prop(dg, v2, v1, :mark)==mark1 && get_prop(dg, v1, v2, :mark)==mark2
        else
            result = get_prop(dg, v1, v2, :mark)==mark2
        end
    else
        result = get_prop(dg, v2, v1, :mark)==mark1
    end
    
    return result
end

function set_marks!(dg, v1, v2, s::String)
 symbols = ['*', 'o', '>', '<', '-']
    symDict = Dict('o' => :circle,
                   '>' => :arrow,
                   '<' => :arrow,
                   '-' => :tail)
    
    @assert length(s)==3
    @assert s[1] ∈ symbols
    @assert s[2] == '-'
    @assert s[3] ∈ symbols

    if s[1]!='*'
        set_prop!(dg, v2, v1, :mark, symDict[s[1]])
    end

    if s[3]!='*'
        set_prop!(dg, v1, v2, :mark, symDict[s[3]])
    end
end

function is_collider(dg, v1, v2, v3)
    return has_marks(dg, v1, v2, "*->") && has_marks(dg, v2, v3, "<-*")
end

function is_parent(dg, v1, v2)
    return has_edge(dg, v1, v2) && has_marks(dg, v1, v2, "-->")
end

function is_triangle(dg, v1, v2, v3)
    return isadjacent(dg, v1, v2) && isadjacent(dg, v2, v3) && isadjacent(dg, v3, v1)
end

function is_discriminating_path(dg, path)
    if length(path)<4 || isadjacent(dg, path[1], path[end])
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

function isUncoveredCirclePath(dg, path)
    if(length(path)<3)
        return has_marks(dg, path[1], path[2], "o-o") 
    end
    
    edges = collect(zip(path[1:end-1], path[2:end]))
    triples = collect(zip(path[1:end-2], path[2:end-1], path[3:end]))
    unshielded = map(t->!isadjacent(dg, t[1], t[3]), triples)
    circles = map(e->has_marks(dg, e[1], e[2], "o-o"), edges)

    return all(unshielded) && all(circles)
end

function isUncoveredPDPath(dg, path)
    if(length(path)<3)
        return (!has_marks(dg, path[1], path[2], "<-*") &&
                !has_marks(dg, path[1], path[2], "*--"))
    end
    
    edges = collect(zip(path[1:end-1], path[2:end]))
    triples = collect(zip(path[1:end-2], path[2:end-1], path[3:end]))
    unshielded = map(t->!isadjacent(dg, t[1], t[3]), triples)
    directions = map(e->(!has_marks(dg, e[1], e[2], "<-*") &&
                         !has_marks(dg, e[1], e[2], "*--")), edges)

    return all(unshielded) && all(directions)
end

function fcialg(n::V, I, par...; kwargs...) where {V<:Integer}

    # Step F1 and F2
    g, S = skeleton(n, I, par...; kwargs...)

    # Apply R0 once
    Z = unshielded(g, S)
    dg = MetaDiGraph(g) # use g to keep track of unoriented edges

    # construct initial PAG
    for e in edges(dg)
        set_prop!(dg, e, :mark, :circle)
    end
        
    for (u, v, w) in Z
        if has_edge(dg, (u, v))
            set_marks!(dg, u, v, "*->")
        end
        if has_edge(dg, (v, w))
            set_marks!(dg, v, w, "<-*")
        end
    end
    
    # find possible d-separation sets
    pdsep = Dict()
    
    for v in vertices(g)
        pdsep[v] = Set{Int64}()
        for w in vertices(g)
            if w == v
                continue
            end
            paths = yen_k_shortest_paths(g, v, w, LightGraphs.weights(g), 100).paths
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
    Z = unshielded(g, S)
    dg = MetaDiGraph(g)

    for e in edges(dg)
        set_marks!(dg, src(e), dst(e), "o-o")
    end
    
    for (u, v, w) in Z
        if has_edge(dg, (u, v))
            set_marks!(dg, u, v, "*->")
        end
        if has_edge(dg, (v, w))
            set_marks!(dg, v, w, "<-*")
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
                    has_marks(dg, α, β, "*->") &&
                    has_marks(dg, β, γ, "o-*"))

                    set_marks!(dg, β, γ, "-->")
                    loop = true

                end
                
                # R2
                if (has_edge(dg, α, γ) &&
                    has_marks(dg, α, γ, "*-o") &&
                    ((has_marks(dg, α, β, "-->") && has_marks(dg, β, γ, "o-*")) ||
                     (has_marks(dg, α, β, "o->") && has_marks(dg, β, γ, "-->"))))
                    
                    set_marks!(dg, α, γ, "*->")
                    loop = true
                end
                
                #R3
                if !isadjacent(dg, α, γ)
                    for θ in inneighbors(dg, γ)
                        if (θ ∈ inneighbors(dg, α) &&
                            θ ∈ inneighbors(dg, β) &&
                            has_marks(dg, α, β, "*->") &&
                            has_marks(dg, β, γ, "<-*") &&
                            has_marks(dg, α, θ, "*-o") &&
                            has_marks(dg, γ, θ, "o-*") &&
                            has_marks(dg, θ, β, "*-o"))
                            set_marks!(dg, θ, β, "*->")
                            loop = true
                        end
                    end
                end

            end

            # R4
            for x in vertices(dg)
                paths = yen_k_shortest_paths(g, x, α, LightGraphs.weights(g), 100).paths
                for path in paths
                    if (is_discriminating_path(dg, path) &&
                        has_marks(dg, path[end-1], path[end], "o-*"))
                        if (haskey(S, Edge(path[1], path[end])) &&
                            path[end-1] ∈ S[Edge(path[1], path[end])])
                            set_marks!(dg, path[end-1], path[end], "-->")
                        else
                            set_marks!(dg, path[end-1], path[end], "<->")
                            set_marks!(dg, path[end-2], path[end-1], "<->")
                        end
                        loop=true
                    end
                end
            end
        end

    end

    # rules R5 to R10

    loop = true
    while loop
        loop = false
        for e in edges(dg)
            (α, β) = Tuple(e)

            # R5
            if has_marks(dg, α, β, "o-o")
                paths = yen_k_shortest_paths(g, α, β, LightGraphs.weights(g), 100).paths
                
                for path in paths
                    if (isUncoveredCirclePath(dg, path) &&
                        !isadjacent(dg, path[1], path[end-1]) &&
                        !isadjacent(dg, path[2], path[end]))
                        
                        set_marks!(dg, α, β, "---")
                        for (e1, e2) in zip(path[1:end-1], path[2:end])
                            set_marks!(dg, e1, e2, "---")
                        end
                        loop = true
                        break
                    end
                end
            end

            for γ in inneighbors(dg, β)
                #R6
                if has_marks(dg, α, β, "---") && has_marks(dg, β, γ, "o-*")
                    set_marks!(dg, β, γ, "--*")
                    loop = true
                end

                # R7
                if(!isadjacent(dg, α, γ) &&
                   has_marks(dg, α, β, "--o") &&
                   has_marks(dg, β, γ, "o-*"))
                    set_marks!(dg, β, γ, "--*")
                    loop = true
                end

                # R8
                if(has_edge(dg, α, γ) &&
                   has_marks(dg, α, γ, "o->") &&
                   (has_marks(dg, α, β, "-->") || has_marks(dg, α, β, "--o")) &&
                   has_marks(dg, β, γ, "-->"))
                    set_marks!(dg, α, γ, "-->")
                    loop = true
                end
            end
            
            # R9
            if has_marks(dg, α, β, "o->")
                paths = yen_k_shortest_paths(g, α, β, LightGraphs.weights(g), 100).paths
                for path in paths
                    if (length(path)>3 &&
                        isUncoveredPDPath(dg, path) &&
                        !isadjacent(dg, path[2], path[end]))
                        set_marks!(dg, α, β, "-->")
                        loop = true
                        break
                    end
                end
            end
            
            #R10
            if has_marks(dg, α, β, "o->")
                for (γ,θ) in combinations(inneighbors(dg, β), 2)
                    if(θ==α || γ==α)
                        continue
                    end
                    if(has_marks(dg, γ, β, "-->") &&
                       has_marks(dg, β, θ, "<--"))

                        p1 = yen_k_shortest_paths(g, α, γ, LightGraphs.weights(g), 100).paths
                        p2 = yen_k_shortest_paths(g, α, θ, LightGraphs.weights(g), 100).paths

                        for path1 in p1
                            if isUncoveredPDPath(dg, p1)
                                for path2 in p2
                                    if isUncoveredPDPath(dg, p2)
                                        μ = path1[2]
                                        ω = path2[2]
                                        if(μ != ω && !isadjacent(dg, μ, ω))
                                            set_marks!(dg, α, β, "-->")
                                            loop=true
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
  
    X = reduce(hcat, map(c->t[c], 1:n))
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
