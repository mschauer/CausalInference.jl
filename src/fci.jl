using LightGraphs, MetaGraphs
using Combinatorics: powerset

function is_collider(dg, v1, v2, v3)
    return has_edge(dg, v1, v2) && has_edge(dg, v3, v2)
end


function is_triangle(dg, v1, v2, v3)
    return isadjacent(dg, v1, v2) && isadjacent(dg, v2, v3) && isadjacent(dg, v3, v1)
end


function fcialg(n::V, I, par...; kwargs...) where {V<:Integer}

    # Step 1
    g, S = skeleton(n, I, par...; kwargs...)

    # Step 2: Apply Rule 0 once
    Z = unshielded(g, S)
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

    println(collect(edges(g)))
    
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
        println("$(v): $(pdsep[v])")
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
        println(seps)
        for s in seps
            if I(src(e), dst(e), s, par...; kwargs...)
                println("Found hidden dsep set:")
                println("$(src(e)) - $(dst(e))")
                println(s)
                rem_edge!(g, e)
                if !(e in keys(S))
                    S[e] = s
                end
                break
            end
        end
    end
    g
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
