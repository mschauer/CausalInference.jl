# http://proceedings.mlr.press/v89/katz19a/katz19a-supp.pdf
# https://arxiv.org/pdf/1302.4972.pdf
function meek_rules!(g)
    while true 
        done = true
        # Loop through all the edges in the graph (u-v)
        for e in collect(edges(g)) # collect iterator as we are deleting
            (u, v) = src(e), dst(e)
            # We only need to update undirected edges
            !has_edge(g, reverse(e)) && continue
            # check only case u->v, well check v->u later
        
            if meek_rule1(g, u, v) || meek_rule2(g, u, v) || meek_rule3(g, u, v) || meek_rule4(g, u, v)
                # Make u→v
                remove!(g, v => u)
                done = false
            end
        end
        done && return 
    end
end

"""
    meek_rule1(dg, v, w)

Rule 1: Orient v-w into v->w whenever there is u->v
such that u and w are not adjacent
(otherwise a new v-structure is created.)
"""
function meek_rule1(dg, v, w)
    for u in inneighbors(dg, v)
        has_edge(dg, v => u) && continue # not directed
        isadjacent(dg, u, w) && continue
        return true
    end
    return false
end

"""
    meek_rule2(dg, v, w)

Rule 2: Orient v-w into v->w whenever there is a chain v->k->w
(otherwise a directed cycle is created.)
"""
function meek_rule2(dg, v, w)
    outs = Int[]
    for k in outneighbors(dg, v)
        !has_edge(dg, k => v) && push!(outs, k)
    end
    ins = Int[]
    for k in inneighbors(dg, w)
        !has_edge(dg, w => k) && push!(ins, k)
    end
    if !disjoint_sorted(ins, outs)
        return true
    end
    return false 
end

"""
    meek_rule3(dg, v, w)

Rule 3 (Diagonal): Orient v-w into v->w whenever there are two chains
v-k->w and v-l->w such that k and l are nonadjacent
(otherwise a new v-structure or a directed cycle is created.)
"""
function meek_rule3(dg, v, w)
    fulls = [] # Find nodes k where v-k
    for k in outneighbors(dg, v)
        has_edge(dg, k => v) || continue 
        # Skip if not k->w (or if not l->w)
        if has_edge(dg, w => k) || !has_edge(dg, k => w)
            continue
        end
        push!(fulls, k)
    end
    for (k, l) in combinations(fulls, 2) # FIXME: 
        isadjacent(dg, k, l) && continue
        return true
    end
    return false
end

"""
    meek_rule4(dg, v, w)

Rule 4: Orient v-w into v→w if v-k→l→w where v--l and not adj(k,w) [check].
Used in GES but not used in PC algorithm?
"""
function meek_rule4(dg, v, w)
    for l in inneighbors(dg, w)
        has_edge(dg, w => l) && continue # undirected
        !isadjacent(dg, v, l) && continue # not adjacent to v      
        for k in inneighbors(dg, l)
            has_edge(dg, l => k) && continue # undirected
            !has_both(dg, v, k) && continue # not undirected to v
            isadjacent(dg, k, w) && continue # adjacent to w
            return true
        end
    end
    return false
end


