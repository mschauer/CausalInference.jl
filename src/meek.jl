# http://proceedings.mlr.press/v89/katz19a/katz19a-supp.pdf
# https://arxiv.org/pdf/1302.4972.pdf
"""
    meek_rules!(g; rule4=false)

Apply Meek's rules 1-3 or 1-4 with rule4=true to orient edges in a 
partially directed graph without creating cycles or new v-structures.
Rule 4 is needed if edges are compelled/preoriented using external knowledge.
"""
function meek_rules!(g; rule4=false)
    while true 
        done = true
        # Loop through all the edges in the graph (u-v)
        for e in collect(edges(g)) # collect iterator as we are deleting
            u, v = Pair(e)
            # We only need to update (still) undirected edges
            !has_both(g, u, v) && continue

            # check only case u->v, we'll check v->u later
            if meek_rule1(g, u, v) || meek_rule2(g, u, v) || meek_rule3(g, u, v) || (rule4 && meek_rule4(g, u, v))
                # Make u→v
                remove!(g, v → u)
                done = false
            end
        end
        done && return g
    end
end

"""
    meek_rule1(g, v, w)

Rule 1: Orient v-w into v->w whenever there is u->v
such that u and w are not adjacent
(otherwise a new v-structure is created.)
"""
function meek_rule1(g, v, w)
    for u in inneighbors(g, v)
        u == w && continue 
        has_edge(g, v => u) && continue # not directed
        isadjacent(g, u, w) && continue
        return true
    end
    return false
end

"""
    meek_rule2(g, v, w)

Rule 2: Orient v-w into v->w whenever there is a chain v->k->w
(otherwise a directed cycle is created.)
"""
function meek_rule2(g, v, w)
    outs = Int[]
    for k in outneighbors(g, v)
        k == w && continue 
        !has_edge(g, k => v) && push!(outs, k)
    end
    ins = Int[]
    for k in inneighbors(g, w)
        k == v && continue 
        !has_edge(g, w => k) && push!(ins, k)
    end
    if !disjoint_sorted(ins, outs)
        return true
    end
    return false 
end

"""
    meek_rule3(g, v, w)

Rule 3 (Diagonal): Orient v-w into v->w whenever there are two chains
v-k->w and v-l->w such that k and l are nonadjacent
(otherwise a new v-structure or a directed cycle is created.)
"""
function meek_rule3(g, v, w)
    fulls = Int[] # Find nodes k where v-k
    for k in outneighbors(g, v)
        has_edge(g, k => v) || continue 
        # Skip if not k->w (or if not l->w)
        if has_edge(g, w => k) || !has_edge(g, k => w)
            continue
        end
        push!(fulls, k)
    end
    for (k, l) in combinations(fulls, 2) # FIXME: 
        isadjacent(g, k, l) && continue
        return true
    end
    return false
end

"""
    meek_rule4(g, v, w)

Rule 4: Orient v-w into v→w if v-k→l→w where adj(v,l) and not adj(k,w) [check].
"""
function meek_rule4(g, v, w)
    for l in inneighbors(g, w)
        has_edge(g, w => l) && continue # undirected
        !isadjacent(g, v, l) && continue # not adjacent to v      
        for k in inneighbors(g, l)
            has_edge(g, l => k) && continue # undirected
            !has_both(g, v, k) && continue # not undirected to v
            isadjacent(g, k, w) && continue # adjacent to w
            return true
        end
    end
    return false
end

"""
    pdag_to_dag_meek!(g, rule4=false)

Complete PDAG to DAG using meek_rules.
"""
function pdag_to_dag_meek!(g, rule4=false)
    while true
        # find unoriented edge
        for e in edges(g) # go through edges (bad to start in the beginning?)
            x, y = Pair(e)
            if has_edge(g, y → x)
                orientedge!(g, rand(Bool) ? y → x : x → y)
                @goto orient
            end
        end
        return g # through without break 
        # orient implied edges
        @label orient
        meek_rules!(g; rule4)
    end
    g
end
"""
Deprecated alias for `pdag_to_dag_meek!`.
"""
const pdag2dag! = pdag_to_dag_meek!