using LightGraphs
import LightGraphs.rem_edge!

Base.start(e::LightGraphs.SimpleGraphs.SimpleEdge) = 1
Base.next(e::LightGraphs.SimpleGraphs.SimpleEdge, state) = state == 1 ? (src(e), 2) : (dst(e), 3)
Base.done(e::LightGraphs.SimpleGraphs.SimpleEdge, state) = state == 3

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
    su = start(u)
    sv = start(v)
    done(u, su) && return true
    done(v, sv) && return true
    x, su = next(u, su)
    y, sv = next(v, sv)
    while true
        x == y && return false
        if x > y
            done(v, sv) && return true
            y, sv = next(v, sv)
        else
            done(u, su) && return true
            x, su = next(u, su)
        end
    end
end


"""
    unshielded(g, S)

Find unshielded triples in the skeleton. Triples are connected vertices v-w-z where
z is not a neighbour of v. Uses that `edges` iterates in lexicographical order.
""" 
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
function unshielded(g, S)
    Z = Tuple{Int64,Int64,Int64}[]
    for e in edges(g)
        v, w = (src(e), dst(e))
        assert(v < w)
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

adjacent(dg::DiGraph, v, w) = has_edge(dg, (v, w)) || has_edge(dg, (w, v))
has_both(dg::DiGraph, v, w) = has_edge(dg, (v, w)) && has_edge(dg, (w, v))

rem_edge!(dg::DiGraph, e::Pair) = rem_edge!(dg, Edge(e))
rem_edge!(dg::Graph, e::Tuple) = rem_edge!(dg, Edge(e))

"""
    vskel(g)

Skeleton and `v`-structures. (Currently from the first step of the pc-Alg.)
"""
vskel(g) = _vskel(nv(g), dseporacle, g)

function _vskel(n::V, I, par...) where {V}
    # Step 1
    g, S = skeleton(n, I, par...)

    # Step 2: Apply Rule 0 once
    Z = unshielded(g, S)
    dg = DiGraph(g) # use g to keep track of unoriented edges

    for (u, v, w) in Z
        rem_edge!(dg, v => u)
        rem_edge!(dg, v => w)

        rem_edge!(g, (v, u))
        rem_edge!(g, (v, w))
    end
    dg
end

"""
    pcalg(n::V, I, par...)

Perform the PC skeleton algorithm for a set of 1:n variables using the tests

    I(u, v, [s1, ..., sn], par...)

Returns the CPDAG as DiGraph.   
"""
function pcalg(n::V, I, par...) where {V}
    const VERBOSE = false

    # Step 1
    g, S = skeleton(n, I, par...)

    # Step 2: Apply Rule 0 once
    Z = unshielded(g, S)
    dg = DiGraph(g) # use g to keep track of unoriented edges

    for (u, v, w) in Z
        rem_edge!(dg, v => u)
        rem_edge!(dg, v => w)

        rem_edge!(g, (v, u)) 
        rem_edge!(g, (v, w))
    end

    # Step 3: Apply Rule 1-3 consecutively
    removed = Tuple{Int64,Int64}[]
    while true
        for e in edges(g)
            for (v, w) in (e, reverse(e))
                # Rule 1: Orient v-w into v->w whenever there is u->v
                # such that u and w are not adjacent
                for u in in_neighbors(dg, v)
                    has_edge(dg, v => u) && continue # not directed
                    adjacent(dg, u, w) && continue
                    VERBOSE && println("rule 1: ", v => w) 
                    rem_edge!(dg, w => v)
                    push!(removed, (w, v))
                    @goto ende
                end
            

                # Rule 2: Orient v-w into v->w whenever there is a chain
                # v->k->w.
                outs = Int[]
                for k in out_neighbors(dg, v)
                    !has_edge(dg, k => v) && push!(outs, k)
                end
                ins = Int[]
                for k in in_neighbors(dg, w)
                    !has_edge(dg, w => k) && push!(ins, k)
                end
                
                if !disjoint_sorted(ins, outs)
                    VERBOSE && println("rule 2: ", v => w) 
                    rem_edge!(dg, w => v)
                    push!(removed, (w, v))
                    @goto ende
                end

                # Rule 3: Orient v-w into v->w whenever there are two chains
                # v-k->w and v-l->w such that k and l are nonadjacent
                fulls = [] # Find nodes k where v-k
                for k in out_neighbors(dg, v)
                    has_edge(dg, k => v) && push!(fulls, k)
                end
                for (k, l) in combinations(fulls, 2) # FIXME: 
                    adjacent(dg, k, l) && continue
                    
                    # Skip if not k->w or if not l->w
                    if has_edge(dg, w => k) || !has_edge(dg, k => w)
                        continue
                    end
                    if has_edge(dg, w => l) || !has_edge(dg, l => w)
                        continue
                    end
                    VERBOSE && println("rule 3: ", v => w) 
                    rem_edge!(dg, w => v)
                    push!(removed, (w, v))
                    @goto ende
                end
            end
        end
        
        @label ende
        for e in removed
            rem_edge!(g, e)
        end
        isempty(removed) && break 
        empty!(removed)
    end
    dg
end
