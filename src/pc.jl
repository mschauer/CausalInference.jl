using LightGraphs

Base.start(e::LightGraphs.SimpleGraphs.SimpleEdge) = 1
Base.next(e::LightGraphs.SimpleGraphs.SimpleEdge, state) = state == 1 ? (src(e), 2) : (dst(e), 3)
Base.done(e::LightGraphs.SimpleGraphs.SimpleEdge, state) = state == 3

insorted(a, x) = !isempty(searchsorted(a, x))

"""
    unshielded(g, S)

Find unshielded triples in a graph. Triples are connected vertices v-w-z where
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

function pcalg(n::V, I, par...) where {V}
    g, S = skeleton(n, I, par...)
    z = unshielded(g, S)
    dg = DiGraph(g)
    for (u, v, w) in Z
        rem_edge!(g, (v, u))
        rem_edge!(g, (v, w))        
    end
    g

    #todo
end
