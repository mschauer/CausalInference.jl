using Base.Test
using LightGraphs

Base.start(e::LightGraphs.SimpleGraphs.SimpleEdge) = 1
Base.next(e::LightGraphs.SimpleGraphs.SimpleEdge, state) = state == 1 ? (src(e), 2) : (dst(e), 3)
Base.done(e::LightGraphs.SimpleGraphs.SimpleEdge, state) = state == 3

insorted(a, x) = !isempty(searchsorted(a, x))

function orient_unshielded(g, S)
    Z = []
    for e in edges(g)
        v, w = (src(e), dst(e))
        assert(v < w)
        for z in neighbors(g, w)
            z <= w && continue
            insorted(neighbors(g, z), v) && continue
            w in S[Edge(v,z)] || push!(Z, (v, w, z))
        end
        for z in neighbors(g, v)
            (z <= v || z == w) && continue
            insorted(neighbors(g, z), w) && continue
            w in S[Edge(minmax(z,w)...)] || push!(Z, (z, v, w))
        end
    end
    Z
end

@time Z = orient_unshielded(h, s)
g = h

for z in Z
    u, v, w = z
    n = neighbors(g, v)
    @test u in n
    @test w in n
    @test !(u in neighbors(g, w))
end