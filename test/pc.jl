using CausalInference
using LightGraphs
using Base.Test

g = DiGraph(5)
d = nv(g)
for (i,j) in [(1,2), (2,3), (2,4),(4,5), (3,5)]
   add_edge!(g,i,j)
end

h, s = skeleton(d, dseporacle, g)
@test Graph(g) == h

Z = unshielded(h, s)

for z in Z
    u, v, w = z
    n = in_neighbors(h, v)
    @test u in n
    @test w in n
    @test !(u in neighbors(h, w))
end

dg = pcalg(d, dseporacle, g)
