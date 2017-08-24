using LightGraphs
using LightGraphs.SimpleGraphs
using  Combinatorics


using PlotRecipes
using Base.Test
gplot(g) = graphplot(g, names=vertices(g), nodesize=2.5, fontsize=20, nodeshape=:circle)


 
"""
    dsep(g::AbstractGraph, u, v, s; verbose = false)
    
    
Algorithm:

"""
function dsep(g::AbstractGraph, u::Integer, v::Integer, S; verbose = false)
    T = eltype(g)
    in_seen = falses(nv(g)) # nodes reached earlier backwards
    out_seen = falses(nv(g)) # nodes reached earlier forwards
    descendant = falses(nv(g)) # descendant in s
    blocked = falses(nv(g))
#    given = falses(nv(g)) # descendant in s
    
    
    for ve in S 
        in_seen[ve] = true
        blocked[ve] = true
    end
   
    (in_seen[u] || in_seen[v]) && throw(ArgumentError("S should not contain u or v"))
    
    u == v && throw(ArgumentError("u == v"))
    

    next = Vector{T}()
    
    # mark vertices with descendants in S
    next = Vector{T}()
    for w in S
        push!(next, w)
        descendant[w] = true
    end    
    
    while !isempty(next)
        for w in in_neighbors(g, shift!(next))
            if !descendant[w]
                push!(next, w) # push onto queue
                descendant[w] = true
            end
        end
    end
 
    verbose && println(descendant)
    
    in_next = Vector{T}()
    out_next = Vector{T}()
    
    push!(in_next, u) # treat u as vertex reached backwards
    in_seen[u] = true

    while true
        sin = isempty(in_next) 
        sout = isempty(out_next) 
        sin && sout && return true # no vertices in the queue
        
   #    for i in in_next
   #         print("<-$i ")
   #     end
   #     for i in out_next
   #         print("->$i ")
   #     end
   #     println()
        
        if !sin # some vertex reach backwards in the queue
            src = shift!(in_next)  
            for w in out_neighbors(g, src) # possible collider at destination
                if !out_seen[w] && (!blocked[w] || descendant[w])
                    verbose && println("<- $src -> $w")
                    w == v && return false
                    push!(out_next, w)
                    out_seen[w] = true
                end
            end
            for w in in_neighbors(g, src)
                if !in_seen[w]
                    verbose && println("<- $src <- $w")
                    w == v && return false
                    push!(in_next, w)
                    in_seen[w] = true
                end
            end    
        end
        if !sout # some vertex reach forwards in the queue
            src = shift!(out_next)  
            for w in out_neighbors(g, src) # possible collider at destination
                if !out_seen[w] && !blocked[src] && (!blocked[w] || descendant[w])
                    verbose && println("-> $src -> $w")
                    w == v && return false
                    push!(out_next, w)
                    out_seen[w] = true
                end
            end
            for w in in_neighbors(g, src) # collider at source
                if !in_seen[w] && descendant[src] # shielded collider
                    verbose && println("-> $src <- $w")
                    w == v && return false
                    push!(out_next, w)
                    in_seen[w] = true
                end
            end    
        end

    end
    return true
end


"""
    removesorted(collection, item) -> contains(collection, item)
    
Remove item from sorted collection. 
"""
function removesorted!(n, v)
    i = searchsorted(n, v)
    isempty(i) && return false   # not found
    deleteat!(n, first(i))
    true
end

"""
    skeleton(n, I) -> g, S

Perform the undirected PC skeleton algorithm for a set of 1:n variables using the test I.
Returns skeleton graph and separating set  
"""
function skeleton(n::V, I, par...) where {V}
    g = CompleteGraph(n)
    S = Dict{edgetype(g),Vector{Vector{V}}}()
    d = 0 # depth
    while true
        isdone = true
        for e in collect(edges(g)) # cannot remove edges while iterating
            n = copy(neighbors(g, src(e)))::Vector{V}
            if length(n) > d  # i.e. |n\{dst(e)}| >= d 
                removesorted!(n, dst(e))
                isdone = false
                for s in combinations(n, d)
                    if I(src(e), dst(e), s, par...) 
                        rem_edge!(g, e)
                        if e in keys(S)
                            push!(S[e], s)
                        else
                            S[e] = [s]
                        end
                    end
                end
            end
        end
        d = d + 1
        if isdone
            return g, S
        end
    end    
end

          
function orient_unshielded(g, S)
    for e in edges(g) #loop over triples v - w - z with v < z
        v, w = src(e), dst(e)
        for z in neighbors(g, w)
            z <= v && continue 
            v in neighbors(g, z) && continue
            vz = SimpleEdge(minmax(v,z))
            [w] in S[vz] || println("$v -> $w <- $z")
            #println("$v - $w - $z")
        end
    end
end

function oracle(i, j, s, g)
    #r = !has_path(g, i, j; exclude_vertices=s)
    r = dsep(g, i, j, s)
    if r 
        println("$i ⫫ $j | $s")
    else
#        println("$i - $j | $s")
    end
    r
end        

g1 = g = DiGraph(7)
d = nv(g)
for (i,j) in [(1,2), (2,3), (2,4), (4,5), (3,5), (5,6), (7,5)]
    add_edge!(g,i,j)
end

@testset "g1" begin
@test !dsep(g, 1, 2, [])
@test !dsep(g, 2, 5, [])

@test !dsep(g, 1, 5, [3])
@test !dsep(g, 1, 5, [4])
@test dsep(g, 1, 5, [3,4])

@test dsep(g, 3, 4, [2])
@test dsep(g, 3, 4, [2, 7])

@test !dsep(g, 3, 4, [5])
@test !dsep(g, 3, 4, [2, 6])


@test !dsep(g, 3, 4, [2, 5])
@test !dsep(g, 3, 4, [2, 6])

@test !dsep(g, 3, 5, [6])
@test !dsep(g, 3, 5, [7])

end


g2 = g = DiGraph(7)
d = nv(g)
for (i,j) in [(1,3), (2,3), (3,4),(3,5), (4,6), (6, 7)]
    add_edge!(g,i,j)
end

@testset "g2" begin
@test !dsep(g, 1, 2, [4, 5])
@test dsep(g, 1, 2, [])
@test !dsep(g, 1, 2, [3])
@test dsep(g, 4, 5, [3])
@test !dsep(g, 4, 5, [])   
@test !dsep(g, 4, 5, [1,2])
end


println("Skeleton g1")
h1, S1 = skeleton(nv(g1), oracle, g1)

println("Skeleton g2")
h2, S2 = skeleton(nv(g2), oracle, g2)
#@test h == g

orient_unshielded(g1, S1)