using Graphs
using CausalInference
using Test

import CausalInference.tails_and_adj_neighbors
import CausalInference.adj_neighbors
using CausalInference: precompute_semidirected, isadjacent, neighbors_undirected,
        isclique, InsertIterator, sorted_intersect_, DeleteIterator, next_CPDAG, →, has_both


isblocked(g, x, y, nodesRemoved) = !has_a_path(g, [x], y, nodesRemoved)
function tails_and_adj_neighbors(g, x, y)
    Nb = neighbors_undirected(g, y)
    a = Bool[isadjacent(g, t, x) for t in Nb]
    Nb[.~ a], Nb[a]
end 
function adj_neighbors(g, x, y)
#    a = intersect(inneighbors(g,y), outneighbors(g,y), all_neighbors(g,x))
    sorted_intersect_(neighbors_undirected(g,y), all_neighbors(g,x))
end 

"""
    Insert!(g, x, y, T)

Inserts x->y and directs previously undirected edges t->y, t ∈ T.
Here x and y are not adjacent and T are undirected-neighbors of y 
that are not adjacent to x.
"""
function Insert!(g, x, y, T)
    add_edge!(g, x, y)
    
    # Orient all edges in T incident into child node
    for t ∈ T
        orientedge!(g, t, y)
    end
    return g
end

"""
    Delete!(g, x, y, H)

Deletes x-y or x->y and directs previously undirected edges x->h and y->h
for h in H.
"""
function Delete!(g, x, y, H)

    # Remove directed or undirected edges (x→y and x-y)
    rem_edge!(g, x → y)
    rem_edge!(g, y → x)
    
    # Orient all vertices in H toward x and y
    for h ∈ H
        if has_both(g, x, h) # reading literally Definition 13, Chickering
            orientedge!(g, x → h) 
        end
        orientedge!(g, y → h)
    end

    return nothing
end

@testset "operators randgraphs" begin
    for rep = 1:10 
        for alpha in [0.01, 0.05, 0.1]
            d = randdag(50, 0.05)
            g = cpdag(d)
            for y in vertices(g)
                semidirected = precompute_semidirected(g, y)
                for x in vertices(g) 
                    insertoperators = Set{Set{Int64}}()
                    insertoperatorsit = Set{Set{Int64}}()
                    deleteoperators = Set{Set{Int64}}()
                    deleteoperatorsit = Set{Set{Int64}}()
                    if !isadjacent(g, x, y) 
                        Tyx, NAyx = tails_and_adj_neighbors(g, x, y)
                        for i in 0:UInt128(2)^length(Tyx) - 1
                            T = Tyx[digits(Bool, i, base=2, pad=length(Tyx))]
                            NAyxT = CausalInference.sorted_union_(NAyx, T)
                            valid = (isclique(g, NAyxT) && isblocked(g, y, x, NAyxT))
                            if valid
                                push!(insertoperators, Set(T)) 
                            end
                        end 
                        for T in InsertIterator(g, x, y, semidirected)
                            push!(insertoperatorsit, Set(T))
                        end
                    end
                    @test insertoperators == insertoperatorsit
                    if has_edge(g, x, y)
                        Hyx = adj_neighbors(g, x, y)
                        for i in 0:UInt128(2)^length(Hyx) - 1
                            mask = digits(Bool, i, base=2, pad=length(Hyx))
                            H = Hyx[mask] 
                            NAyx_H = Hyx[map(~, mask)] 
                            valid = isclique(g, NAyx_H)
                            if valid
                                push!(deleteoperators, Set(H)) 
                            end
                        end
                        for T in DeleteIterator(g, x, y)
                            push!(deleteoperatorsit, Set(T))
                        end
                    end
                    @test deleteoperators == deleteoperatorsit
                end
            end
        end
    end
end

@testset "nextCPDAG randgraphs" begin  
    for rep = 1:10 
        for alpha in [0.01, 0.05, 0.1]
            d = randdag(50, 0.05)
            g = cpdag(d)
            x = -1
            y = -1
            while true 
                x = rand(vertices(g))
                y = rand(vertices(g))
                semidirected = precompute_semidirected(g, y)
                if x != y && !isadjacent(g, x, y) && !last(semidirected)[x]
                    break
                end
            end
            semidirected = precompute_semidirected(g, y)
            T = rand(collect(InsertIterator(g, x, y, semidirected)))
            wien = next_CPDAG(g, :up, x, y, T)
            meek = copy(g)
            Insert!(meek, x, y, T)
            vskel!(meek)
            meek_rules!(meek)
            @test wien == meek
            x = -1
            y = -1
            while x == y || !has_edge(g, x, y)
                x = rand(vertices(g))
                y = rand(vertices(g))
            end
            H = rand(collect(DeleteIterator(g, x, y)))
            wien = next_CPDAG(g, :down, x, y, H)
            meek = copy(g)
            Delete!(meek, x, y, H)
            vskel!(meek)
            meek_rules!(meek)
            @test wien == meek
        end
    end
end
