using CausalInference
using Graphs
using Test

@testset "chickering" begin
    global g = randdag(200)
    @test collect(CausalInference.ordered_edges(g)) == map(Edge,  CausalInference.chickering_order(g))
end

# http://www.multimedia-computing.de/mediawiki/images/5/55/SS08_BN-Lec2-BasicProbTheory_3.pdf
E1a = [1=>2, 1=>3, 2=>4, 3=>4]
E1b = [2=>1, 1=>3, 2=>4, 3=>4]
E1c = [1=>2, 3=>1, 2=>4, 3=>4]
C1 = [1=>2, 2=>1, 1=>3, 3=>1, 2=>4, 3=>4]

# Corrected from 
E2a = [1=>2, 2=>3, 3=>5, 1=>4, 4=>5]
C2 = [1=>2, 2=>1, 2=>3, 3=>2, 3=>5, 1=>4, 4=>1, 4=>5]

E3 = [6=>4, 4=>3, 3=>2, 3=>1, 2=>1, 3=>5]
C3 = copy(E3)
append!(C3, reverse.(E3))

h, S = skel_oracle(digraph(E1a))
@test vpairs(h) == E1a
@test length(S) == 2
@test S[Edge(2, 3)] == [1]
@test S[Edge(1, 4)] == [2, 3]

h = pc_oracle(digraph(E1a))
@test sort(C1) == sort(vpairs(h))
@test digraph(C1) == cpdag(digraph(E1a))
@test digraph(C1) == cpdag(digraph(E1b))
@test digraph(C1) == cpdag(digraph(E1c))


for (C, Es) in [(C1, [E1a, E1b, E1c]),
                (C2, [E2a]),
                (C3, [E3]),
    ]
    for E in Es
        global h = pc_oracle(digraph(E)) 
        @test sort(C) == sort(vpairs(h))
        @test digraph(C) == cpdag(digraph(E))        
    end
end

for n in 0:10
    
    alpha = 0.1
    @testset "randdag($n)" begin
        for k = 1:1000
            global g = randdag(n, alpha)
            h1 = pc_oracle(g) 
            h2 = cpdag(g)    
            h1 == h2  || println(vpairs(g))
            @test vpairs(h1) ⊆ vpairs(h2) 
            @test vpairs(h2) ⊆ vpairs(h1) 
            
            #@test h1 == h2     
        end
    end
end

@testset "randdag(50)" begin # up to level 10 or so
    n = 50
    global g = randdag(n, 0.07)
    h1 = pc_oracle(g) 
    h2 = cpdag(g)    
    h1 == h2  || println(vpairs(g))
    @test vpairs(h1) ⊆ vpairs(h2) 
    @test vpairs(h2) ⊆ vpairs(h1) 
end

@testset "orient after v structre" begin
    global g = digraph([1=>4, 2=>1, 3=>1])
    @test sort(vpairs(vskel(g))) ==  [1=>4, 2=>1, 3=>1, 4=>1]

    h1 = pc_oracle(g) 
    h2 = cpdag(g)    
    @test h1 == h2 
end

@testset "check that edge is still undirected" begin
    global g = digraph([1=>4, 2=>1, 3=>1, 3=>4])
    @test sort(vpairs(vskel(g))) == [1=>4, 2=>1, 3=>1, 3=>4, 4=>1, 4=>3]

    h1 = pc_oracle(g) 
    h2 = cpdag(g)    
    @test h1 == h2 
end
