# Solving the problem in https://bsky.app/profile/p-hunermund.com/post/3lci6xojlmt25

using CausalInference, Graphs, Test
@testset "z-intervention" begin
    # Define the causal Dag?
    U, V, X, Y, Z = 1:5
    dag = digraph([U=>X, U=>Z, V=>Y, V=>Z, Z=>X, X=>Y])
    
    ∅ = Set{Int}()
    observed = [X, Y, Z]
    
    # Can we estimate the causal effect from purely observational data, perhaps by a covariate adjustement?
    adjustments = collect(list_covariate_adjustment(dag, X, Y, ∅, observed))
    @test isempty(adjustments)
    
    # We are not out of luck, because we can intervene on Z.
    dag2 = copy(dag)
    CausalInference.do!(dag2, Z)
        
    # In the DAG with the intervation, all incoming edges to Z=5 are removed, e.g. 1=>5.
    # Enumerate again the possible adjustment?
    adjustments2 = collect(list_covariate_adjustment(dag2, X, Y, ∅, observed))

    @test Set(adjustments2) == Set([∅, Set([Z])])

end