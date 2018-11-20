using CausalInference, LightGraphs, MetaGraphs

@testset "FCI utils" begin
    dg = MetaDiGraph(4)
    add_edge!(dg, 1, 2)
    set_prop!(dg, 1, 2, :mark, :arrow)
    add_edge!(dg, 3, 2)
    set_prop!(dg, 3, 2, :mark, :arrow)
    add_edge!(dg, 3, 4)
    add_edge!(dg, 4, 2)
    @test is_collider(dg, 1, 2, 3)
    @test is_triangle(dg, 2, 3, 4)
    @test !(is_triangle(dg, 1, 2, 3))
end

@testset "FCI Skeleton" begin
    N = 10000
    t1 = 2*randn(N)
    t2 = 2*randn(N)
    c = randn(N)
    b = t1 + c + randn(N)*0.25
    d = t2 + c + randn(N)*0.25
    a = t1 + d + randn(N)*0.5
    e = t2 + b + randn(N)*0.5
    df = (a=a,b=b,c=c,d=d,e=e)

    true_g = DiGraph(7)
    add_edge!(true_g,6,1)
    add_edge!(true_g,6,2)
    add_edge!(true_g,3,2)
    add_edge!(true_g,3,4)
    add_edge!(true_g,7,4)
    add_edge!(true_g,7,5)
    add_edge!(true_g,4,1)
    add_edge!(true_g,2,5)
    
    g_oracle = fcialg(5, dseporacle, true_g)
    g_gauss = fcialg(df, 0.05, gausscitest)

    # test that hidden dsep has been found
    @test !has_edge(g_oracle, 1 ,5)
    @test !has_edge(g_gauss, 1 ,5)
    @test has_edge(g_oracle, 2 ,5)
    @test has_edge(g_gauss, 2 ,5)
    @test has_edge(g_oracle, 1 ,4)
    @test has_edge(g_gauss, 1 ,4)
end
