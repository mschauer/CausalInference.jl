using CausalInference, Graphs, MetaGraphs, Random

@testset "FCI utils" begin
    dg = MetaDiGraph(4)
    add_edge!(dg, 1, 2)
    add_edge!(dg, 2, 1)
    set_prop!(dg, 1, 2, :mark, :arrow)
    set_prop!(dg, 2, 1, :mark, :tail)

    @test has_marks(dg, 1, 2, arrow"-->")
    @test has_marks(dg, 2, 1, arrow"*--")

    add_edge!(dg, 3, 2)
    add_edge!(dg, 2, 3)
    set_marks!(dg, 2, 3, arrow"<-*")
    add_edge!(dg, 3, 4)
    add_edge!(dg, 4, 2)

    @test is_collider(dg, 1, 2, 3)
    @test is_triangle(dg, 2, 3, 4)
    @test !(is_triangle(dg, 1, 2, 3))
    @test is_parent(dg, 1, 2)

    dg = MetaDiGraph(5)
    add_edge!(dg, 4, 5)

    add_edge!(dg, 3, 4)
    add_edge!(dg, 4, 3)
    set_marks!(dg, 3, 4, arrow"<->")

    add_edge!(dg, 3, 5)
    add_edge!(dg, 5, 3)
    set_marks!(dg, 3, 5, arrow"-->")
    
    add_edge!(dg, 2, 3)
    add_edge!(dg, 3, 2)
    set_marks!(dg, 2, 3, arrow"<->")
    
    add_edge!(dg, 2, 5)
    add_edge!(dg, 5, 2)
    set_marks!(dg, 2, 5, arrow"-->")
    
    add_edge!(dg, 1, 2)
    set_marks!(dg, 1, 2, arrow"*->")

    @test is_discriminating_path(dg, collect(1:5))
end

@testset "FCI Skeleton" begin
    Random.seed!(123)

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

    # test that edges remain
    @test has_edge(g_oracle, 2 ,5)
    @test has_edge(g_gauss, 2 ,5)
    @test has_edge(g_oracle, 1 ,4)
    @test has_edge(g_gauss, 1 ,4)
end

@testset "FCI Orientation" begin
    true_g = DiGraph(4)
    add_edge!(true_g,1,3)
    add_edge!(true_g,2,3)
    add_edge!(true_g,3,4)
    g_oracle = fcialg(4, dseporacle, true_g, verbose=true)

    @test has_marks(g_oracle, 1, 3, arrow"o->")
    @test has_marks(g_oracle, 3, 4, arrow"-->")

    true_g = DiGraph(5)
    add_edge!(true_g,1,2)
    add_edge!(true_g,3,4)
    add_edge!(true_g,5,2)
    add_edge!(true_g,5,4)
    g_oracle = fcialg(4, dseporacle, true_g, verbose=true)

    @test has_marks(g_oracle, 2, 4, arrow"<->")
    @test has_marks(g_oracle, 1, 2, arrow"o->")

    # test graph from Figure 11-iv in Richardson & Sprites, 2002
    true_g = DiGraph(7)
    add_edge!(true_g, 1, 5)
    add_edge!(true_g, 5, 2)
    add_edge!(true_g, 6, 2)
    add_edge!(true_g, 2, 7)
    add_edge!(true_g, 7, 3)
    add_edge!(true_g, 3, 4)
    add_edge!(true_g, 6, 4)
    g_oracle = fcialg(4, dseporacle, true_g, sel=[7], verbose=true)

    # test for that weird edge I don't understand...
    @test has_marks(g_oracle, 1, 4, arrow"o->")
    
    # test graph Figure 6 in Zhang, 2008
    true_g = DiGraph(5)
    add_edge!(true_g, 1, 2)
    add_edge!(true_g, 2, 3)
    add_edge!(true_g, 4, 3)
    add_edge!(true_g, 5, 1)
    add_edge!(true_g, 5, 4)
    g_oracle = fcialg(4, dseporacle, true_g, verbose=true)

    @test has_marks(g_oracle, 1, 4, arrow"o-o")
    @test has_marks(g_oracle, 4, 3, arrow"-->")
    @test has_marks(g_oracle, 2, 3, arrow"-->")
    @test has_marks(g_oracle, 1, 2, arrow"o-o")
end
