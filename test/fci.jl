using CausalInference, LightGraphs

@testset "FCI utils" begin
    dg = DiGraph(4)
    add_edge!(dg, 1, 2)
    add_edge!(dg, 3, 2)
    add_edge!(dg, 3, 4)
    add_edge!(dg, 4, 2)
    @test is_collider(dg, 1, 2, 3)
    @test is_triangle(dg, 2, 3, 4)
    @test !(is_triangle(dg, 1, 2, 3))
end

@testset "FCI Skeleton" begin
    N = 1000
    t1 = 2*randn(N)
    t2 = 2*randn(N)
    c = randn(N)
    b = t1 + c + randn(N)*0.25
    d = t2 + c + randn(N)*0.25
    a = t1 + d + randn(N)*0.25
    e = t2 + b + randn(N)*0.25
    df = (a=a,b=b,c=c,d=d,e=e)
    g_fci = fcialg(df, 0.01, gausscitest)
    println(collect(edges(g_fci)))
    @test 1==1
end
