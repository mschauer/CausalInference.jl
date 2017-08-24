using Base.Test, LightGraphs, Combinatorics


gplot(g) = graphplot(g, names=vertices(g), nodesize=2.5, fontsize=20, nodeshape=:circle)


import LightGraphs: discover_vertex!


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
    S = Dict{edgetype(g),Vector{V}}()
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
                        S[e] = s
                        break 
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

Base.start(e::LightGraphs.SimpleGraphs.SimpleEdge) = 1
Base.next(e::LightGraphs.SimpleGraphs.SimpleEdge, state) = state == 1 ? (src(e), 2) : (dst(e), 3)
Base.done(e::LightGraphs.SimpleGraphs.SimpleEdge, state) = state == 3



function orient_unshielded(g, S)
    for e in edges(g)
        v, w = src(e), dst(e)
        for z in neighbors(g, w)
            z == v && continue
            v in neighbors(g, z) && continue
            println("$v - $w - $z")
        end
    end
end


function skeleton_stable(n, I, par...)
    
    g = CompleteGraph(n)
    S = Dict{edgetype(g),Vector{Int}}()
    
    d = 0 # depth
    r_edges = Vector{Int}()
    while true
        isdone = true
        for e in edges(g) # cannot remove edges while iterating
            n = setdiff(neighbors(g, src(e)), dst(e))
            if length(n) >= d  # more than d neighbors != dst(e)
                isdone = false
                for s in combinations(n, d)
                    if I(src(e), dst(e), s, par...) 
                        append!(r_edges, e)
                        S[e] = s
                        break 
                    end
                end
            end # break
        end
        d = d + 1
        for e in r_edges
            rem_edge!(g, e)
        end
        empty!(redges)
        if isdone
            return g, S
        end
    end    
end



"""
    partialcor(i, j, s, C)

Compute conditional correlation of variables i and j given variables in s with 
C is Gaussian covariance.
"""
function partialcor(i, j, s, C)
    n = length(s)
    if n == 0
        C[i,j]
    elseif n == 1
        k = s[1]
        (C[i, j] - C[i, k]*C[j, k])/sqrt((1 - C[j, k]^2)*(1 - C[i, k]^2))
    else 
        C0 = C[[i;j;s],[i;j;s]]
        P = cholfact(C0*C0')\C0
        -P[1, 2]/sqrt(P[1, 1]*P[2, 2])
    end
    
end


function partialcor2(i, j, s, C)
    n = length(s)
    if n == 0
        C[i,j]
    elseif n == 1
        k = s[1]
        (C[i, j] - C[i, k]*C[j, k])/sqrt((1 - C[j, k]^2)*(1 - C[i, k]^2))
    else 
        C0 = C[[i;j;s],[i;j;s]]
        P = pinv(C0, 1.5e-8)
        -P[1, 2]/sqrt(P[1, 1]*P[2, 2])
    end
    
end

"""
    gausscitest(i, j, s, (C,n), c)

Test for conditional independence of variable no i and j given variables in s with 
Gaussian test at the critical value c. C is covariance of n observations.

"""
@inline function gausscitest(i, j, s, stat, c)
    C, n = stat
    r = partialcor2(i, j, s, C)
    r = clamp(r, -1, 1)
    n - length(s) - 3 <= 0 && return true # remove edges which cannot be tested for
    t = sqrt(n - length(s) - 3)*atanh(r)
    abs(t) < c
end 

truetest(i, j, s) = true
falsetest(i, j, s) = false

randtest(i, j, s, c) = rand() < c


function oracle2(i, j, s, g)
    dsep(g, i, j, s)
end        

function oracle(i, j, s, g)
    !has_path(g, i, j; exclude_vertices=s)
end        


srand(5)
let d = 30 # 40 disconnected
    n = 10000
    alpha = 0.05
    qu(x) = x*x'
    L = chol(full(qu(I + sprandn(d, d, alpha))))'
    X = L'\randn(d, n)
    Σ = cov(X,2)

    @test norm(inv(Σ) - L*L')  < .2*d*d/sqrt(n)

    C = cor(X,2)
    g = Graph((L*L' .!= 0)-I)
    println("Test Seth's version")
    @time h, s = skeleton(d, oracle, g)
    @test g == h
    println("Test old version")
    @time h, s = skeleton(d, oracle2, g)
    @test g == h
    @time h, s = skeleton(d, truetest)
    @test ne(h) == 0

    @time h, s = skeleton(d, gausscitest, (C,n), 0.8)

    O = full.(adjacency_matrix.(g)) 
    a = O +  2*full.(adjacency_matrix.(h))
    println("num edges ", div(sum(a .== 3),2), " of ", ne(g), ", fp ", div(sum(a .== 2),2) )
end
    

PLOT = false
if PLOT
    using PlotRecipes
    
    graphplot(g, node_weights=40*ones(d), names=1:d)
    graphplot(g2, node_weights=40*ones(d), names=1:d)
end


if !isfile(joinpath("data","NCI-60.csv"))
    cd("data") do
        run(`wget http://nugget.unisa.edu.au/ParallelPC/data/real/NCI-60.csv`)
    end
end 
let d = 20
    data = readcsv("data/NCI-60.csv")
    d = min(size(data,2),d)
    X = data[:,1:d]'
    d, n = size(X)
    C = Symmetric(cor(X,2))
    @time h, s = skeleton(d, gausscitest, (C,n), 2.575829) #abs(qnorm(0.01/2))
    println("inferred edges ", ne(h))
    #O = full.(adjacency_matrix.(h)) 
end


let g = Graph(5)
    d = nv(g)
    for (i,j) in [(1,2), (2,3), (2,4),(4,5), (3,5)]
        add_edge!(g,i,j)
    end
    
    h, s = skeleton(d, oracle, g)
    @test g == h
    h, s = skeleton(d, truetest)
    @test ne(h) == 0
    
end    

let d = 15
    @time h, s = skeleton(d, falsetest)
    @test ne(h) == div(d*(d-1),2)
end    

let d = 100
    @time h, s = skeleton(d, truetest)
    @test ne(h) == 0
end    

let d = 100
    @time h, s = skeleton(d, randtest, 0.01)
end    
