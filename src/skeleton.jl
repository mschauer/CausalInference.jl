using LightGraphs
using LightGraphs.SimpleGraphs
using Combinatorics
using Base.Test

 
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
                        if !(e in keys(S))
                            S[e] = s
                        end
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


function oracle(i, j, s, g)
    r = dsep(g, i, j, s)
    if r 
        println("$i ⫫ $j | $s")
    else
#        println("$i - $j | $s")
    end
    r
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
