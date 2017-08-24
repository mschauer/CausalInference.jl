
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


