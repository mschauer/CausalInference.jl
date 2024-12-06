using Graphs
using Graphs.SimpleGraphs
using Combinatorics
using LinearAlgebra
using Tables

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

#=
function remove_edges!(g::SimpleDiGraph, es)
    bs = falses(nv(g))
    fs = falses(nv(g))
    
    for e in es
        u, v = Pair(e)
        push!(g.fadjlist[u], v)
        push!(g.badjlist[v], u)
        fs[u] = bs[v] = true 
    end      
end
=#

struct IClosure{F<:Function,S,T}
    I::F
    args::S
    kwargs::T
end
IClosure(I, args) = IClosure(I, args, (;))
(I::IClosure)(x, y, s) = I.I(x, y, s, I.args...; I.kwargs...)

"""
    skeleton(n::Integer, I) -> g, S
    skeleton(g, I) -> g, S

Perform the undirected PC skeleton algorithm for a set of `1:n` variables using the test `I`.
Start with a subgraph `g` or the complete undirected graph on `n` vertices.
Returns skeleton graph and separating set.
"""
function skeleton(n::Integer, I, par...; stable=true)
    if stable
        skeleton_stable(complete_graph(n), IClosure(I, par))
    else 
        skeleton(complete_graph(n), IClosure(I, par))
    end
end
function skeleton(g::SimpleGraph, I)
    V = eltype(g)
    n = nv(g)
    S = Dict{edgetype(g), Vector{V}}()
    d = 0 # depth

    while true
        isdone = true
        for e0 in collect(edges(g)) # cannot remove edges while iterating
            for e in (e0, reverse(e0))
                nb = neighbors(g, src(e))
                if length(nb) > d  # i.e. |nb\{dst(e)}| >= d
                    r = searchsorted(nb, dst(e)) # neighbors are sorted in SimpleGraph
                    isempty(r) && continue # the edge does not exist anymore
                    isdone = false
                    for s in combinations_without(nb, d, first(r)) # do not modify s!
                        if I(src(e), dst(e), s)
                            @debug "Removing edge $(e0) given $(s)"
                            rem_edge!(g, e0)
                            if !(e0 in keys(S))
                                S[e0] = s
                            end
                            @goto nextedge
                        end
                    end
                end
            end
            @label nextedge
        end
        d = d + 1
        if isdone || d > n-2 
            return g, S
        end
    end
end

"""
    skeleton_stable(n::Integer, I) -> g, S
    skeleton_stable(g, I) -> g, S

Perform the undirected stable PC skeleton algorithm for a set of `1:n` variables using the test `I`.
Start with a subgraph `g` or the complete undirected graph on `n` vertices.
Returns skeleton graph and separating set.
"""
function skeleton_stable(n::Integer, I, par...)
    skeleton_stable(complete_graph(n), IClosure(I, par))
end
function skeleton_stable(g, I)
    l = ReentrantLock()
    V = eltype(g)
    n = nv(g)
    S = Dict{edgetype(g), Vector{V}}()
    d = 0 # depth
    remove = Vector{edgetype(g)}()
    while true
        isdone = true
        Threads.@threads for e0 in collect(edges(g)) # cannot remove edges while iterating
            for e in (e0, reverse(e0))
                nb = neighbors(g, src(e))
                if length(nb) > d  # i.e. |nb\{dst(e)}| >= d
                    r = searchsorted(nb, dst(e))
                    isempty(r) && continue # the edge does not exist anymore
                    isdone = false
                    for s in combinations_without(nb, d, first(r)) # do not modify s!                
                        if I(src(e), dst(e), s)
                            lock(l) do                            
                                push!(remove, e)
                                if !(e0 in keys(S))
                                    S[e0] = s
                                end
                            end
                            @goto nextedge
                        end
                    end
                end
            end
            @label nextedge
        end
        for r in remove
            rem_edge!(g, r) # could be efficienter
        end
        empty!(remove)
        d = d + 1
        if isdone || d > n-2 
            return g, S
        end
    end
end

"""
    dseporacle(i, j, s, g)

Oracle for the `skeleton` and `pcalg` functions using [`dsep`](@ref) on the true graph `g`
"""
function dseporacle(i, j, s, g; sel = [])
    dsep(g, i, j, vcat(s, sel))
end

"""
    partialcor(i, j, s, C)

Compute the partial correlation of nodes `i` and `j` given list of nodes `s`
using the correlation matrix `C`.
"""
function partialcor(i, j, s, C)
    n = length(s)
    if n == 0
        C[i, j]
    elseif n == 1
        k = s[1]
        (C[i, j] - C[i, k] * C[j, k]) / sqrt(max(eps(), (1 - C[j, k]^2) * (1 - C[i, k]^2)))
    else
        is = zeros(Int, n + 2)
        is[1] = i
        is[2] = j
        for k in 1:n
            is[k + 2] = s[k]
        end
        C0 = C[is, is] # @view not supported on older versions
        P = pinv(C0, rtol=1.5e-8)
        -P[1, 2] / sqrt(P[1, 1] * P[2, 2])
    end
end

"""
    gausscitest(i, j, s, (C,n), c)

Test for conditional independence of variable no i and j given variables in s with
Gaussian test at the critical value c. C is covariance of n observations.

"""
@inline function gausscitest(i, j, s, stat, c)
    C, n = stat
    r = partialcor(i, j, s, C)
    r = clamp(r, -1, 1)
    n - length(s) - 3 <= 0 && return true # remove edges which cannot be tested for
    t = sqrt(n - length(s) - 3) * atanh(r)
    #@debug "testing $(i)-$(j) given $(s): $(abs(t)) -- $(c)"
    abs(t) < c
end

"""
    cmitest(i,j,s,data,crit; kwargs...)

Test for conditional independence of variables i and j given variables in s with
permutation test using nearest neighbor conditional mutual information estimates
at p-value crit.

keyword arguments:
kwargs...: keyword arguments passed to independence tests
"""
function cmitest(i, j, s, data, crit; kwargs...)
    columns = Tables.columns(data)

    x = collect(transpose(convert(Array, Tables.getcolumn(columns, i))))
    y = collect(transpose(convert(Array, Tables.getcolumn(columns, j))))

    if length(s) == 0
        res = kl_perm_mi_test(x, y; kwargs...)
    else
        z = reduce(vcat,
                   map(c -> collect(transpose(convert(Array, Tables.getcolumn(columns, c)))),
                       s))
        res = kl_perm_cond_mi_test(x, y, z; kwargs...)
    end

    #@debug "CMI test for $(i)-$(j) given $(s): $(res) compared to $(crit)"
    return res > crit
end

truetest(i, j, s) = true
falsetest(i, j, s) = false
