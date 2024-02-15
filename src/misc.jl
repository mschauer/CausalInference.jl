using Random
"""
    digraph(E)

Create `DiGraph` from edge-list.
"""
function digraph(E, d = maximum(flatten(E)))
    g = DiGraph(d)
    for (i, j) in E
        add_edge!(g, i, j)
    end
    g
end

"""
    graph(E)

Create `Graph` from edge-list.
"""
function graph(E, d = maximum(flatten(E)))  
    g = Graph(d)
    for (i, j) in E
        add_edge!(g, i, j)
    end
    g
end

"""
    arrows(g)

Return the edge-list as `Pair`s.
"""
arrows(g::SimpleDiGraph{T}) where T = nv(g) > 0 ? map(Pair, edges(g)) : Pair{T,T}[]
const vpairs = arrows

"""
    skel_oracle(g; stable=true)

Compute the `skeleton` using the `dseporacle` for the DAG `g`.
"""
skel_oracle(g; stable=true) = skeleton(nv(g), dseporacle, g; stable)

"""
    pc_oracle(g; stable=true)

Compute CPDAG using the PC algorithm using the `dseporacle` on the DAG `g`. 
"""
pc_oracle(g; stable=true) = pcalg(nv(g), dseporacle, g; stable)

"""
    randdag(n, alpha = 0.1)

Create Erdős–Rényi random DAG from randomly permuted random triangular matrix with
edge probability `alpha`.
"""
function randdag(n, alpha = 0.1)
    g = DiGraph(n)
    p = randperm(n) 
    for i in 1:n
        for j in 1:n
            i == j && continue
            rand() < alpha && add_edge!(g, p[min(i,j)], p[max(i,j)])
        end
    end
    g
end


"""
    keyedreduce(op, key::AbstractVector{T}, a, init=0.0) where T

Similar to `countmap` returning a dictionary mapping unique key in `key` to the 
reduction the given collection itr with the given binary operator `op`.

```
julia> keyedreduce(+, [:a, :b, :a], [7, 3, 2])
Dict{Symbol, Float64} with 2 entries:
  :a => 9.0
  :b => 3.0
```
"""
function keyedreduce(op, key::AbstractVector{T}, a, init=0.0) where T
    cm = Dict{T,typeof(init)}()
    for (v, τ) in zip(key, a)
        index = Base.ht_keyindex2!(cm, v)
        if index > 0
            @inbounds cm.vals[index] = op(cm.vals[index] , τ)
        else
            @inbounds Base._setindex!(cm, op(init, τ), v, -index)
        end
    end
    return cm
end