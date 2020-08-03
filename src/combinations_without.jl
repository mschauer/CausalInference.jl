# Combinations iterator with tabu
using Combinatorics
import Base: iterate, eltype, length

struct CombinationsWithout{T}
     a::T
     t::Int
     w::Int
     l::Int
     L::Int
 end

 function Base.iterate(c::CombinationsWithout)
     s = zeros(Int, c.t)
     if 0 < c.t <= c.l
         s[1] = 1
         if c.w == 1
             s[1] = 2
         end
         i = 1
         while i < c.t
             i += 1
             s[i] = s[i-1] + 1
             if s[i] == c.w
                 s[i] += 1
             end
         end
     end
     Base.iterate(c, (s, Vector{eltype(c.a)}(undef, c.t), 1))
 end

 function Base.iterate(c::CombinationsWithout, state)
     s, comb, k =  state
     k > c.L && return nothing
     k += 1
     for i in 1:c.t
         comb[i] = c.a[s[i]]
     end
     for i in length(s):-1:1
         s[i] += 1
         if s[i] > c.l - (length(s)-i) - (s[i] <= c.w)
             continue
         end
         if s[i] == c.w
             s[i] += 1
         end
         for j = i+1:length(s)
             s[j] = s[j-1]+1
             if s[j] == c.w
                 s[j] += 1
             end
         end
         break
     end
     comb, (s, comb, k)
 end

 length(c::CombinationsWithout) = c.L
 eltype(::Type{CombinationsWithout{T}}) where {T} = Vector{eltype(T)}

 """
     combinations_without(a, n::Integer, w)

 Generate all combinations of `n` elements from an indexable object except those with index `w`. Note that the combinations
 are modified inplace.
 """
 function combinations_without(a, t::Integer, w)
     l = length(a)
     ins = 1 <= w <= l
     L = binomial(l - ins, t)
     if ins && length(a) == w # w is last index
         CombinationsWithout(a, t, 0, l - 1, L)
     elseif ins && l == t # t too big
         CombinationsWithout(a, t, w, l, 0)
     elseif ins
         CombinationsWithout(a, t, w, l, L)
     else
         CombinationsWithout(a, t, 0, l, L)
     end
 end
