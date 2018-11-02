# CausalInference.jl

A Julia package for causal inference, graphical models and structure learning with the PC algorithm. This package contains for now the classical (unstable) PC algorithm [`pcalg`](@ref), tested on random DAGs by comparing the result of the PC algorithm using the *d*-separation oracle with the CPDAG computed with Chickering's DAG->CPDAG conversion algorithm (implemented as [`dsep`](@ref) and [`cpdag`](@ref) in this package).

See the [Library](https://mschauer.github.io/CausalInference.jl/latest/library.html) for other implemented functionality.

The algorithms use the `SimpleGraph` and `SimpleDiGraph` graph representation of the Julia package [LightGraphs](https://github.com/JuliaGraphs/LightGraphs.jl).
Both types of graphs are represented by sorted adjacency lists (vectors of vectors in the LightGraphs implementation).

CPDAGs are just modeled as `SimpleDiGraph`s, where unoriented edges are represented by a forward and a backward directed edge.

## Example


Compute skeleton graph `h` with separating sets `S` and CPDAG `g` from
the 47x1190 data set NCI-60 on expression profiles of miRNAs and mRNAs.

```julia
using Distributions
using CausalInference
using LightGraphs

using DelimitedFiles, LinearAlgebra, Statistics

p = 0.01

# Download data
run(`wget http://nugget.unisa.edu.au/ParallelPC/data/real/NCI-60.csv`)

# Read data and compute correlation matrix
X = readdlm("NCI-60.csv", ',')
d, n = size(X)
C = Symmetric(cor(X, dims=2))

# Compute skeleton `h` and separating sets `S`
h, S = skeleton(d, gausscitest, (C, n), quantile(Normal(), 1-p/2))

# Compute the CPDAG `g`
g = pcalg(d, gausscitest, (C, n), quantile(Normal(), 1-p/2))
```

(Using data from T. D. Le, L. Liu et al.: Inferring microRNA–mRNA causal regulatory relationships from expression data, *Bioinformatics*, vol. 29, no. 6, 765–771, 2013.)


## Performance

The speed of the algorithm is comparable with the C++ code of the R package `pcalg` after some pending optimizations.

## Contribution
See [issue #1 (Roadmap/Contribution)](https://github.com/mschauer/CausalInference.jl/issues/1) for questions and coordination of the development.

## References

* D. M. Chickering: Learning Equivalence Classes of Bayesian-Network Structures. *Journal of Machine Learning Research* 2 (2002), 445-498.
* D. Colombo, M. H. Maathuis: Order-Independent Constraint-Based Causal Structure Learning. *Journal of Machine Learning Research* 15 (2014), 3921-3962.
