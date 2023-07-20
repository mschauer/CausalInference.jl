# CausalInference.jl

A Julia package for causal inference, graphical models and structure learning. This package contains the stable PC algorithm [`pcalg`](@ref) and the extended FCI algorithm, as well as functionality to compute adjustment sets which can be used as control variables to compute causal effects by regression. 

## Introduction

The first aim of this package is to provide Julia implementations of popular algorithms for causal structure identification, the PC algorithm and the FCI algorithm. The aim of these algorithms is to identify causal relationships in observational data alone, in circumstances where running experiments or A/B tests is impractical or even impossible. While identification of all causal relationships in observational data is not always possible, both algorithms clearly indicate which causal can and which cannot be determined from observational data. Secondly, similarly to DAGitty, covariate adjustment sets for estimating causal effects for example by regression can be computed. 

Causal inference is by no means an easy subject. Readers without any prior exposure to these topics are encouraged to go over the following resources in order to get a basic idea of what's involved in causal inference:

- [Causal Inference in Statistics: A Primer](https://www.wiley.com/en-gb/Causal+Inference+in+Statistics%3A+A+Primer-p-9781119186847)
- [Review of Causal Discovery Methods Based on Graphical Models](https://www.frontiersin.org/articles/10.3389/fgene.2019.00524/full)
- [On Pearl’s Hierarchy and the Foundations of Causal Inference](https://causalai.net/r60.pdf)

There are also tutorials and examples linked in the navigation bar of this package.

## Implementation Details

The PC algorithm was tested on random DAGs by comparing the result of the PC algorithm using the *d*-separation oracle with the CPDAG computed with Chickering's DAG->CPDAG conversion algorithm (implemented as [`dsep`](@ref) and [`cpdag`](@ref) in this package).

See the [Library](https://mschauer.github.io/CausalInference.jl/latest/library/) for other implemented functionality.

The algorithms use the `SimpleGraph` and `SimpleDiGraph` graph representation of the Julia package [Graphs](https://github.com/JuliaGraphs/Graphs.jl).
Both types of graphs are represented by sorted adjacency lists (vectors of vectors in the Graphs implementation).

CPDAGs are just modeled as `SimpleDiGraph`s, where unoriented edges are represented by a forward and a backward directed edge.

The listing algorithms for adjustment sets are implemented from scratch using an memority efficient iterator protocol to handle large problems.

## Performance
The speed of the PC algorithm is comparable with the C++ code of the R package `pcalg`.

## Plotting
Main package provides a text-based output describing all identified edges for PC and FCI algorithm ([`plot_pc_graph_text`](@ref) and [`plot_fci_graph_text`](@ref), respectively).

In addition, additional plotting backends are supported with lazy code loading orchestrated by [Requires.jl](https://github.com/JuliaPackaging/Requires.jl). Upon importing of [TikzGraphs.jl](https://github.com/JuliaTeX/TikzGraphs.jl), additional plotting methods [`plot_pc_graph_tikz`](@ref) and [`plot_fci_graph_tikz`](@ref) will be loaded (these are also aliased as [`plot_pc_graph`](@ref) and [`plot_fci_graph`](@ref) for backward compatibility). Similarly, upon importing of both [GraphRecipes.jl](https://github.com/JuliaPlots/GraphRecipes.jl) and [Plots.jl](https://github.com/JuliaPlots/Plots.jl), additional plotting methods [`plot_pc_graph_recipes`](@ref) and [`plot_fci_graph_recipes`](@ref) will be loaded.

At the time of writing (December 2022), TikzGraphs.jl cannot be installed on ARM-based systems, so GraphRecipes.jl + Plots.jl is the recommended plotting backend in such cases.

## Contribution
See [issue #1 (Roadmap/Contribution)](https://github.com/mschauer/CausalInference.jl/issues/1) for questions and coordination of the development.

## References

* P. Spirtes, C. Glymour, R. Scheines, R. Tillman: Automated search for causal relations: Theory and practice. *Heuristics, Probability and Causality: A Tribute to Judea Pearl* 2010
* P. Spirtes, C. Glymour, R. Scheines: Causation, Prediction, and Search. *MIT Press* 2000
* J. Zhang: On the completeness of orientation rules for causal discovery in the presence of latent confounders and selection bias. *Artificial Intelligence* 16-17 (2008), 1873-1896
* T. Richardson, P. Spirtes: Ancestral Graph Markov Models. *The Annals of Statistics* 30 (2002), 962-1030
* D. M. Chickering: Learning Equivalence Classes of Bayesian-Network Structures. *Journal of Machine Learning Research* 2 (2002), 445-498.
* D. Colombo, M. H. Maathuis: Order-Independent Constraint-Based Causal Structure Learning. *Journal of Machine Learning Research* 15 (2014), 3921-3962.
* B. van der Zander, M. Liśkiewicz, J. Textor: Separators and Adjustment Sets in Causal Graphs: Complete Criteria and an Algorithmic Framework. (https://doi.org/10.48550/arXiv.1803.00116)

