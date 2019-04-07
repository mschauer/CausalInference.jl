[![Build Status](https://travis-ci.org/mschauer/CausalInference.jl.svg?branch=master)](https://travis-ci.org/mschauer/CausalInference.jl)
<!-- [![Coverage Status](https://coveralls.io/repos/github/mschauer/CausalInference.jl/badge.svg?branch=master)](https://coveralls.io/github/mschauer/CausalInference.jl?branch=master)
[![codecov.io](http://codecov.io/github/mschauer/CausalInference.jl/coverage.svg?branch=master)](http://codecov.io/github/mschauer/CausalInference.jl?branch=master)-->

[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://mschauer.github.io/CausalInference.jl/latest/)
	
# Causalinference.jl

Julia package for causal inference, graphical models and structure learning.

This package contains code for the stable version of the PC algorithm and the FCI algorithm as described in Zhang's article.

The algorithms use the Julia packagea [LightGraphs](https://github.com/JuliaGraphs/LightGraphs.jl) and [MetaGraphs](https://github.com/JuliaGraphs/MetaGraphs.jl). Graphs are represented by sorted adjacency lists (vectors in the implemention). CPDAGs are just `DiGraph`s where unoriented edges are represented by both a forward and a backward directed edge. PAGs are `MetaDiGraph`s where every edge is represented by a forward and a backward edge. The marks of the endpoint of an edge are stored in the `:marks` property. Marks can be checked and set using the `has_marks` and `set_marks!` functions (see documentation for details).

The PC algorithm is tested on random DAGs by comparing the result of the PC algorithm using the *d*-separation oracle with CPDAGs computed with Chickering's DAG->CPDAG conversion algorithm (implemented as `dsep` and `cpdag` in this package).

See the [documentation](https://mschauer.github.io/CausalInference.jl/latest/) for other implemented functionality and [issue #1 (Roadmap/Contribution)](https://github.com/mschauer/CausalInference.jl/issues/1) for coordination of the development.

# References

* P. Spirtes, C. Glymour, R. Scheines, R. Tillman: Automated search for causal relations: Theory and practice. *Heuristics, Probability and Causality: A Tribute to Judea Pearl* 2010
* P. Spirtes, C. Glymour, R. Scheines: Causation, Prediction, and Search. *MIT Press* 2000
* J. Zhang: On the completeness of orientation rules for causal discovery in the presence of latent confounders and selection bias. *Artificial Intelligence* 16-17 (2008), 1873-1896
* T. Richardson, P. Spirtes: Ancestral Graph Markov Models. *The Annals of Statistics* 30 (2002), 962-1030
* D. M. Chickering: Learning Equivalence Classes of Bayesian-Network Structures. *Journal of Machine Learning Research* 2 (2002), 445-498.
* D. Colombo, M. H. Maathuis: Order-Independent Constraint-Based Causal Structure Learning. *Journal of Machine Learning Research* 15 (2014), 3921-3962.


