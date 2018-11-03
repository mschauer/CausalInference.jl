
[![Build Status](https://travis-ci.org/mschauer/CausalInference.jl.svg?branch=master)](https://travis-ci.org/mschauer/CausalInference.jl)
<!-- [![Coverage Status](https://coveralls.io/repos/github/mschauer/CausalInference.jl/badge.svg?branch=master)](https://coveralls.io/github/mschauer/CausalInference.jl?branch=master)
[![codecov.io](http://codecov.io/github/mschauer/CausalInference.jl/coverage.svg?branch=master)](http://codecov.io/github/mschauer/CausalInference.jl?branch=master)-->

[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://mschauer.github.io/CausalInference.jl/latest/)

# CausalInference.jl

Julia package for causal inference, graphical models and structure learning with the PC algorithm. This package contains for now the classical PC algorithm, tested on random DAGs by comparing the result of the PC algorithm using the *d*-separation oracle with CPDAGs computed with Chickering's DAG->CPDAG conversion algorithm (implemented as `dsep` and `cpdag` in this package).

See the [documentation](https://mschauer.github.io/CausalInference.jl/latest/) for other implemented functionality and [issue #1 (Roadmap/Contribution)](https://github.com/mschauer/CausalInference.jl/issues/1) for coordination of the development.

The algorithms use the Julia package [LightGraphs](https://github.com/JuliaGraphs/LightGraphs.jl). Graphs are represented by sorted adjacency lists (vectors in the implemention). CPDAGs are just `DiGraph`s where unoriented edges are represented by both a forward and a backward directed edge.

# References

* D. M. Chickering: Learning Equivalence Classes of Bayesian-Network Structures. *Journal of Machine Learning Research* 2 (2002), 445-498.
* D. Colombo, M. H. Maathuis: Order-Independent Constraint-Based Causal Structure Learning. *Journal of Machine Learning Research* 15 (2014), 3921-3962.


