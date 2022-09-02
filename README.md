
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://mschauer.github.io/CausalInference.jl/latest/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1005091.svg)](https://doi.org/10.5281/zenodo.1005091)

# Causalinference.jl

Julia package for causal inference, graphical models and structure learning.

This package contains code for the stable version of the PC algorithm and the extended FCI algorithm.

See the [documentation](https://mschauer.github.io/CausalInference.jl/latest/) for implemented functionality and [issue #1 (Roadmap/Contribution)](https://github.com/mschauer/CausalInference.jl/issues/1) for coordination of the development.

![Example output of PC algorithm](assets/pc_graph_linear.png)


# FAQ

**Q:** I looked for "causal inference" and found [CausalInference.jl](.) and [Omega.jl](http://www.zenna.org/Omega.jl/latest/causal/)... 
**A:** CausalInference.jl is about causal discovery; you observe the data and want to infer the causal structure. Omega lets you reason what happens then: when you intervene ("do calculus") and want to cause changes.

