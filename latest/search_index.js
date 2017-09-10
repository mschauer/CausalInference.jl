var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#CausalInference.jl-1",
    "page": "Home",
    "title": "CausalInference.jl",
    "category": "section",
    "text": "A Julia package for causal inference, graphical models and structure learning with the PC algorithm. This package contains for now the classical (unstable) PC algorithm pcalg, tested on random DAGs by comparing the result of the PC algorithm using the d-separation oracle with the CPDAG computed with Chickering's DAG->CPDAG conversion algorithm (implemented as dsep and cpdag in this package).See the Library for other implemented functionality.The algorithms use the SimpleGraph and SimpleDiGraph graph representation of the Julia package LightGraphs. Both types of graphs are represented by sorted adjacency lists (vectors of vectors in the LightGraphs implemention).CPDAGs are just modelled as SimpleDiGraphs, where unoriented edges are represented by a forward and a backward directed edge."
},

{
    "location": "index.html#Performance-1",
    "page": "Home",
    "title": "Performance",
    "category": "section",
    "text": "The speed of the algorithm is comparable with the C++ code of the R package pcalg after some pending optimizations. "
},

{
    "location": "index.html#Contribution-1",
    "page": "Home",
    "title": "Contribution",
    "category": "section",
    "text": "See issue #1 (Roadmap/Contribution) for questions and coordination of the development."
},

{
    "location": "index.html#References-1",
    "page": "Home",
    "title": "References",
    "category": "section",
    "text": "D. M. Chickering: Learning Equivalence Classes of Bayesian-Network Structures. Journal of Machine Learning Research 2 (2002), 445-498.\nD. Colombo, M. H. Maathuis: Order-Independent Constraint-Based Causal Structure Learning. Journal of Machine Learning Research 15 (2014), 3921-3962."
},

{
    "location": "library.html#",
    "page": "Library",
    "title": "Library",
    "category": "page",
    "text": ""
},

{
    "location": "library.html#Library-1",
    "page": "Library",
    "title": "Library",
    "category": "section",
    "text": ""
},

{
    "location": "library.html#CausalInference.dsep",
    "page": "Library",
    "title": "CausalInference.dsep",
    "category": "Function",
    "text": "dsep(g::AbstractGraph, u, v, s; verbose = false)\n\nCheck  whether u and v are d-separated given set s. Algorithm: unrolled https://arxiv.org/abs/1304.1505\n\n\n\n"
},

{
    "location": "library.html#CausalInference.cpdag",
    "page": "Library",
    "title": "CausalInference.cpdag",
    "category": "Function",
    "text": "cpdag(skel::DiGraph)\n\nReference: M. Chickering: Learning equivalence classes of Bayesian network structures. Journal of Machine Learning Research 2 (2002). M. Chickering: A Transformational Characterization of Equivalent Bayesian Network Structures. (1995).\n\nNote that the edge order defined there is already partly encoded into the representation of a DiGraph.\n\n\n\n"
},

{
    "location": "library.html#CausalInference.vskel",
    "page": "Library",
    "title": "CausalInference.vskel",
    "category": "Function",
    "text": "vskel(g)\n\nSkeleton and v-structures. (Currently from the first step of the pc-Alg.)\n\n\n\n"
},

{
    "location": "library.html#Directed-acyclic-graphs-(DAGs)-1",
    "page": "Library",
    "title": "Directed acyclic graphs (DAGs)",
    "category": "section",
    "text": "dsep\ncpdag\nvskel"
},

{
    "location": "library.html#CausalInference.pcalg",
    "page": "Library",
    "title": "CausalInference.pcalg",
    "category": "Function",
    "text": "pcalg(n::V, I, par...)\n\nPerform the PC skeleton algorithm for a set of 1:n variables using the tests\n\nI(u, v, [s1, ..., sn], par...)\n\nReturns the CPDAG as DiGraph.   \n\n\n\n"
},

{
    "location": "library.html#CausalInference.skeleton",
    "page": "Library",
    "title": "CausalInference.skeleton",
    "category": "Function",
    "text": "skeleton(n, I) -> g, S\n\nPerform the undirected PC skeleton algorithm for a set of 1:n variables using the test I. Returns skeleton graph and separating set  \n\n\n\n"
},

{
    "location": "library.html#CausalInference.dseporacle",
    "page": "Library",
    "title": "CausalInference.dseporacle",
    "category": "Function",
    "text": "dseporacle(i, j, s, g)\n\nOracle for the skeleton and pcalg functions using dsep on the true graph g     \n\n\n\n"
},

{
    "location": "library.html#CausalInference.unshielded",
    "page": "Library",
    "title": "CausalInference.unshielded",
    "category": "Function",
    "text": "unshielded(g, S)\n\nFind unshielded triples in the skeleton. Triples are connected vertices v-w-z where z is not a neighbour of v. Uses that edges iterates in lexicographical order.\n\n\n\n"
},

{
    "location": "library.html#PC-algorithm-1",
    "page": "Library",
    "title": "PC algorithm",
    "category": "section",
    "text": "pcalg\nskeleton\ndseporacle\nunshielded"
},

{
    "location": "library.html#CausalInference.gausscitest",
    "page": "Library",
    "title": "CausalInference.gausscitest",
    "category": "Function",
    "text": "gausscitest(i, j, s, (C,n), c)\n\nTest for conditional independence of variable no i and j given variables in s with  Gaussian test at the critical value c. C is covariance of n observations.\n\n\n\n"
},

{
    "location": "library.html#CausalInference.partialcor",
    "page": "Library",
    "title": "CausalInference.partialcor",
    "category": "Function",
    "text": "partialcor(i, j, s, C)\n\nCompute the partial correlation of nodes i and j given list of nodes s using the correlation matrix C.\n\n\n\n"
},

{
    "location": "library.html#Statistics-1",
    "page": "Library",
    "title": "Statistics",
    "category": "section",
    "text": "gausscitest\npartialcor"
},

{
    "location": "library.html#CausalInference.digraph",
    "page": "Library",
    "title": "CausalInference.digraph",
    "category": "Function",
    "text": "digraph(E)\n\nCreate DiGraph from edge-list.\n\n\n\n"
},

{
    "location": "library.html#CausalInference.pairs",
    "page": "Library",
    "title": "CausalInference.pairs",
    "category": "Function",
    "text": "pairs(g)\n\nReturn the edge-list as Pairs.\n\n\n\n"
},

{
    "location": "library.html#CausalInference.pc_oracle",
    "page": "Library",
    "title": "CausalInference.pc_oracle",
    "category": "Function",
    "text": "pc_oracle(g)\n\nCompute CPDAG using the PC algorithm using the dseporacle on the DAG g. \n\n\n\n"
},

{
    "location": "library.html#CausalInference.skel_oracle",
    "page": "Library",
    "title": "CausalInference.skel_oracle",
    "category": "Function",
    "text": "skel_oracle(g)\n\nCompute the skeleton using the dseporacle for the DAG g.\n\n\n\n"
},

{
    "location": "library.html#CausalInference.randdag",
    "page": "Library",
    "title": "CausalInference.randdag",
    "category": "Function",
    "text": "randdag(n, alpha = 0.1)\n\nCreate random DAG from randomly permuted random triangular matrix with edge probability alpha.\n\n\n\n"
},

{
    "location": "library.html#CausalInference.disjoint_sorted",
    "page": "Library",
    "title": "CausalInference.disjoint_sorted",
    "category": "Function",
    "text": "disjoint_sorted(u, v)\n\nCheck if the intersection of sorted collections is empty. The intersection of empty collectios is empty.\n\n\n\n"
},

{
    "location": "library.html#CausalInference.ordered_edges",
    "page": "Library",
    "title": "CausalInference.ordered_edges",
    "category": "Function",
    "text": "ordered_edges(dag)\n\nIterator of edges of a dag, ordered in Chickering order:\n\nPerform a topological sort on the NODES\nwhile there are unordered EDGES in g\n    Let y be the lowest ordered NODE that has an unordered EDGE incident into it\n    Let x be the highest ordered NODE for which x => y is not ordered\n    return x => y \nend\n\n\n\n"
},

{
    "location": "library.html#CausalInference.insorted",
    "page": "Library",
    "title": "CausalInference.insorted",
    "category": "Function",
    "text": "insorted(a, x)\n\nCheck if x is in the sorted collection a\n\n\n\n"
},

{
    "location": "library.html#CausalInference.removesorted!",
    "page": "Library",
    "title": "CausalInference.removesorted!",
    "category": "Function",
    "text": "removesorted(collection, item) -> contains(collection, item)\n\nRemove item from sorted collection. \n\n\n\n"
},

{
    "location": "library.html#Miscellaneous-1",
    "page": "Library",
    "title": "Miscellaneous",
    "category": "section",
    "text": "digraph\npairs\npc_oracle\nskel_oracle\nranddag\nCausalInference.disjoint_sorted \nCausalInference.ordered_edges\nCausalInference.insorted\nCausalInference.removesorted!"
},

]}
