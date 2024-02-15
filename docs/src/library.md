# Library

## (Partially) directed acyclic graphs (PDAGs and DAGs)
```@docs
has_a_path
CausalInference.graph
CausalInference.isclique
CausalInference.parents
CausalInference.children
CausalInference.isundirected
CausalInference.isparent
CausalInference.neighbors_undirected
CausalInference.isoriented
CausalInference.orientedge!
CausalInference.neighbors_adjacent
CausalInference.isadjacent
CausalInference.ischild
```

## Causal graphs
```@docs
dsep
cpdag
alt_cpdag
vskel
has_recanting_witness
backdoor_criterion
meek_rules!
CausalInference.meek_rule1
CausalInference.meek_rule2
CausalInference.meek_rule3
CausalInference.meek_rule4
pdag2dag!
pdag_to_dag_dortarsi!
pdag_to_dag_meek!
```

## PC algorithm

```@docs
pcalg
orient_unshielded
apply_pc_rules
skeleton
CausalInference.skeleton_stable
dseporacle
unshielded
orientable_unshielded
plot_pc_graph_tikz
plot_pc_graph_recipes
plot_pc_graph_text
CausalInference.prepare_pc_graph
```

## Statistics

```@docs
gausscitest
partialcor
cmitest
```

## KL Entropy Estimators
```@docs
n_ball
kl_entropy
kl_renyi
kl_mutual_information
kl_cond_mi
kl_perm_mi_test
kl_perm_cond_mi_test
```

## FCI algorithm
```@docs
has_marks
set_marks!
is_collider
is_parent
is_triangle
is_discriminating_path
is_uncovered_circle_path
is_uncovered_PD_path
fcialg
plot_fci_graph_tikz
plot_fci_graph_recipes
plot_fci_graph_text
CausalInference.prepare_fci_graph
```

## Miscellaneous
```@docs
digraph
vpairs
pc_oracle
skel_oracle
randdag
CausalInference.disjoint_sorted
CausalInference.ordered_edges
CausalInference.removesorted!
graph_to_text
kwargs_pdag_graphmakie
CausalInference.combinations_without
```
## Bayes ball
```@docs
CausalInference.bayesball
CausalInference.bayesball_graph
```

## Adjustment
```@docs
ancestors
descendants
alt_test_dsep
test_covariate_adjustment
alt_test_backdoor
find_dsep
find_min_dsep
find_covariate_adjustment
find_backdoor_adjustment
find_frontdoor_adjustment
find_min_covariate_adjustment
find_min_backdoor_adjustment
find_min_frontdoor_adjustment
list_dseps
list_covariate_adjustment
list_backdoor_adjustment
list_frontdoor_adjustment
```

## GES
```@docs
ges
CausalInference.local_score
CausalInference.Insert!
CausalInference.Delete!
```

## Causal Zig-Zag
```@docs
CausalInference.nup
CausalInference.ndown
CausalInference.count_moves_uniform
CausalInference.causalzigzag
CausalInference.keyedreduce
CausalInference.ne_total
CausalInference.arrows
CausalInference.remove!
CausalInference.count_moves
CausalInference.ne_undirected
CausalInference.operator_mcs
CausalInference.InsertIterator
CausalInference.DeleteIterator
CausalInference.count_mcs
CausalInference.precompute_semidirected
CausalInference.next_CPDAG
```

## DAG Zig-Zag
```@docs
CausalInference.dagzigzag
CausalInference.count_dag_moves 
```

## Exact
```@docs
CausalInference.exactscorebased
```
