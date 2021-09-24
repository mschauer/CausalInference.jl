# Library

## Directed acyclic graphs (DAGs)

```@docs
dsep
cpdag
vskel
has_recanting_witness
backdoor_criterion
has_a_path
```

## PC algorithm

```@docs
pcalg
orient_unshielded
apply_pc_rules
skeleton
dseporacle
unshielded
orientable_unshielded
plot_pc_graph
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
plot_fci_graph
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
CausalInference.insorted
CausalInference.removesorted!
```
