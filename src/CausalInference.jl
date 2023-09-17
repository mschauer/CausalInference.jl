module CausalInference

using Graphs
using Graphs.SimpleGraphs
using Combinatorics
using Base.Iterators
using Memoization, LRUCache
using ThreadsX
using LinkedLists

import Base: ==, show

export ancestors, descendants, alt_test_dsep, test_covariate_adjustment, alt_test_backdoor, find_dsep, find_min_dsep, find_covariate_adjustment, find_backdoor_adjustment, find_frontdoor_adjustment, find_min_covariate_adjustment, find_min_backdoor_adjustment, find_min_frontdoor_adjustment, list_dseps, list_covariate_adjustment, list_backdoor_adjustment, list_frontdoor_adjustment
export dsep, skeleton, gausscitest, dseporacle, partialcor
export unshielded, pcalg, vskel, vskel!, alt_vskel
export cpdag, alt_cpdag, meek_rules!
export digraph, vpairs, skel_oracle, pc_oracle, randdag
export cmitest, kl_entropy, kl_renyi, kl_mutual_information
export kl_cond_mi, kl_perm_mi_test, kl_perm_cond_mi_test
export n_ball
export fcialg, is_collider, is_triangle, is_parent
export is_discriminating_path, has_marks, set_marks!, is_uncovered_circle_path
export is_uncovered_PD_path, @arrow_str
export plot_pc_graph_text, plot_fci_graph_text, kwargs_pdag_graphmakie
export plot_pc_graph_recipes, plot_fci_graph_recipes # if GraphRecipes is loaded
export plot_pc_graph_tikz, plot_fci_graph_tikz # if TikzGraphs is loaded
export orient_unshielded, orientable_unshielded, apply_pc_rules
export ges
export pdag2dag!, pdag_to_dag_meek!, pdag_to_dag_dortasi!

#include("pinv.jl")
include("graphs.jl")
include("combinations_without.jl")
include("klentropy.jl")
include("skeleton.jl")
include("pdag.jl")
include("dsep.jl")
include("meek.jl")
include("pc.jl")
include("cpdag.jl")
include("fci.jl")
include("misc.jl")
include("recantingwitness.jl")
include("backdoor.jl")
include("ges.jl")
include("gensearch.jl")
include("workloads.jl")
#include("sampler.jl")
#include("mcs.jl")

# Compatibility with the new "Package Extensions" (https://github.com/JuliaLang/julia/pull/47695)
const EXTENSIONS_SUPPORTED = isdefined(Base, :get_extension)

if !EXTENSIONS_SUPPORTED
    using Requires: @require
end

function __init__()
    @static if !EXTENSIONS_SUPPORTED
        # requires both GraphRecipes.jl and Plots.jl
        @require GraphRecipes="bd48cda9-67a9-57be-86fa-5b3c104eda73" begin @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("../ext/GraphRecipesExt.jl") end
        @require TikzGraphs="b4f28e30-c73f-5eaf-a395-8a9db949a742" include("../ext/TikzGraphsExt.jl")
    end
end
end # end of module
