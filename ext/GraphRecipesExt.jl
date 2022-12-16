module GraphRecipesExt

import CausalInference

CausalInference.EXTENSIONS_SUPPORTED ? (using GraphRecipes) : (using ..GraphRecipes)

function CausalInference.plot_pc_graph_recipes(g,
                                               node_labels::AbstractVector{<:AbstractString
                                                                           } = String[])
    objs = CausalInference.prepare_pc_graph(g, node_labels)
    # use the original graph to plot the edges as we cannot change arrow direction
    GraphRecipes.graphplot(g, names = objs.node_labels,
                           curves = true, nodeshape = :circle)
end

function CausalInference.plot_fci_graph_recipes(g,
                                                node_labels::AbstractVector{<:AbstractString
                                                                            } = String[])
    objs = CausalInference.prepare_fci_graph(g, node_labels)
    # use the original graph to plot the edges as we cannot change arrow direction
    GraphRecipes.graphplot(g, names = objs.node_labels,
                           curves = true, nodeshape = :circle)
end

end # end of module
