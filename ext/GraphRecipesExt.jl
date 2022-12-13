module GraphRecipesExt

import CausalInference

CausalInference.EXTENSIONS_SUPPORTED ? (using GraphRecipes) : (using ..GraphRecipes)

"""
    plot_pc_graph_recipes(g, node_labels::AbstractVector{<:AbstractString}=String[])

plot the output of the PC algorithm (GraphRecipes backend)
"""
function CausalInference.plot_pc_graph_recipes(g, node_labels::AbstractVector{<:AbstractString}=String[])
    objs = CausalInference.prepare_pc_graph(g, node_labels)
    GraphRecipes.graphplot(objs.plot_g, objs.node_labels)
    # TODO: add edge styles // communicate directedness
    # edge_styles=objs.edge_styles,node_style=objs.node_style, options=objs.options)
end

"""
    plot_fci_graph_recipes(g, node_labels::AbstractVector{<:AbstractString}=String[])

plot the output of the FCI algorithm (GraphRecipes backend)
"""
function CausalInference.plot_fci_graph_recipes(g, node_labels::AbstractVector{<:AbstractString}=String[])
    objs = CausalInference.prepare_fci_graph(g, node_labels)
    GraphRecipes.graphplot(objs.plot_g, objs.node_labels)
    # TODO: add edge styles // communicate directedness
    # edge_styles=objs.edge_styles,node_style=objs.node_style, options=objs.options)
end

end # end of module