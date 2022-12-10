module TikzGraphsExt

import CausalInference

CausalInference.EXTENSIONS_SUPPORTED ? (using TikzGraphs) : (using ..TikzGraphs)

"""
    plot_pc_graph(g, df)

plot the output of the PC algorithm (TikzGraphs backend)
"""
function CausalInference.plot_pc_graph(g, node_labels::Array=[])
    objs = CausalInference.prepare_pc_graph(g, node_labels)
    TikzGraphs.plot(objs.plot_g, objs.node_labels, edge_styles=objs.edge_styles,
        node_style=objs.node_style, options=objs.options)
end

"""
    plot_fci_graph(g, node_labels)

plot the output of the FCI algorithm (TikzGraphs backend)
"""
function plot_fci_graph(g, node_labels::Array=[])
    objs = CausalInference.prepare_fci_graph(g, node_labels)
    TikzGraphs.plot(objs.plot_g, objs.node_labels, edge_styles=objs.edge_styles,
        node_style=objs.node_style, options=objs.options)
end

end # end of module