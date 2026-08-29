module TikzGraphsExt

import CausalInference

CausalInference.EXTENSIONS_SUPPORTED ? (using TikzGraphs) : (using ..TikzGraphs)

"""Reformat an input string to escape characters reserved for latex"""
function escape_latex_characters(input_string::String)

    simple_replacements = replace(input_string, r"([{}&#_%]{1})"=>s"\\\1")
    caret = replace(simple_replacements, r"\^"=>s"\\^{}")
    dollar = replace(caret, r"\$"=>s"\\$")
    backslash = replace(dollar, r"\\(?![$%^&{}_#])"=>s"")
    return backslash
end

function CausalInference.plot_pc_graph_tikz(g,
                                            node_labels::AbstractVector{<:AbstractString} = String[])
    objs = CausalInference.prepare_pc_graph(g, node_labels)
    TikzGraphs.plot(objs.plot_g,
                    [escape_latex_characters(label) for label in objs.node_labels],
                    edge_styles = objs.edge_styles, node_style = objs.node_style,
                    options = objs.options)
end

function CausalInference.plot_fci_graph_tikz(g,
                                             node_labels::AbstractVector{<:AbstractString} = String[])
    objs = CausalInference.prepare_fci_graph(g, node_labels)
    TikzGraphs.plot(objs.plot_g,
                    [escape_latex_characters(label) for label in objs.node_labels],
                    edge_styles = objs.edge_styles, node_style = objs.node_style,
                    options = objs.options)
end

end # end of module
