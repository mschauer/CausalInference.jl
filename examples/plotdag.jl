using TikzGraphs
using TikzPictures

"""
    plot_dag(g, node_labels)

Plot DAG or the output of the PC algorithm.
"""
function plot_dag(g, node_labels::Array=[])
    plot_g = DiGraph(nv(g))

    if length(node_labels) != nv(g)
        node_labels = map(string, 1:nv(g))
    end

    node_style = "draw, rounded corners, fill=blue!10"
    options = "scale=2"

    styles_dict = Dict()

    for e in edges(g)
       if has_edge(plot_g, e.dst, e.src)
            styles_dict[(e.dst, e.src)] = "o-o"
       else
           add_edge!(plot_g, e.src, e.dst)
           push!(styles_dict, (e.src, e.dst)=>"->")
       end
    end

    TikzGraphs.plot(plot_g, node_labels, edge_styles=styles_dict,
                    node_style=node_style, options=options)
end
