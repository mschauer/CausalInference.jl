module GraphMakieExt

import CausalInference

CausalInference.EXTENSIONS_SUPPORTED ? (using GraphMakie) : (using ..GraphMakie)

import CausalInference.kwargs_pdag_graphmakie


function CausalInference.kwargs_pdag_graphmakie(g::GraphMakie.Observable; ilabels=1:CausalInference.nv(g[]), arrowsize=25, ilabels_fontsize=25)
    
    kwargs = @GraphMakie.lift kwargs_pdag_graphmakie($g; ilabels, arrowsize, ilabels_fontsize)
    arrow_size = @GraphMakie.lift ($kwargs).arrow_size
    edge_width = @GraphMakie.lift ($kwargs).edge_width
    (; arrow_shift=:end, arrow_size, edge_width)
end

end