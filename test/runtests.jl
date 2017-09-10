using CausalInference
using LightGraphs
using Base.Test

include(joinpath("..", "docs", "make.jl"))


include("skeleton.jl")
include("dsep.jl")
include("pc.jl")
include("cpdag.jl")