push!(LOAD_PATH,"../src/")
using Documenter
using CausalInference

makedocs(modules = [CausalInference],
         sitename = "CausalInference.jl",
         authors = "Moritz Schauer and contributors",
         pages = Any["Home" => "index.md",
                     "Examples" => Any["examples/pc_basic_examples.md",
                                       "examples/pc_cmi_examples.md",
                                       "examples/pc_real_example.md",
                                       "examples/backdoor_example.md"],
                     "Library" => "library.md"])

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
# deploydocs(repo = "github.com/mschauer/CausalInference.jl.git")
