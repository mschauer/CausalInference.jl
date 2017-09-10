using Documenter
using CausalInference

makedocs(
    modules = [CausalInference],
    format = :html,
    sitename = "CausalInference.jl",
    authors = "Moritz Schauer and contributors",
    pages = Any[ # Compat: `Any` for 0.4 compat
        "Home" => "index.md",
#        "Manual" => "manual.md",
        "Library" => "library.md",
        ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/mschauer/CausalInference.jl.git",
    julia  = "0.6",
    target = "build",
    deps = nothing,
    make = nothing,
)