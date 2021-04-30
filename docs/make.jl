push!(LOAD_PATH, "../src/")

using Documenter
using Arrhenius

makedocs(
    format = Documenter.HTML(),
    sitename = "Arrhenius.jl",
    modules  = [Arrhenius],
    pages = [
        "Home" => "index.md"
    ]
)

deploydocs(;
    repo = "github.com/SuXY15/Arrhenius.jl",
    target = "build",
    branch = "gh-pages",
    devbranch = "main",
    devurl = "dev",
)