push!(LOAD_PATH, "../src/")

using Documenter
using Arrhenius

makedocs(sitename="Arrhenius.jl",
         modules  = [Arrhenius],
         pages=[
                "Home" => "index.md"
               ])

deploydocs(;
    repo="github.com/DENG-MIT/Arrhenius.jl",
)