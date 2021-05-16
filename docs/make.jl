push!(LOAD_PATH, "../src/")

using Documenter
using Arrhenius

makedocs(
    format=Documenter.HTML(),
    sitename="Arrhenius.jl",
    modules=[Arrhenius],
    pages=[
        ##############################################
        ## MAKE SURE TO SYNC WITH docs/src/index.md ##
        ##############################################
        "Basics" => [
            "index.md",
            "install.md",
            "get_started.md",
            "concepts.md"
           ],
        "tutorial.md",
        "faq.md",
        "api.md"
    ],
)

deploydocs(;
    repo="github.com/DENG-MIT/Arrhenius.jl",
    devbranch="main"
)