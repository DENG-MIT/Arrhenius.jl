using Documenter, Arrhenius

makedocs(sitename="Arrhenius.jl",
         modules  = [Arrhenius],
         pages=[
                "Home" => "index.md"
               ])

deploydocs(;
    repo="github.com/USERNAME/Arrhenius.jl",
)