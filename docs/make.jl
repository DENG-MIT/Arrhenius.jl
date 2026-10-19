push!(LOAD_PATH, "../src/")

using Documenter
using Arrhenius

makedocs(
    repo=Documenter.Remotes.GitHub("DENG-MIT", "Arrhenius.jl"),
    format=Documenter.HTML(
        repolink="https://github.com/DENG-MIT/Arrhenius.jl",
        size_threshold=300 * 1024,
        size_threshold_warn=250 * 1024,
    ),
    sitename="Arrhenius.jl",
    modules=[Arrhenius],
    checkdocs=:exports,
    pages=[
        "Overview" => [
            "index.md",
            "install.md",
            "mechanisms.md",
            "concepts.md",
            "design_rationale.md",
        ],
        "reactors.md",
        "Flames" => ["flames.md", "solver_options.md", "counterflow_continuation.md"],
        "examples.md",
        "validation.md",
        "api.md",
    ],
)

deploydocs(;
    repo="github.com/DENG-MIT/Arrhenius.jl",
    devbranch="main"
)
