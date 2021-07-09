using Documenter, Tulip

const _FAST = findfirst(isequal("--fast"), ARGS) !== nothing

makedocs(
    sitename = "Tulip.jl",
    authors = "Mathieu Tanneau",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    doctest = !_FAST,
    pages = [
        "Home" => "index.md",
        "Tutorials" => Any[
            "tutorials/lp_example.md",
        ],
        "User manual" => Any[
            "Problem formulation" => "manual/formulation.md",
            "Solving linear systems" => "manual/linear_systems.md",
            "Algorithms" => Any[
                "Homogeneous Self-Dual" => "manual/IPM/HSD.md",
                "Predictor-Corrector" => "manual/IPM/MPC.md"
            ],
            "Setting options" => "manual/options.md"
        ],
        "Reference" => Any[
            "Presolve" => "reference/Presolve/presolve.md",
            "KKT" => [
                "reference/KKT/kkt_systems.md",
                "reference/KKT/kkt_solvers.md",
            ],
            "Julia API" => "reference/API.md",
            "Attributes" => "reference/attributes.md",
            "Options" => "reference/options.md"
        ],
    ]
)

deploydocs(
    repo = "github.com/ds4dm/Tulip.jl.git"
)
