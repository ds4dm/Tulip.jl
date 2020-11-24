using Documenter, Tulip

makedocs(
    sitename = "Tulip.jl",
    authors = "Mathieu Tanneau",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
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
            "KKT solvers" => "reference/KKT/kkt_solvers.md",
            "Julia API" => "reference/API.md",
            "Attributes" => "reference/attributes.md",
            "Options" => "reference/options.md"
        ],
    ]
)

deploydocs(
    repo = "github.com/ds4dm/Tulip.jl.git"
)