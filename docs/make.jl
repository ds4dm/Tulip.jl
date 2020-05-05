using Documenter, Tulip

makedocs(
    sitename = "Tulip.jl",
    authors = "Mathieu Tanneau",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "Home" => "index.md",
        "User manual" => Any[
            "Problem formulation" => "manual/formulation.md",
            "Solving linear systems" => "manual/linear_systems.md"
        ],
        "Reference" => Any[
            "Julia API" => "reference/API.md",
            "Attributes" => "reference/attributes.md",
            "Parameters" => "reference/parameters.md",
            "KKT solvers" => "reference/kkt_solvers.md"
        ],
    ]
)

deploydocs(
    repo = "github.com/ds4dm/Tulip.jl.git"
)