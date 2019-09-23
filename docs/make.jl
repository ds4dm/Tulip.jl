using Documenter, Tulip

makedocs(
    sitename = "Tulip.jl",
    authors = "Mathieu Tanneau",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "Home" => "index.md",
        "Problem formulation" => "formulation.md"
    ]
)

deploydocs(
    repo = "github.com/ds4dm/Tulip.jl.git"
)