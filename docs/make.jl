using Documenter, Tulip

makedocs(
    sitename = "Tulip.jl",
    authors = "Mathieu Tanneau",
    pages = [
        "Home" => "index.md"
    ]
)

deploydocs(
    repo = "github.com/ds4dm/Tulip.jl.git"
)