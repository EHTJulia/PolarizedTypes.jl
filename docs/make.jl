using PolarizedTypes
using Documenter

makedocs(;
    modules = [PolarizedTypes],
    authors = "Paul Tiede <ptiede91@gmail.com> and contributors",
    repo = "https://github.com/EHTJulia/PolarizedTypes.jl/blob/{commit}{path}#{line}",
    sitename = "PolarizedTypes.jl",
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo = "github.com/EHTJulia/PolarizedTypes.jl",
    push_preview = false,
    devbranch = "main",
)
