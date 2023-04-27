using PolarizedTypes
using Documenter

DocMeta.setdocmeta!(PolarizedTypes, :DocTestSetup, :(using PolarizedTypes); recursive=true)

makedocs(;
    modules=[PolarizedTypes],
    authors="Paul Tiede <ptiede91@gmail.com> and contributors",
    repo="https://github.com/ptiede/PolarizedTypes.jl/blob/{commit}{path}#{line}",
    sitename="PolarizedTypes.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ptiede.github.io/PolarizedTypes.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ptiede/PolarizedTypes.jl",
    devbranch="main",
)