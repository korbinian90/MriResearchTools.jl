using MriResearchTools
using Documenter

DocMeta.setdocmeta!(MriResearchTools, :DocTestSetup, :(using MriResearchTools); recursive=true)

makedocs(;
    modules=[MriResearchTools],
    authors="Korbinian Eckstein korbinian90@gmail.com",
    repo="https://github.com/korbinian90/MriResearchTools.jl/blob/{commit}{path}#{line}",
    sitename="MriResearchTools.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://korbinian90.github.io/MriResearchTools.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/korbinian90/MriResearchTools.jl",
    devbranch="master",
)
