using Documenter
push!(LOAD_PATH, "..")
using MFDecoupling


makedocs(;
    modules=[MFDecoupling],
    authors="Julian Stobbe <Atomtomate@gmx.de> and contributors",
    repo="https://github.com/Atomtomate/MFDecoupling.jl/blob/{commit}{path}#L{line}",
    sitename="MFDecoupling",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://Atomtomate.github.io/MFDecoupling.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Dependencies" => "deps.md",
    ],
)

deploydocs(;
    branch="gh-pages",
    devbranch = "main",
    devurl = "stable",
    repo="github.com/Atomtomate/MFDecoupling.jl.git",
)
