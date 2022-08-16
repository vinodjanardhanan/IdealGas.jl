using IdealGas
using Documenter

DocMeta.setdocmeta!(IdealGas, :DocTestSetup, :(using IdealGas); recursive=true)

makedocs(;
    modules=[IdealGas],
    authors="Vinod Janardhanan",
    repo="https://github.com/vinodjanardhanan/IdealGas.jl/blob/{commit}{path}#{line}",
    sitename="IdealGas.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vinodjanardhanan.github.io/IdealGas.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/vinodjanardhanan/IdealGas.jl",
    devbranch="main",
)
