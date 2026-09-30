using SampledDomains
using Documenter

makedocs(;
    modules=[SampledDomains],
    authors="Oleg Soloviev",
    repo="https://github.com/JuliaPhase/SampledDomains.jl/blob/{commit}{path}#L{line}",
    sitename="SampledDomains.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaphase.github.io/SampledDomains.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaPhase/SampledDomains.jl",
)
