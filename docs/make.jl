using SampledDomains
using Documenter, Literate

# The examples in `examples/` are Literate scripts; they become pages in `docs/src/examples/`
# (generated, not in git) and are run during the build, so their figures are always current.
examples_folder = joinpath(@__DIR__, "..", "examples")
docs_examples_folder = joinpath(@__DIR__, "src", "examples")
for f in readdir(examples_folder; join=true)
    endswith(f, ".jl") && Literate.markdown(f, docs_examples_folder)
end

makedocs(;
    modules=[SampledDomains],
    authors="Oleg Soloviev",
    repo=Remotes.GitHub("JuliaPhase", "SampledDomains.jl"),
    sitename="SampledDomains.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaphase.github.io/SampledDomains.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Arrays and pictures" => "examples/ArraysAndPictures.md",
        "Grids" => "examples/Grids.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaPhase/SampledDomains.jl",
)
