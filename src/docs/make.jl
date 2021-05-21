using Documenter, PowerModelsDistributionRDT

makedocs(
    modules = [PowerModelsDistributionRDT],
    format = Documenter.HTML(analytics = ""),
    sitename = "PowerModelsDistributionRDT",
    authors = "Russell Bent, David Fobes, Jose Tabarez, and contributors.",
    pages = [
        "Home" => "index.md",
    ]
)

deploydocs(
    repo = "github.com/lanl-ansi/PowerModelsDistributionRDT.jl.git",
)
