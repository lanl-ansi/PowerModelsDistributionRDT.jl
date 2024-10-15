using Documenter
using PowerModelsDistributionRDT


pages = [
    "Introduction" => "index.md",
    "installation.md",
    "Manual" => [
        "Getting Started" => "manual/quickguide.md",
#        "The RDT Workflow" => "manual/rdt_workflow.md",
#        "Resilient Design Mathematical Model" => "manual/rdt_model.md",
    ],
    "Tutorials" => [
#        "Beginners Guide" => "tutorials/Beginners Guide.md",
#        "RDT Basic Example" => "tutorials/Block MLD Basic Example.md",
#        "JuMP Model by Hand - RDT Example" => "tutorials/JuMP Model by Hand - MLD-Block.md",
#        "Use Case Examples" => "tutorials/Use Case Examples.md",
    ],
#    "API Reference" => [
#        "Base functions" => "reference/base.md",
#        "Data Handling" => "reference/data.md",
#        "GraphML Functions" => "reference/graphml.md",
#        "Main Entrypoint" => "reference/entrypoint.md",
#        "Internal Functions" => "reference/internal.md",
#        "IO Functions" => "reference/io.md",
#        "Logging" => "reference/logging.md",
#        "Optimization Problems" => "reference/prob.md",
#        "Solution Statistics" => "reference/stats.md",
#        "Variables and Constraints" => "reference/variable_constraint.md",
#        "Types" => "reference/types.md",
#    ],
#    "Schemas" => schema_pages,
#    "Developer Docs" => [
#        "Contributing Guide" => "developer/contributing.md",
#        "Style Guide" => "developer/style.md",
#        "Roadmap" => "developer/roadmap.md",
#    ],
]


# build documents
makedocs(
    format=Documenter.HTML(
        analytics="",
        mathengine=Documenter.MathJax(),
        prettyurls=false,
        collapselevel=2,
    ),
    warnonly=true,
    sitename="PowerModelsDistributionRDT",
    authors="Russell Bent, Jose Tabarez, David M Fobes, and contributors",
    pages = pages
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo="github.com/lanl-ansi/PowerModelsDistributionRDT.jl.git",
    push_preview=false,
    devbranch="main",
)
