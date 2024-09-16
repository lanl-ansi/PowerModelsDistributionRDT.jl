using Documenter
using PowerModelsDistributionRDT


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
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo="github.com/lanl-ansi/PowerModelsDistributionRDT.jl.git",
    push_preview=false,
    devbranch="main",
)
