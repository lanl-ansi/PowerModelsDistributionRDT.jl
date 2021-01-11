using Documenter
using PowerModelsDistributionRDT

makedocs(
    sitename = "PowerModelsDistributionRDT",
    format = Documenter.HTML(),
    modules = [PowerModelsDistributionRDT]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
