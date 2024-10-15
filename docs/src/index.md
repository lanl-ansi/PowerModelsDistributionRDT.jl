# PowerModelsDistributionRDT.jl

```@meta
CurrentModule = PowerModelsDistributionRDT
```

## What is PowerModelsDistributionRDT?

[PowerModelsDistributionRDT.jl](https://github.com/lanl-ansi/PowerModelsDistributionRDT.jl) is a Julia/JuMP-based library for designing distribution networks, in particular unbalanced (i.e., multiconductor) power distribution networks.

## Resources for Getting Started

Read the [Installation Guide](@ref Installation-Guide)

Read the [Quickstart Guide](@ref Quick-Start-Guide)

Read the introductory tutorial [Introduction to PowerModelsDistributionRDT](@ref Introduction-to-PowerModelsDistributionRDT)

## How the documentation is structured

The following is a high-level overview of how our documentation is structured. There are three primary sections:

- The **Manual** contains detailed documentation for certain aspects of PowerModelsDistributionRDT, such as

- **Tutorials** contains working examples of how to use PowerModelsDistributionRDT. Start here if you are new to PowerModelsDistributionRDT.

- The **API Reference** contains a complete list of the functions you can use in PowerModelsDistributionRDT. Look here if you want to know how to use a particular function.

## PowerModelsDistributionRDT Analyses Packages

PowerModelsDistribution depends on several other PowerModels(...) packages from the InfrastructureModels ecosystem. The packages in blue below are created and maintained by the core InfrastructureModels developer team, and the other packages are those that are built as extensions or rely on one of the core InfrastructureModels packages in some way.

![InfrastructureModels Ecosystem](assets/infrastructuremodels_ecosystem.png)

### PowerModelsDistribution

[PowerModelsDistribution.jl](https://github.com/lanl-ansi/PowerModelsDistribution.jl) is a Julia/JuMP-based package for modeling unbalanced (i.e., multiconductor) power networks. This is the primary modeling framework utilized in PowerModelsDistributionRDT, and contains the primary logic for optimization and parsing of network data.

### PowerModelsONM

[PowerModelsONM.jl](https://github.com/lanl-ansi/PowerModelsONM.jl) is a Julia/JuMP-based package for operation and restoration of electric power distribution feeders featuring networked microgrids. Many of the formulations for optimal operations under extreme events are integrated into the design model of PowerModelDistributionRDT.


## License

This code is provided under a BSD license as part of the Multi-Infrastructure Control and Optimization Toolkit (MICOT) project, C15024.
