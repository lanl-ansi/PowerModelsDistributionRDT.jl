# Installation Guide

From Julia, PowerModelsDistributionRDT is installed using the built-in package manager:

```julia
]add PowerModelsDistributionRDT
```

or equivalently,

```julia
import Pkg
Pkg.add("PowerModelsDistributionRDT")
```

## Developer Installation

To install PowerModelsDistributionRDT as a developer,

```julia
import Pkg
Pkg.develop(Pkg.PackageSpec(; name="PowerModelsDistributionRDT", url="https://github.com/lanl-ansi/PowerModelsDistributionRDT.jl"))
```

From the command-line, outside Julia, one could download the repository, either via Github.com, or using git, _i.e._,

```sh
git clone https://github.com/lanl-ansi/PowerModelsDistributionRDT.jl.git
git checkout tags/v1.0.0
```

Then to install PowerModelsDistributionRDT and its required packages

```sh
julia --project="path/to/PowerModelsDistributionRDT" -e 'using Pkg; Pkg.instantiate(); Pkg.precompile();'
```
