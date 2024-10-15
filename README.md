# PowerModelsDistributionRDT

|                                      **Documentation**                                       |                                          **Build Status**                                          |
| :------------------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------: |
| [![docs-stable][docs-stable-img]][docs-stable-url] [![docs-dev][docs-dev-img]][docs-dev-url] | [![github-actions][github-actions-img]][github-actions-url] [![codecov][codecov-img]][codecov-url] |

This package combines various packages in the [InfrastructureModels.jl](https://github.com/lanl-ansi/InfrastructureModels.jl) optimization library ecosystem, particularly those related to electric power distribution.

PowerModelsDistributionRDT focuses on optimal design of phase unbalanced (multiconductor) distribution feeders, primarily for meeting resilience metrics. Phase unbalanced modeling is achieved using [PowerModelsDistribution](https://github.com/lanl-ansi/PowerModelsDistribution.jl). This library features a custom implementation of network design for resilience. See [documentation][docs-stable-url] for more details.

## Installation

To install PowerModelsDistributionRDT, use the built-in Julia package manager

```julia
pkg> add PowerModelsDistributionRDT
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("PowerModelsDistributionRDT")
```

or to develop the package,

```julia
julia> import Pkg; Pkg.develop(Pkg.PackageSpec(; name="PowerModelsDistributionRDT", url="https://github.com/lanl-ansi/PowerModelsDistributionRDT.jl"))
```

## Questions and contributions

Usage questions can be posted on the [Github Discussions forum][discussions-url].

Contributions, feature requests, and suggestions are welcome; please open an [issue][issues-url] if you encounter any problems. The [contributing page][contrib-url] has guidelines that should be followed when opening pull requests and contributing code.

This software was supported by the LPNORM project funded by the U.S. Department of Energy's Microgrid Research and Development Program.

## Citing PowerModelsDistributionRDT

If you find PowerModelsDistributionRDT useful for your work, we kindly request that you cite the following [publication](https://doi.org/10.1287/ijoc.2019.0899):

```bibtex
@article{byeon2020communication,
  title={Communication-constrained expansion planning for resilient distribution systems},
  author={Byeon, Geunyeong and Van Hentenryck, Pascal and Bent, Russell and Nagarajan, Harsha},
  journal={INFORMS Journal on Computing},
  volume={32},
  number={4},
  pages={968--985},
  year={2020},
  publisher={INFORMS}
}
```

which first introduced the problem specification implemented in this software.

## License

This code is provided under a BSD license as part of the Multi-Infrastructure Control and Optimization Toolkit (MICOT) project, C15024.

[docs-dev-img]: https://github.com/lanl-ansi/PowerModelsDistributionRDT.jl/workflows/Documentation/badge.svg
[docs-dev-url]: https://lanl-ansi.github.io/PowerModelsDistributionRDT.jl/dev
[docs-stable-img]: https://github.com/lanl-ansi/PowerModelsDistributionRDT.jl/workflows/Documentation/badge.svg
[docs-stable-url]: https://lanl-ansi.github.io/PowerModelsDistributionRDT.jl/stable
[github-actions-img]: https://github.com/lanl-ansi/PowerModelsDistributionRDT.jl/workflows/CI/badge.svg
[github-actions-url]: https://github.com/lanl-ansi/PowerModelsDistributionRDT.jl/actions/workflows/ci.yml
[codecov-img]: https://codecov.io/gh/lanl-ansi/PowerModelsDistributionRDT.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/lanl-ansi/PowerModelsDistributionRDT.jl
[contrib-url]: https://lanl-ansi.github.io/PowerModelsDistributionRDT.jl/stable/developer/contributing.html
[discussions-url]: https://github.com/lanl-ansi/PowerModelsDistributionRDT.jl/discussions
[issues-url]: https://github.com/lanl-ansi/PowerModelsDistributionRDT.jl/issues
