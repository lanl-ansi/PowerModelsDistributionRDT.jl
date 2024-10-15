# Quick Start Guide

Once PowerModelsDistributionRDT is [installed](@ref Installation-Guide), to operate PowerModelsDistributionRDT several other things are required, at a minimum, a distribution data set in LPNORM [format](https://github.com/lanl-ansi/micot/wiki/Resilient-Design-Executable) and an optimization solver like Juniper is needed.  The basic RDT problems can be executed with


```julia
using PowerModelsDistributionRDT

import PowerModelsDistribution as _PMD
import Ipopt
import Juniper
import HiGHS


ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes", "max_iter" => 1000, "acceptable_tol" => 1.0e-2)
highs_solver = JuMP.optimizer_with_attributes(HiGHS.Optimizer, "small_matrix_value" => 1e-12, "output_flag"=>false)
juniper_solver = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt_solver, "mip_solver" => highs_solver, "log_levels" => [],)

result = solve_rdt("test/data/small.json", _PMD.ACPUPowerModel, juniper_solver)
```

where ACPUPowerModel can be replaced with a desired distribution network formulation, typically defined in the PowerModelsDistribution package.

## Getting Results

The solve commands in PowerModelsDistributionRDT return detailed results data in the form of a dictionary.  This dictionary can be saved for further processing as follows,

```julia
result = olve_rdt("test/data/small.json", _PMD.ACPUPowerModel, juniper_solver)
```

For example, the algorithm's runtime and final objective value can be accessed with,

```
result["solve_time"]
result["objective"]
```

The `"solution"` field contains detailed information about the solution produced by the solve method.
For example, the following dictionary comprehension can be used to inspect the choice to build a new branch in the solution,

```
Dict(name => data["ze"] for (name, data) in result["solution"]["branch_ne"])
```

The `print_summary(result["solution"])` function can be used show an table-like overview of the solution data.  
