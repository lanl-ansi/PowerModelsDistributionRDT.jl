using PowerModelsDistributionRDT


include("../src/microgrid.jl")
const _RDT = RDT



import InfrastructureModels
import Memento

import PowerModels
const _PMs = PowerModels

import PowerModelsDistribution
const _PMD = PowerModelsDistribution

# Suppress warnings during testing.
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
PowerModels.logger_config!("error")

import Cbc
import Ipopt
import SCS
import Juniper

import JuMP
import JSON

import LinearAlgebra
using Test

# import Gurobi

# const GRB_ENV = Gurobi.Env()
# solver = JuMP.with_optimizer(Gurobi.Optimizer, GRB_ENV)

using Test

ipopt_solver = JuMP.with_optimizer(Ipopt.Optimizer, tol=1e-6, print_level=0)
cbc_solver = JuMP.with_optimizer(Cbc.Optimizer, logLevel=0)
juniper_solver = JuMP.with_optimizer(Juniper.Optimizer, nl_solver=_PMs.with_optimizer(Ipopt.Optimizer, tol=1e-4, print_level=0), mip_solver=cbc_solver, log_levels=[])
# solver = JuMP.with_optimizer(Gurobi.Optimizer, GRB_ENV)


@testset "microgrid" begin
    include("rdt.jl")
end