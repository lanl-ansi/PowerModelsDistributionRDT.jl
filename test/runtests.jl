using PowerModelsDistributionRDT

const _RDT = PowerModelsDistributionRDT

import InfrastructureModels
const _IM = InfrastructureModels

import Memento

import PowerModels
const _PM = PowerModels

import PowerModelsDistribution
const _PMD = PowerModelsDistribution

# Suppress warnings during testing.
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
PowerModels.logger_config!("error")

import Ipopt
import SCS
import Juniper
import HiGHS
import SCIP

import JuMP
import JSON

import LinearAlgebra
using Test


using Test

_PMD.silence!()

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer,"print_level" => 0,"sb" => "yes","max_iter" => 1000,"acceptable_tol" => 1.0e-2)
juniper_solver = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt_solver,"log_levels" => [],)
highs_solver = JuMP.optimizer_with_attributes(HiGHS.Optimizer)
scip_solver = JuMP.optimizer_with_attributes(SCIP.Optimizer)

@testset "microgrid" begin
    include("data.jl")
    include("transform.jl")
    include("rdt.jl")
end
