using PowerModelsDistributionRDT

const _RDT = PowerModelsDistributionRDT

import InfrastructureModels
const _IM = InfrastructureModels

import Memento

import PowerModelsDistribution as _PMD
import PowerModelsONM as _ONM

# Suppress warnings during testing.
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")

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
_ONM.silence!()

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes", "max_iter" => 1000, "acceptable_tol" => 1.0e-2)
juniper_solver = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt_solver, "log_levels" => [],)
highs_solver = JuMP.optimizer_with_attributes(HiGHS.Optimizer)
scip_solver = JuMP.optimizer_with_attributes(SCIP.Optimizer)

example_data = "../test/data/example.json"
small_data = "../test/data/small.json"
transform_test_data = "../test/data/transform_test.json"


@testset "microgrid" begin
    include("data.jl")
    include("transform.jl")
    include("rdt.jl")
end
