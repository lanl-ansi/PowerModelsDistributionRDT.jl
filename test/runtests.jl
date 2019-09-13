using PowerModelsDistributionRDT

import InfrastructureModels
import PowerModels
import PowerModelsDistribution

const _PMD = PowerModelsDistribution
const _PMs = PowerModels

using Test
using LinearAlgebra

import Ipopt
ipopt_solver = with_optimizer(Ipopt.Optimizer, tol=1e-6, print_level=0)

@testset "PowerModelsDistributionRDT" begin
    include("data_json.jl")

end
