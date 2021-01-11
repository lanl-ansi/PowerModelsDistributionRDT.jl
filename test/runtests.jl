using PowerModelsDistributionRDT

import InfrastructureModels
import PowerModels
import PowerModelsDistribution

import Ipopt
ipopt_solver = optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>0)

using Test

@testset "PowerModelsDistributionRDT" begin

end
