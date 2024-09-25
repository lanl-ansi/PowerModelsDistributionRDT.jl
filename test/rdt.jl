
@testset "Small Test" begin
    data = parse_file(small_data)
    result = solve_rdt(data, ACPUPowerModel, juniper_solver)
    @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
    @test isapprox(result["objective"], 167, atol=1e-1)
end

@testset "Small Test - Cost Adjust 1" begin
    # adjust the cost to build a DER, so the better solution is to build the DER and harden the line between 3 and 4 (allowing line 1 to fail in one scenario)
    data = parse_file(small_data)
    data["nw"]["0"]["branch"]["1"]["harden_cost"] = 1000.0
    data["nw"]["0"]["gen_ne"]["2"]["construction_cost"] = 10.0
    result = solve_rdt(data, ACPUPowerModel, juniper_solver)
    @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
    @test isapprox(result["objective"], 313, atol=1e-1)
end

@testset "Small Test - Cost Adjust 2" begin
    #adjust the cost, so the better solution is to build the expansion branch between 1 and 4, and needs a new switch for radiol operations under the base scenario
    data = parse_file(small_data)
    data["nw"]["0"]["branch_ne"]["4"]["construction_cost"] = 7.0
    result = solve_rdt(data, ACPUPowerModel, juniper_solver)
    @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
    @test isapprox(result["objective"], 167, atol=1e-1)
end
