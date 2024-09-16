
@testset "Small Test" begin
    data = parse_file(small_data)
    result = solve_rdt(data, ACPUPowerModel, juniper_solver)
    @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
    @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
    @test isapprox(result["objective"], 167, atol=1e-1)
end
