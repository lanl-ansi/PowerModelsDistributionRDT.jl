
@testset "Small Test" begin
    data = _RDT.parse_file(small_data)
    result = _RDT.solve_rdt(data, _PMD.ACPUPowerModel, juniper_solver)
    @test result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED
    println(result["termination_status"])
end
