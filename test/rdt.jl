
@testset "data input" begin
    data = _RDT.parse_file("../test/data/example.json")
    # ["bus", "source_type", "name", "dcline", "source_version", "branch", "gen", "storage", "switch", "disabled_lines", "baseMVA", "hardened_disabled_lines", "conductors", "data_model", "shunt", "transformer", "load"]
    result = _RDT.solve_rdt(data, _PMD.LPUBFDiagPowerModel, juniper_solver)
    # ipopt_solver
    @test result["termination_status"] == _PMs.LOCALLY_SOLVED || result["termination_status"] == _PMs.ALMOST_LOCALLY_SOLVED
end
