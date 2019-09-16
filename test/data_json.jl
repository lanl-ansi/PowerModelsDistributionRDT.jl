@testset "test json data import" begin
    files = readdir("../test/data/")
    counter = 1
    for file in files
        @testset "$file" begin
            mp_data = parse_file(string("../test/data/",file))
            # checks data required for PM was imported 
            result = _PMD.run_tp_pf(mp_data, _PMs.ACPPowerModel, ipopt_solver)
            @test result["termination_status"] == _PMs.LOCALLY_SOLVED

            # checking for RDT required data was imported 
            @test haskey(mp_data, "phase_variation")
            @test haskey(mp_data, "total_load_met")
            @test haskey(mp_data, "critical_load_met")
            println("percentage completed: ",counter/length(files)*100, "%")
            counter += 1  
        end
    end
end