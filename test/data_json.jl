@testset "test json data import" begin
    files = readdir("../test/data/")
    counter = 1
    for file in files
        @testset "$file" begin
            mp_data = parse_file(string("../test/data/",file))
            println("percentage completed: ",counter/length(files)*100, "%")
            counter += 1  
        end
    end
    @testset "Ice_Harden_Rural_1" begin
        mp_data = parse_file("../test/data/Ice_Harden_Rural_10.json")
        result = _PMD.run_tp_pf(mp_data, _PMs.ACPPowerModel, ipopt_solver)
        @test result["termination_status"] == _PMs.LOCALLY_SOLVED
    end
    @testset "network123_0." begin
        mp_data = parse_file("../test/data/network123_0.json")
        result = _PMD.run_tp_pf(mp_data, _PMs.ACPPowerModel, ipopt_solver)
        @test result["termination_status"] == _PMs.LOCALLY_SOLVED
    end
    @testset "Ice_Harden_Rural_10" begin
        mp_data = parse_file("../test/data/Ice_Harden_Rural_10.json")
        result = _PMD.run_tp_pf(mp_data, _PMs.ACPPowerModel, ipopt_solver)
        @test result["termination_status"] == _PMs.LOCALLY_SOLVED
    end
end