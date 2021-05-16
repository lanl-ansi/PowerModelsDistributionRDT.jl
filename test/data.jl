
@testset "test data input" begin
    @testset "json parser" begin
        data = _RTD.parse_file("../test/data/example.json")
        println(keys(data)) # ["total_load_met", "name", "critical_load_met", "multinetwork", "phase_variation", "nw", "chance_constraint", "per_unit"]
        println(keys(data["nw"]["1"])) # ["bus", "source_type", "name", "dcline", "source_version", "branch", "gen", "storage", "switch", "disabled_lines", "baseMVA", "hardened_disabled_lines", "conductors", "data_model", "shunt", "transformer", "load"]
        println(keys(data["nw"]["1"]["branch"]["1"]))
        # @test data["multinetwork"]
        # @test haskey(data["nw"]["1"], "disabled_lines")
        # @test !haskey(data, "scenarios")
    end
    # @testset "data input" begin
    #     data = _PMs.parse_file("../test/data/case30.m")
    #     # checks for keys
    #     @test "vm_vuf_max" in keys(data)
    #     @test "beta_e" in keys(data)
    #     @test "critical_load" in keys(data)
    #     @test "non_critical_load" in keys(data)
    #     # checks data
    #     @test data["vm_vuf_max"] == .03
    #     @test data["beta_e"] == .15
    #     @test data["critical_load"] == .9
    #     @test data["non_critical_load"] == .5
    #     @test data["critical_bus"]["1"]["critical"] == 0
    #     @test data["critical_bus"]["2"]["critical"] == 1
    # end
    
    # @testset "ref input" begin
    #     data = _PMs.parse_file("../test/data/case30.m")
    #     pm = _PMs.InitializePowerModel(_PMs.ACPPowerModel, data)
    #     _PMs.ref_add_core!(pm)
    #     ref_extensions = [_RTD.ref_add_load_weights!, _RTD.ref_add_vm_imbalance!]
    #     for ref_ext in ref_extensions
    #         ref_ext(pm)
    #     end
    # end
    # @testset "scenario input" begin
    #     
    #     dataScen = _RTD.gen_scenarios("../test/data/scenData_ori.json",data);
    #     @test !(1 in keys(dataScen[1]["branch"]))
    #     @test !(4 in keys(dataScen[1]["branch"]))
    #     @test !(2 in keys(dataScen[2]["branch"]))
    #     @test !(3 in keys(dataScen[2]["branch"]))
    # end
   
end
