
@testset "test data input" begin

    @testset "json parser" begin
        data = _RDT.parse_file("../test/data/example.json")

        @test data["multinetwork"]
        @test haskey(data["nw"]["2"], "damaged_branch")
        @test !haskey(data, "scenarios")
    end

#    @testset "data input" begin
#         data = _PMs.parse_file("../test/data/case30.m")
         # checks for keys
#         @test "vm_vuf_max" in keys(data)
#         @test "beta_e" in keys(data)
#         @test "critical_load" in keys(data)
#         @test "non_critical_load" in keys(data)
         # checks data
#         @test data["vm_vuf_max"] == .03
#         @test data["beta_e"] == .15
#         @test data["critical_load"] == .9
#         @test data["non_critical_load"] == .5
#         @test data["critical_bus"]["1"]["critical"] == 0
#         @test data["critical_bus"]["2"]["critical"] == 1
#     end

     @testset "test ref input" begin
         data = _RDT.parse_file("../test/data/example.json")
         pm = _PMD.instantiate_mc_model(data, _PMD.LinDist3FlowPowerModel, build_mc_rdt; ref_extensions=[ref_add_rdt!],  eng2math_extensions=[_RDT.transform_switch_inline_ne!,_RDT.transform_switch_inline!], multinetwork=true)

         @test length(_PMD.ref(pm, 0, :damaged_branch)) == 0
         @test length(_PMD.ref(pm, 0, :undamaged_branch)) == length(_PMD.ref(pm, 0, :branch))

         @test length(_PMD.ref(pm, 1, :damaged_branch)) == 4
         @test length(_PMD.ref(pm, 1, :damaged_transformer)) == 0
         @test length(_PMD.ref(pm, 1, :damaged_branch_ne)) == 5
         @test length(_PMD.ref(pm, 1, :undamaged_branch)) == length(_PMD.ref(pm, 1, :branch)) - 4

         @test length(_PMD.ref(pm, 2, :damaged_branch)) == 2
         @test length(_PMD.ref(pm, 2, :damaged_transformer)) == 0
         @test length(_PMD.ref(pm, 2, :damaged_branch_ne)) == 8
         @test length(_PMD.ref(pm, 2, :undamaged_branch)) == length(_PMD.ref(pm, 2, :branch))  - 2

         @test length(_PMD.ref(pm, 0, :branch_ne)) == 28

     end

#     @testset "scenario input" begin
#         dataScen = _RTD.gen_scenarios("../test/data/scenData_ori.json",data);
#         @test !(1 in keys(dataScen[1]["branch"]))
#         @test !(4 in keys(dataScen[1]["branch"]))
#         @test !(2 in keys(dataScen[2]["branch"]))
#         @test !(3 in keys(dataScen[2]["branch"]))
#     end

end
