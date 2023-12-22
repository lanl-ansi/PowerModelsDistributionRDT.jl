
@testset "test data input" begin

    @testset "json parser" begin
        data = _RDT.parse_file(example_data)

        @test data["multinetwork"]
        @test haskey(data["nw"]["2"], "damaged_branch")
        @test !haskey(data, "scenarios")
    end

     @testset "test ref input" begin
         data = _RDT.parse_file(example_data)
         pm = _PMD.instantiate_mc_model(data, _PMD.ACPUPowerModel, build_mc_rdt; ref_extensions=[ref_add_rdt!],  eng2math_extensions=[_RDT.transform_switch_inline_ne!,_RDT.transform_switch_inline!], multinetwork=true)

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

end
