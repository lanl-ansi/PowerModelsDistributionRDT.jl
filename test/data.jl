
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


     @testset "test load block" begin
         data = _RDT.parse_file(small_data)
         pm = _PMD.instantiate_mc_model(data, _PMD.ACPUPowerModel, build_mc_rdt; ref_extensions=[ref_add_rdt!],  eng2math_extensions=[_RDT.transform_switch_inline_ne!,_RDT.transform_switch_inline!], multinetwork=true)

         @test length(_PMD.ref(pm, 0, :blocks)) == 5
         @test length(_PMD.ref(pm, 1, :blocks)) == 6
         @test length(_PMD.ref(pm, 2, :blocks)) == 6

        found_8   = false
        found_1_5 = false
        found_2_6 = false
        found_3_7 = false
        found_4   = false
        for block in values(_PMD.ref(pm, 0, :blocks))
            if length(block) == 1 && 8 in block
                found_8 = true
            end
            if length(block) == 1 && 4 in block
                found_4 = true
            end
            if length(block) == 2 && 1 in block && 5 in block
                found_1_5 = true
            end
            if length(block) == 2 && 2 in block && 6 in block
                found_2_6 = true
            end
            if length(block) == 2 && 3 in block && 7 in block
                found_3_7 = true
            end
        end
        @test found_8   == true
        @test found_1_5 == true
        @test found_2_6 == true
        @test found_3_7 == true
        @test found_4   == true


     end

end
