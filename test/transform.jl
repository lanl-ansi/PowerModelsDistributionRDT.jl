
@testset "transform inline switch (existing line)" begin
    data = parse_file(transform_test_data)
    pm = _PMD.instantiate_mc_model(data, _PMD.ACPUPowerModel, build_mc_rdt; ref_extensions=[ref_add_rdt!], multinetwork=true)

    @test length(ref(pm, :bus)) == 5
    @test length(ref(pm, :switch_inline_ne)) == 1 # adding a serial switch option for a line where we could add a switch to it
    @test length(ref(pm, :branch_harden)) == 1
    @test length(ref(pm, :branch_ne)) == 1
    @test length(ref(pm, :gen_ne)) == 1
    @test length(ref(pm, :switch)) == 2 # one switch added for an existing line that "has_switch = true" and one implict switch added for an expansion edge

    # TODO: Add some tests for the topology
end

# TODO: Add some tests for connected components with damaged lines, test it for inclusion of inline_switches_ne and branches_ne
