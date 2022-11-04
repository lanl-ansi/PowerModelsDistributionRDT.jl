
@testset "transform inline switch (existing line)" begin
    data = _RDT.parse_file("../test/data/transform_test.json")
    pm = _PMD.instantiate_mc_model(data, _PMD.LPUBFDiagPowerModel, build_mc_rdt; ref_extensions=[ref_add_rdt!], multinetwork=true)

    @test length(_PMD.ref(pm, :bus)) == 5
    @test length(_PMD.ref(pm, :switch_inline_ne)) == 1
    @test length(_PMD.ref(pm, :branch_harden)) == 1
    @test length(_PMD.ref(pm, :branch_ne)) == 1
    @test length(_PMD.ref(pm, :gen_ne)) == 1
    @test length(_PMD.ref(pm, :switch)) == 2
end
