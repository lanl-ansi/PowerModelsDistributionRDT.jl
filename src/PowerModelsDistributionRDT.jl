module PowerModelsDistributionRDT

    import JuMP
    import Memento

    import InfrastructureModels
    import PowerModels
    import PowerModelsDistribution

    const _PMs = PowerModels
    const _PMD = PowerModelsDistribution
    const _INs = InfrastructureModels

    function __init__()
        global _LOGGER = Memento.getlogger(PowerModels)
    end

    const TESTLOG = Memento.getlogger(PowerModels)
    Memento.setlevel!(TESTLOG, "error")

    include("core/ref.jl")
    include("core/variable.jl")
    include("core/constraint_template.jl")
    include("core/constraint.jl")
    include("core/data.jl")
    include("core/objective.jl")

    include("io/common.jl")
    include("io/json_parse.jl")

    include("prob/rdt.jl")

    include("core/export.jl")  # must be last include to properly export functions
end # module
