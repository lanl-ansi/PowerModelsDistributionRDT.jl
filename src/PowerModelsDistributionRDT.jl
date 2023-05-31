module PowerModelsDistributionRDT

    import JuMP
    import Memento

    import InfrastructureModels
    import PowerModels
    import PowerModelsDistribution

    const _PM  = PowerModels
    const _PMD = PowerModelsDistribution
    const _IM  = InfrastructureModels

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
    include("core/eng2math.jl")

    include("io/common.jl")
    include("io/json_parse.jl")

    include("prob/rdt.jl")

    include("form/dcp.jl")
    include("form/shared.jl")
    include("form/acr.jl")
    include("form/acp.jl")

    include("core/export.jl")  # must be last include to properly export functions
end # module
