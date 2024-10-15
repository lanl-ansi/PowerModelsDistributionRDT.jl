module PowerModelsDistributionRDT
    import JuMP

    # InfrastructureModels ecosystem
    import InfrastructureModels as _IM
    import InfrastructureModels: ismultinetwork, ismultiinfrastructure

    import PowerModelsDistribution as _PMD
    import PowerModelsDistribution: ref, var, con, sol, ids, nw_ids, nw_id_default, nws, pmd_it_sym
    import PowerModelsDistribution: AbstractUnbalancedPowerModel, ACRUPowerModel, ACPUPowerModel, IVRUPowerModel, LPUBFDiagPowerModel, LinDist3FlowPowerModel, NFAUPowerModel

    import PowerModelsONM as _PMONM

    import SHA

    import Logging
    import LoggingExtras

    function __init__()
        global _LOGGER = Logging.ConsoleLogger(; meta_formatter=_PMD._pmd_metafmt)
        global _DEFAULT_LOGGER = Logging.current_logger()

        Logging.global_logger(_LOGGER)

        _PMD.set_logging_level!(:Info)
        _PMONM.set_log_level!(:Info)
    end

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
    include("form/apo.jl")
    include("form/lindistflow.jl")

    include("core/export.jl")  # must be last include to properly export functions
end # module
