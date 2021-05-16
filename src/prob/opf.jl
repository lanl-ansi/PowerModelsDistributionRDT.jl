

function run_opf(data::Dict{String,Any}, model_type, solver; kwargs...)
#    return _PMs.run_model(data, model_type, solver, build_mc_opf_org; multiconductor=true, ref_extensions=[_PMD.ref_add_arcs_trans!, ref_add_critical_leve!, ref_add_demand_level!, ref_add_vm_imbalance!], kwargs...)  #ref_add_load_weights!, ref_add_vm_imbalance!, ref_add_trans_imbalance!
    return _PMs.run_model(data, model_type, solver, build_mc_opf_org; multiconductor=true, ref_extensions=[_PMD.ref_add_arcs_trans!], kwargs...)
end

function build_mc_opf(pm::_PMs.AbstractPowerModel)

    _PMD.variable_mc_indicator_bus_voltage(pm)
    _PMD.variable_mc_bus_voltage_on_off(pm)

    _PMD.variable_mc_branch_flow(pm)
    _PMD.variable_mc_transformer_flow(pm)

    _PMD.variable_mc_indicator_generation(pm; relax=true)  
    _PMD.variable_mc_generation_on_off(pm)

    _PMD.variable_mc_indicator_demand(pm)
    _PMD.variable_mc_indicator_shunt(pm; relax=true)
    
    _PMD.constraint_mc_model_voltage(pm)
    # variable_v_imbalance(pm)

    for i in _PMs.ids(pm, :bus)
        constraint_mc_vm_vuf_mod(pm, i) 
    end

    for i in _PMs.ids(pm, :ref_buses)
        _PMD.constraint_mc_theta_ref(pm, i)
    end

    _PMD.constraint_mc_bus_voltage_on_off(pm)
    
    for i in _PMs.ids(pm, :gen)
        _PMs.constraint_generation_on_off(pm, i)
    end

    for i in _PMs.ids(pm, :bus)
        _PMD.constraint_mc_power_balance_shed(pm, i)
    end

    for i in _PMs.ids(pm, :branch)
        _PMD.constraint_mc_ohms_yt_from(pm, i)
        _PMD.constraint_mc_ohms_yt_to(pm, i)

        _PMD.constraint_mc_voltage_angle_difference(pm, i)

        _PMD.constraint_mc_thermal_limit_from(pm, i)
        _PMD.constraint_mc_thermal_limit_to(pm, i)
    end

    for i in _PMs.ids(pm, :transformer)
        _PMD.constraint_mc_trans(pm, i)
    end

    # constraint_critical_load(pm)
    constraint_non_critical_load(pm)

    _PMD.objective_mc_min_load_delta(pm)
end

function build_mc_opf_org(pm::_PMs.AbstractPowerModel)
    _PMD.variable_mc_voltage(pm)
    _PMD.variable_mc_branch_flow(pm)
    _PMD.variable_mc_transformer_flow(pm)
    _PMD.variable_mc_generation(pm)
    _PMD.variable_mc_load(pm)
    _PMD.variable_mc_storage(pm)

    _PMD.constraint_mc_model_voltage(pm)

    for i in _PMs.ids(pm, :ref_buses)
        _PMD.constraint_mc_theta_ref(pm, i)
    end

    # generators should be constrained before KCL, or Pd/Qd undefined
    for id in _PMs.ids(pm, :gen)
        _PMD.constraint_mc_generation(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in _PMs.ids(pm, :load)
        _PMD.constraint_mc_load(pm, id)
    end

    for i in _PMs.ids(pm, :bus)
        _PMD.constraint_mc_power_balance_load(pm, i)
    end

    for i in _PMs.ids(pm, :storage)
        _PMs.constraint_storage_state(pm, i)
        _PMs.constraint_storage_complementarity_nl(pm, i)
        _PMD.constraint_mc_storage_loss(pm, i)
        _PMD.constraint_mc_storage_thermal_limit(pm, i)
    end

    for i in _PMs.ids(pm, :branch)
        _PMD.constraint_mc_ohms_yt_from(pm, i)
        _PMD.constraint_mc_ohms_yt_to(pm, i)

        _PMD.constraint_mc_voltage_angle_difference(pm, i)

        _PMD.constraint_mc_thermal_limit_from(pm, i)
        _PMD.constraint_mc_thermal_limit_to(pm, i)
    end

    for i in _PMs.ids(pm, :transformer)
        _PMD.constraint_mc_trans(pm, i)
    end

    _PMs.objective_min_fuel_cost(pm)
end


 