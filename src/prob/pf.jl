


function run_pf(data::Dict{String,Any}, model_type, solver; kwargs...)
    # return _PMs.run_model(data, model_type, solver, build_mc_pf; multiconductor=true, ref_extensions=[_PMD.ref_add_arcs_trans!, ref_add_transformer_imbalance!, ref_add_vm_imbalance!, ref_add_transformer_imbalance_arcs!], kwargs...) 
    return _PMs.run_model(data, model_type, solver, build_mc_pf; multiconductor=true, ref_extensions=[_PMD.ref_add_arcs_transformer!, ref_add_pq_imbalance!, ref_add_vm_imbalance!], multinetwork=false, kwargs...) 
end

function run_pf_bf(data::Dict{String,Any}, model_type, solver; kwargs...)
    # return _PMs.run_model(data, model_type, solver, build_mc_pf; multiconductor=true, ref_extensions=[_PMD.ref_add_arcs_trans!, ref_add_transformer_imbalance!, ref_add_vm_imbalance!, ref_add_transformer_imbalance_arcs!], kwargs...) 
    return _PMs.run_model(data, model_type, solver, build_mc_pf; multiconductor=true, kwargs...) 
end

function build_mc_pf(pm::_PMs.AbstractPowerModel)
    _PMD.variable_mc_bus_voltage(pm; bounded=false)
    _PMD.variable_mc_branch_power(pm; bounded=false)
    _PMD.variable_mc_transformer_power(pm; bounded=false)
    _PMD.variable_mc_gen_power_setpoint(pm; bounded=false)
    _PMD.variable_mc_load_setpoint(pm; bounded=false)
    variable_branch_be(pm)

    _PMD.constraint_mc_model_voltage(pm)

    for i in _PMs.ids(pm, :bus_bal)
        constraint_mc_vm_vuf(pm, i)
        # _PMD.constraint_mc_bus_voltage_balance(pm, i)
    end

    for (i,bus) in _PMs.ref(pm, :ref_buses)
        @assert bus["bus_type"] == 3

        _PMD.constraint_mc_theta_ref(pm, i)
        _PMD.constraint_mc_voltage_magnitude_only(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in _PMs.ids(pm, :gen)
        _PMD.constraint_mc_gen_setpoint(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in _PMs.ids(pm, :load)
        _PMD.constraint_mc_load_setpoint(pm, id)
    end

    for (i,bus) in _PMs.ref(pm, :bus)
        _PMD.constraint_mc_load_power_balance(pm, i)

        # PV Bus Constraints
        if length(_PMs.ref(pm, :bus_gens, i)) > 0 && !(i in _PMs.ids(pm,:ref_buses)) 
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2
            _PMD.constraint_mc_voltage_magnitude_only(pm, i)
            for j in _PMs.ref(pm, :bus_gens, i)
                _PMD.constraint_mc_gen_power_setpoint_real(pm, j)
            end
        end
    end

    for i in _PMs.ids(pm, :branch)
        _PMD.constraint_mc_ohms_yt_from(pm, i)
        _PMD.constraint_mc_ohms_yt_to(pm, i)
    end

    for i in _PMs.ids(pm, :transformer)
        _PMD.constraint_mc_transformer_power(pm, i)
    end
    
    for i in _PMs.ids(pm, :arcs_bal)
        constraint_branch_be(pm, i)
        constraint_balance_flow(pm, i)
    end
    
end
# function build_mc_pf(pm::_PMs.AbstractPowerModel)
#     # pm.data
#     _PMD.variable_mc_bus_voltage(pm; bounded=false)
#     _PMD.variable_mc_branch_flow(pm; bounded=false)
#     _PMD.variable_mc_transformer_flow(pm; bounded=false)
#     _PMD.variable_mc_generation(pm; bounded=false)
#     _PMD.variable_mc_load(pm; bounded=false)
#     _PMD.variable_mc_indicator_demand(pm; relax=true) # this may change to load in PMD 

#     # variable_branch_be(pm)
#     # variable_vm_imbalance(pm)

#     _PMD.constraint_mc_model_voltage(pm)

#     for (i,bus) in _PMs.ref(pm, :ref_buses)
#         @assert bus["bus_type"] == 3

#         _PMD.constraint_mc_theta_ref(pm, i)
#         _PMD.constraint_mc_voltage_magnitude_setpoint(pm, i)
#     end

#     # gens should be constrained before KCL, or Pd/Qd undefined
#     for id in _PMs.ids(pm, :gen)
#         _PMD.constraint_mc_generation(pm, id)
#     end

#     # loads should be constrained before KCL, or Pd/Qd undefined
#     for id in _PMs.ids(pm, :load)
#         _PMD.constraint_mc_load(pm, id)
#     end

#     for (i,bus) in _PMs.ref(pm, :bus)
#         _PMD.constraint_mc_power_balance_load(pm, i)

#         # PV Bus Constraints
#         if length(_PMs.ref(pm, :bus_gens, i)) > 0 && !(i in _PMs.ids(pm,:ref_buses))
#             # this assumes inactive generators are filtered out of bus_gens
#             @assert bus["bus_type"] == 2

#             _PMD.constraint_mc_voltage_magnitude_setpoint(pm, i)
#             for j in _PMs.ref(pm, :bus_gens, i)
#                 _PMD.constraint_mc_active_gen_setpoint(pm, j)
#             end
#         end
#     end

#     # for i in _PMs.ids(pm, :bus_bal)
#     #     # _PMD.constraint_mc_voltage_balance(pm, i)
#     #     constraint_mc_vm_vuf_mod(pm, i)
#     # end

#     for i in _PMs.ids(pm, :branch)
#         _PMD.constraint_mc_ohms_yt_from(pm, i)
#         _PMD.constraint_mc_ohms_yt_to(pm, i)
#     end

#     for i in _PMs.ids(pm, :transformer)
#         _PMD.constraint_mc_trans(pm, i)
#     end
#     # for i in _PMs.ids(pm, :arcs_bal)
#     #     constraint_branch_be(pm, i)
#     #     constraint_balance_flow(pm, i)
#     # end

#     # constraint_critical_load(pm)
#     # constraint_non_critical_load(pm)
# end

# function build_mc_pf(pm::_PMD.AbstractLPUBFModel)  
#     # Variables
#     _PMD.variable_mc_bus_voltage(pm; bounded=true) # TODO should be false
#     _PMD.variable_mc_branch_current(pm)
#     _PMD.variable_mc_branch_power(pm)
#     _PMD.variable_mc_transformer_power(pm; bounded=false)
#     _PMD.variable_mc_gen_power_setpoint(pm; bounded=false)
#     _PMD.variable_mc_load_setpoint(pm)

#     # Constraints
#     _PMD.constraint_mc_model_current(pm)

#     for (i,bus) in _PMs.ref(pm, :ref_buses)
#         if !(typeof(pm)<:_PMD.LPUBFDiagPowerModel)
#             _PMD.constraint_mc_theta_ref(pm, i)
#         end

#         @assert bus["bus_type"] == 3
#         _PMD.constraint_mc_voltage_magnitude_only(pm, i)
#     end

#     for id in _PMs.ids(pm, :gen)
#         _PMD.constraint_mc_gen_setpoint(pm, id)
#     end

#     for id in _PMs.ids(pm, :load)
#         _PMD.constraint_mc_load_setpoint(pm, id)
#     end

#     for (i,bus) in _PMs.ref(pm, :bus)
#         _PMD.constraint_mc_load_power_balance(pm, i)

#         # PV Bus Constraints
#         if length(_PMs.ref(pm, :bus_gens, i)) > 0 && !(i in _PMs.ids(pm,:ref_buses))
#             # this assumes inactive generators are filtered out of bus_gens
#             @assert bus["bus_type"] == 2

#             _PMD.constraint_mc_voltage_magnitude_only(pm, i)
#             for j in _PMs.ref(pm, :bus_gens, i)
#                 _PMD.constraint_mc_gen_power_setpoint_real(pm, j)
#             end
#         end
#     end

#     for i in _PMs.ids(pm, :branch)
#         _PMD.constraint_mc_power_losses(pm, i)
#         _PMD.constraint_mc_model_voltage_magnitude_difference(pm, i)
#         _PMD.constraint_mc_voltage_angle_difference(pm, i)

#         _PMD.constraint_mc_thermal_limit_from(pm, i)
#         _PMD.constraint_mc_thermal_limit_to(pm, i)
#     end

#     for i in _PMs.ids(pm, :transformer)
#         _PMD.constraint_mc_transformer_power(pm, i)
#     end
#     println(stop)
# end