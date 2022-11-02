
function solve_rdt(data::Dict{String,Any}, model_type, solver; kwargs...)
    #return _PMs.run_model(data, model_type, solver, build_mc_rdt; multiconductor=true, multinetwork=true, ref_extensions=[ref_add_rdt!], kwargs...)
    return _PMD.solve_mc_model(data, model_type, solver, build_mc_rdt; multinetwork=true, ref_extensions=[ref_add_rdt!], eng2math_extensions=[transform_switch_inline_ne!,transform_switch_inline!], kwargs...)
end

function build_mc_rdt(pm::_PMD.AbstractUnbalancedPowerModel)
    variable_he(pm); # 1d h_e variables
    variable_te(pm); # 1d t_e variables
    variable_xe(pm); # 1d x_e variables
    variable_ue(pm); # 1d u_e variables

    for n in _INs.nw_ids(pm, _PMD.pmd_it_sym)
        _PMD.variable_mc_bus_voltage_indicator(pm; nw=n, relax=true);
        _PMD.variable_mc_bus_voltage_on_off(pm; nw=n);  # is this necessary ?
        _PMD.variable_mc_branch_power(pm; nw=n,);
        _PMD.variable_mc_transformer_power(pm; nw=n,);
        _PMD.variable_mc_switch_power(pm;nw=n)
        _PMD.variable_mc_gen_indicator(pm; nw=n, relax=true);
        #_PMD.variable_mc_gen_power_setpoint_on_off(pm; nw=n);
        _PMD.variable_mc_load_indicator(pm; nw=n, relax=true);
        _PMD.variable_mc_shunt_indicator(pm; nw=n, relax=true);

#        variable_branch_be(pm) # b_e variables

         variable_xe_s(pm; nw=n) # x_e variables constraint 1c
#        variable_te_s(pm; nw=n) # t_e variables constraint 1c
#        variable_he_s(pm; nw=n) # h_e variables constraint 1c
         _PMD.variable_mc_switch_state(pm; nw=n) # t_e variables
#        variable_ye_s(pm; nw=n) # # y_e variables

        variable_mc_switch_inline_ne_state(pm; nw=n) # t_e variables

#        constraint_variable(pm; nw=n); # constraint 1b
#        constraint_switch(pm; nw=n); # constraint 6a
#        constraint_active_line(pm; nw=n); # constraint 6b

        _PMD.constraint_mc_model_voltage(pm; nw=n);   # Some forms of the power flow equations have special constraints to link voltages together.  Most power flow models don't use this

# need an accurate constraint for "switch_inline_ne" to interpret them correction

#        for i in _PMs.ids(pm, :bus_bal)
#            constraint_mc_vm_vuf(pm, i) # voltage imbalance constraint
#        end

        for i in _PMD.ids(pm, :ref_buses; nw=n)
            _PMD.constraint_mc_theta_ref(pm, i; nw=n) # slack bus constraint
        end

        _PMD.constraint_mc_bus_voltage_on_off(pm; nw=n) # allow voltages to go to zero if the bus is turned off.  Not in the paper but the analgoue of constraint (2e)

        for i in _PMD.ids(pm, :gen; nw=n)
#            _PMD.constraint_mc_gen_power_on_off(pm, i; nw=n)
        end

        for i in _PMD.ids(pm, :bus; nw=n)
#            _PMD.constraint_mc_shed_power_balance(pm, i; nw=n)
        end

        for i in _PMD.ids(pm, :branch; nw=n)
#            constraint_mc_ohms_yt_from(pm, i; nw=n) # defines pij on z_e
#            constraint_mc_ohms_yt_to(pm, i; nw=n) # defines pji on z_e

            _PMD.constraint_mc_voltage_angle_difference(pm, i; nw=n) # not in paper, but fine to include

            _PMD.constraint_mc_thermal_limit_from(pm, i; nw=n) # not in paper, but fine to include
            _PMD.constraint_mc_thermal_limit_to(pm, i; nw=n) # not in paper, but fine to include
#            constraint_cycle_function(pm, i; nw=n)
        end

        for i in _PMD.ids(pm, :transformer; nw=n)
            _PMD.constraint_mc_transformer_power(pm, i; nw=n) # not in paper, but fine to include
        end

#        for i in _PMs.ids(pm, :arcs_bal)
#            constraint_branch_be(pm, i); # constraint 2f & 2g
#            constraint_balance_flow(pm, i); # constraint 3a & 3b
#        end

        # constraint_critical_load(pm) # constraint 4c
        # constraint_non_critical_load(pm) # constraint

        # cycle elimination constraints
#        for tour in _PMs.ids(pm, :arc_tour)
#            constraint_cycle_elimination(pm, n, tour)
#        end

        for i in _PMD.ids(pm, :switch; nw=n)
            _PMD.constraint_mc_switch_state(pm, i; nw=n)
            _PMD.constraint_mc_switch_thermal_limit(pm, i; nw=n)
            _PMD.constraint_mc_switch_ampacity(pm, i; nw=n)
        end
    end


#    objective(pm)
end
