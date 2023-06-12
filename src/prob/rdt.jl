
function solve_rdt(data::Dict{String,Any}, model_type, solver; kwargs...)
    return _PMD.solve_mc_model(data, model_type, solver, build_mc_rdt; multinetwork=true, ref_extensions=[ref_add_rdt!], eng2math_extensions=[transform_switch_inline_ne!,transform_switch_inline!], kwargs...)
end

function build_mc_rdt(pm::_PMD.AbstractUBFModels)
    variable_he(pm); # 1d h_e variables
    variable_te(pm; relax=true); # 1d t_e variables - can be continous because the combination of the objective and constraint 1b will force them to 0 or 1
    variable_xe(pm); # 1d x_e variables
    variable_ue(pm); # 1d u_e variables

    for n in _IM.nw_ids(pm, _PMD.pmd_it_sym)
        _PMD.variable_mc_bus_voltage_indicator(pm; nw=n, relax=true);
        _PMD.variable_mc_bus_voltage_on_off(pm; nw=n);  # is this necessary ?
        _PMD.variable_mc_branch_power(pm; nw=n);
         variable_mc_branch_ne_power(pm; nw=n);
        _PMD.variable_mc_transformer_power(pm; nw=n);
        variable_mc_transformer_ne_power(pm; nw=n)
        _PMD.variable_mc_switch_power(pm;nw=n)
        variable_mc_switch_inline_ne_power(pm;nw=n)
        _PMD.variable_mc_gen_indicator(pm; nw=n, relax=true);
        variable_mc_gen_ne_indicator(pm; nw=n, relax=true);
        _PMD.variable_mc_generator_power_on_off(pm; nw=n);
        variable_mc_generator_ne_power_on_off(pm; nw=n);
        _PMD.variable_mc_load_indicator(pm; nw=n, relax=true);
        _PMD.variable_mc_shunt_indicator(pm; nw=n, relax=true);

#        variable_branch_be(pm) # b_e variables

         variable_xe_s(pm; nw=n, relax=true) # x_e variables - can be continous because the combination of the objective and constraint 1b will force them to 0 or 1
#         variable_ze_s(pm; nw=n) # z_e variables - may not need because switches are distinct objects now
         variable_he_s(pm; nw=n, relax=true) # h_e variables - can be continous because the combination of the objective and constraint 1b will force them to 0 or 1
         variable_ue_s(pm; nw=n, relax=true) # u_e variables - can be continous because the combination of the objective and constraint 1b will force them to 0 or 1
         _PMD.variable_mc_switch_state(pm; nw=n) # t_e variables
         variable_mc_switch_inline_ne_state(pm; nw=n) # t_e variables
#        variable_ye_s(pm; nw=n) # # y_e variables

        for i in _PMD.ids(pm, n, :gen_ne)
            constraint_ue(pm, i; nw=n) # constraint 1b for u variables
        end

        for i in _PMD.ids(pm, n, :branch_ne)
            constraint_xe(pm, i; nw=n) # constraint 1b for x variables
        end

        for i in _PMD.ids(pm, n, :switch_inline_ne)
            constraint_te(pm, i; nw=n) # constraint 1b for t variables
        end

        for i in _PMD.ids(pm, n, :branch_harden)
            constraint_he(pm, i; nw=n) # constraint 1b for h variables
        end

#        constraint_switch(pm; nw=n); # constraint 6a
#        constraint_active_line(pm; nw=n); # constraint 6b

        _PMD.constraint_mc_model_voltage(pm; nw=n);   # Some forms of the power flow equations have special constraints to link voltages together.  Most power flow models don't use this


#        for i in _PMs.ids(pm, :bus_bal)
#            constraint_mc_vm_vuf(pm, i) # voltage imbalance constraint
#        end

        for i in _PMD.ids(pm, n, :ref_buses)
            _PMD.constraint_mc_theta_ref(pm, i; nw=n) # slack bus constraint
        end

#        _PMD.constraint_mc_bus_voltage_on_off(pm; nw=n) # allow voltages to go to zero if the bus is turned off.  Not in the paper but the analgoue of constraint (2e)

        for i in _PMD.ids(pm, n, :gen)
#            _PMD.constraint_mc_gen_power_on_off(pm, i; nw=n)
        end

        for id in _PMD.ids(pm, n, :load)
#            constraint_mc_load_power(pm, id)
        end

        for i in _PMD.ids(pm, n, :gen_ne)
        #  some stuff
        end

        for i in _PMD.ids(pm, n, :bus)
            constraint_mc_power_balance_shed_ne(pm, i; nw=n) # constraint 2c
        end

        for i in _PMD.ref(pm, n, :undamaged_branch)
            _PMD.constraint_mc_power_losses(pm, i; nw=n)
            _PMD.constraint_mc_model_voltage_magnitude_difference(pm, i, nw=n)

            _PMD.constraint_mc_voltage_angle_difference(pm, i; nw=n) # not in paper, but fine to include

            _PMD.constraint_mc_thermal_limit_from(pm, i; nw=n) # constraint 2d
            _PMD.constraint_mc_thermal_limit_to(pm, i; nw=n) # constraint 2d

            _PMD.constraint_mc_ampacity_from(pm, i; nw=n) # not in paper, but fine to include
            _PMD.constraint_mc_ampacity_to(pm, i; nw=n) # not in paper, but fine to include

#            constraint_cycle_function(pm, i; nw=n)
            # constraint 2b is implict
        end

        for i in _PMD.ref(pm, n, :damaged_branch) # need to break this out into damaged and un damaged branches
            constraint_mc_thermal_limit_from_damaged(pm, i; nw=n) # constraint 2d
            constraint_mc_thermal_limit_to_damaged(pm, i; nw=n) # constraint 2d

            constraint_mc_ampacity_from_damaged(pm, i; nw=n) # not in paper, but fine to include
            constraint_mc_ampacity_to_damaged(pm, i; nw=n) # not in paper, but fine to include


            # constraint 2b is implict
        end

        for i in _PMD.ids(pm, n, :branch_ne)
            constraint_mc_thermal_limit_from_ne(pm, i; nw=n) # constraint 2d
            constraint_mc_thermal_limit_to_ne(pm, i; nw=n) # constraint 2d

            ### some stuff
        end

        for i in _PMD.ids(pm, n, :transformer)
            _PMD.constraint_mc_transformer_power(pm, i; nw=n) # not in paper, but fine to include
        end

#        for i in _PMs.ids(pm, :arcs_bal)
#            constraint_branch_be(pm, i); # constraint 2f & 2g
#            constraint_balance_flow(pm, i); # constraint 3a & 3b
#        end

        constraint_critical_load(pm; nw=n) # constraint 4c
        constraint_total_load(pm; nw=n) # anaologue to constraint 4c

        # cycle elimination constraints
#        for tour in _PMs.ids(pm, :arc_tour)
#            constraint_cycle_elimination(pm, n, tour)
#        end

        for i in _PMD.ids(pm, n, :switch)
            _PMD.constraint_mc_switch_state(pm, i; nw=n)
            _PMD.constraint_mc_switch_thermal_limit(pm, i; nw=n)
            _PMD.constraint_mc_switch_ampacity(pm, i; nw=n)
        end

        for i in _PMD.ids(pm, n, :switch_inline_ne)
            ### some stuff
        end


#        for i in ids(pm, n, :storage)
#               _PMD.constraint_storage_state(pm, i; nw=n)
#               _PMD.constraint_storage_complementarity_nl(pm, i; nw=n)
#               _PMD.constraint_mc_storage_losses(pm, i; nw=n)
#               _PMD.constraint_mc_storage_thermal_limit(pm, i; nw=n)
#           end

    end

    objective_rdt(pm)
end
