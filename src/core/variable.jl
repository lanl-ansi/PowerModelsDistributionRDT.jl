
#function variable_branch_be(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
#    variable_branch_be_p(pm)
#    variable_branch_be_q(pm)
#end

#function variable_branch_be_p(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
#    _PMs.var(pm, nw)[:be_p] = JuMP.@variable(pm.model,
#    [(l,i,j) in _PMs.ref(pm, nw, :arcs_bal)],
#    base_name = "$(nw)_branch_be_p",
#    binary = true,
#    start = 0)
#end

#function variable_branch_be_q(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
#    _PMs.var(pm, nw)[:be_q] = JuMP.@variable(pm.model,
#    [(l,i,j) in _PMs.ref(pm, nw, :arcs_bal)],
#    base_name = "$(nw)_branch_be_q",
#    binary = true,
#    start = 0)
#end

function variable_xe(pm::_PMD.AbstractUnbalancedPowerModel; nw::Int=_PMD.nw_id_default, relax::Bool=false, report::Bool=true)
    if relax
        xe = _PMD.var(pm, nw)[:xe] = JuMP.@variable(pm.model,
                [i in _PMD.ids(pm, nw, :branch_ne)],
                base_name="$(nw)_xe",
                lower_bound = 0,
                upper_bound = 1,
                start=_PMD.comp_start_value(_PMD.ref(pm, nw, :branch_ne, i), "xe_start", i, 0.0)
             )
    else
        xe = _PMD.var(pm, nw)[:xe] = JuMP.@variable(pm.model,
                [i in _PMD.ids(pm, nw, :branch_ne)],
                base_name="$(nw)_xe",
                binary = true,
                start=_PMD.comp_start_value(_PMD.ref(pm, nw, :branch_ne, i), "xe_start", i, 0.0)
             )
     end
    report && _INs.sol_component_value(pm, _PMD.pmd_it_sym, nw, :branch_ne, :xe, _PMD.ids(pm, nw, :branch_ne), xe)
end

function variable_ue(pm::_PMD.AbstractUnbalancedPowerModel; nw::Int=_PMD.nw_id_default, relax::Bool=false, report::Bool=true)
    if relax
        ue = _PMD.var(pm, nw)[:ue] = JuMP.@variable(pm.model,
                [i in _PMD.ids(pm, nw, :gen_ne)],
                base_name="$(nw)_ue",
                lower_bound = 0,
                upper_bound = 1,
                start=_PMD.comp_start_value(_PMD.ref(pm, nw, :gen_ne, i), "ue_start", i, 0.0)
             )
    else
        ue = _PMD.var(pm, nw)[:ue] = JuMP.@variable(pm.model,
                [i in _PMD.ids(pm, nw, :gen_ne)],
                base_name="$(nw)_ue",
                binary = true,
                start=_PMD.comp_start_value(_PMD.ref(pm, nw, :gen_ne, i), "ue_start", i, 0.0)
             )
     end
    report && _INs.sol_component_value(pm, _PMD.pmd_it_sym, nw, :gen_ne, :ue, _PMD.ids(pm, nw, :gen_ne), ue)
end

function variable_he(pm::_PMD.AbstractUnbalancedPowerModel; nw::Int=_PMD.nw_id_default, relax::Bool=false, report::Bool=true)
    if relax
        he = _PMD.var(pm, nw)[:he] = JuMP.@variable(pm.model,
                [i in _PMD.ref(pm, nw, :branch_harden)],
                base_name="$(nw)_he",
                lower_bound = 0,
                upper_bound = 1,
                start=_PMD.comp_start_value(_PMD.ref(pm, nw, :branch, i), "he_start", i, 0.0)
             )
    else
        he = _PMD.var(pm, nw)[:he] = JuMP.@variable(pm.model,
                [i in _PMD.ref(pm, nw, :branch_harden)],
                base_name="$(nw)_he",
                binary = true,
                start=_PMD.comp_start_value(_PMD.ref(pm, nw, :branch, i), "he_start", i, 0.0)
             )
     end

     report && _INs.sol_component_value(pm, _PMD.pmd_it_sym, nw, :branch, :he, _PMD.ref(pm, nw, :branch_harden), he)
end

function variable_te(pm::_PMD.AbstractUnbalancedPowerModel; nw::Int=_PMD.nw_id_default, relax::Bool=false, report::Bool=true)
    if relax
        te = _PMD.var(pm, nw)[:te] = JuMP.@variable(pm.model,
                [i in _PMD.ids(pm, nw, :switch_inline_ne)],
                base_name="$(nw)_te",
                lower_bound = 0,
                upper_bound = 1,
                start=_PMD.comp_start_value(_PMD.ref(pm, nw, :switch_inline_ne, i), "te_start", i, 0.0)
             )
    else
        te = _PMD.var(pm, nw)[:te] = JuMP.@variable(pm.model,
                [i in _PMD.ids(pm, nw, :switch_inline_ne)],
                base_name="$(nw)_te",
                binary = true,
                start=_PMD.comp_start_value(_PMD.ref(pm, nw, :switch_inline_ne, i), "te_start", i, 0.0)
             )
     end

     report && _INs.sol_component_value(pm, _PMD.pmd_it_sym, nw, :switch_inline_ne, :te, _PMD.ids(pm, nw, :switch_inline_ne), te)
end

function variable_xe_s(pm::_PMD.AbstractUnbalancedPowerModel; nw::Int=_PMD.nw_id_default, relax::Bool=false, report::Bool=true)
    if relax
        xe_s = _PMD.var(pm, nw)[:xe_s] = JuMP.@variable(pm.model,
                [i in _PMD.ids(pm, nw, :branch_ne)],
                base_name="$(nw)_xe_s_ne",
                lower_bound = 0,
                upper_bound = 1,
                start=_PMD.comp_start_value(_PMD.ref(pm, nw, :branch_ne, i), "xe_start", i, 0.0)
              )
    else
        xe_s = _PMD.var(pm, nw)[:xe_s] = JuMP.@variable(pm.model,
                     [i in _PMD.ids(pm, nw, :branch_ne)],
                     base_name="$(nw)_xe_s_ne",
                     binary = true,
                     start=_PMD.comp_start_value(_PMD.ref(pm, nw, :branch_ne, i), "xe_start", i, 0.0)
            )
     end
    report && _INs.sol_component_value(pm, _PMD.pmd_it_sym, nw, :branch_ne, :xe_s, _PMD.ids(pm, nw, :branch_ne), xe_s)
end

function variable_ze_s(pm::_PMD.AbstractUnbalancedPowerModel; nw::Int=_PMD.nw_id_default, relax::Bool=false, report::Bool=true)
    if relax
        ze_s = _PMD.var(pm, nw)[:ze_s] = JuMP.@variable(pm.model,
                [i in _PMD.ids(pm, nw, :branch)],
                base_name="$(nw)_ze_s",
                lower_bound = 0,
                upper_bound = 1,
                start=_PMD.comp_start_value(_PMD.ref(pm, nw, :branch, i), "ze_start", i, 0.0)
             )
        ze_s_xfr = _PMD.var(pm, nw)[:ze_s_xfr] = JuMP.@variable(pm.model,
                [i in _PMD.ids(pm, nw, :transformer)],
                base_name="$(nw)_ze_s_xfr",
                lower_bound = 0,
                upper_bound = 1,
                start=_PMD.comp_start_value(_PMD.ref(pm, nw, :transformer, i), "ze_start", i, 0.0)
              )
    else
        ze_s = _PMD.var(pm, nw)[:ze_s] = JuMP.@variable(pm.model,
                [i in _PMD.ids(pm, nw, :branch)],
                base_name="$(nw)_ze_s",
                binary = true,
                start=_PMD.comp_start_value(_PMD.ref(pm, nw, :branch, i), "ze_start", i, 0.0)
             )
        ze_s_xfr = _PMD.var(pm, nw)[:ze_s_xfr] = JuMP.@variable(pm.model,
                    [i in _PMD.ids(pm, nw, :transformer)],
                    base_name="$(nw)_ze_s_xfr",
                    binary = true,
                    start=_PMD.comp_start_value(_PMD.ref(pm, nw, :transformer, i), "ze_start", i, 0.0)
            )
     end
    report && _INs.sol_component_value(pm, _PMD.pmd_it_sym, nw, :branch, :ze_s, _PMD.ids(pm, nw, :branch), ze_s)
    report && _INs.sol_component_value(pm, _PMD.pmd_it_sym, nw, :transformer, :ze_s_xfr, _PMD.ids(pm, nw, :transformer), ze_s_xfr)
end


function variable_he_s(pm::_PMD.AbstractUnbalancedPowerModel; nw::Int=_PMD.nw_id_default, relax::Bool=false, report::Bool=true)
    if relax
        he_s = _PMD.var(pm, nw)[:he_s] = JuMP.@variable(pm.model,
                [i in _PMD.ref(pm, nw, :branch_harden)],
                base_name="$(nw)_he_s",
                lower_bound = 0,
                upper_bound = 1,
                start=_PMD.comp_start_value(_PMD.ref(pm, nw, :branch, i), "he_start", i, 0.0)
             )
    else
        he_s = _PMD.var(pm, nw)[:he_s] = JuMP.@variable(pm.model,
                [i in _PMD.ref(pm, nw, :branch_harden)],
                base_name="$(nw)_he_s",
                binary = true,
                start=_PMD.comp_start_value(_PMD.ref(pm, nw, :branch, i), "he_start", i, 0.0)
             )
     end
    report && _INs.sol_component_value(pm, _PMD.pmd_it_sym, nw, :branch, :he_s, _PMD.ref(pm, nw, :branch_harden), he_s)
end

function variable_ue_s(pm::_PMD.AbstractUnbalancedPowerModel; nw::Int=_PMD.nw_id_default, relax::Bool=false, report::Bool=true)
    if relax
        ue_s = _PMD.var(pm, nw)[:ue_s] = JuMP.@variable(pm.model,
                [i in _PMD.ids(pm, nw, :gen_ne)],
                base_name="$(nw)_ue_s",
                lower_bound = 0,
                upper_bound = 1,
                start=_PMD.comp_start_value(_PMD.ref(pm, nw, :gen_ne, i), "ue_start", i, 0.0)
             )
    else
        ue_s = _PMD.var(pm, nw)[:ue_s] = JuMP.@variable(pm.model,
                [i in _PMD.ids(pm, nw, :gen_ne)],
                base_name="$(nw)_ue_s",
                binary = true,
                start=_PMD.comp_start_value(_PMD.ref(pm, nw, :gen_ne, i), "ue_start", i, 0.0)
             )
     end
    report && _INs.sol_component_value(pm, _PMD.pmd_it_sym, nw, :gen_ne, :ue_s, _PMD.ids(pm, nw, :gen_ne), ue_s)
end


#function variable_he_s(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
#    _PMs.var(pm, nw)[:he_d_s] = JuMP.@variable(pm.model,
#    [(l,i,j) in _PMs.ref(pm, nw, :arcs_damaged)],
#    base_name = "$(nw)_branch_he_damaged_s",
#    binary = true,
#    start = 0)
#    _PMs.var(pm, nw)[:he_n_s] = JuMP.@variable(pm.model,
#    [(l,i,j) in _PMs.ref(pm, nw, :arcs_new)],
#    base_name = "$(nw)_branch_he_new_s",
#    binary = true,
#    start = 0)
#end

#function variable_ye_s(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
#    _PMs.var(pm, nw)[:ye_s] = JuMP.@variable(pm.model,
#    [(l,i,j) in _PMs.ref(pm, nw, :arcs)],
#    base_name = "$(nw)_branch_ye_s",
#    binary = true,
#    start = 0)
#end


"switch_inline_ne state (open/close) variables"
function variable_mc_switch_inline_ne_state(pm::_PMD.AbstractUnbalancedPowerModel; nw::Int=nw_id_default, report::Bool=true, relax::Bool=false)
    if relax
        state = _PMD.var(pm, nw)[:switch_inline_ne_state] = JuMP.@variable(
            pm.model,
            [l in _PMD.ids(pm, nw, :switch_inline_ne)],
            base_name="$(nw)_switch_inline_ne_state_$(l)",
            lower_bound = 0,
            upper_bound = 1,
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :switch_inline_ne, l), "state_start", 0)
        )
    else
        state = _PMD.var(pm, nw)[:switch_inline_ne_state] = JuMP.@variable(
            pm.model,
            [l in _PMD.ids(pm, nw, :switch_inline_ne)],
            base_name="$(nw)_switch_inline_ne_state_$(l)",
            binary = true,
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :switch_inline_ne, l), "state_start", 0)
        )
    end

    report && _INs.sol_component_value(pm, _PMD.pmd_it_sym, nw, :switch_inline_ne, :switch_inline_ne_state, _PMD.ids(pm, nw, :switch_inline_ne), state)
end
