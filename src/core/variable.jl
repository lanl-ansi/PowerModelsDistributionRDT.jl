
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


#function variable_z_branch(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
#    _PMs.var(pm, nw)[:z_branch] = JuMP.@variable(pm.model,
#    [(l,i,j) in _PMs.ref(pm, nw, :arcs)],
#    lower_bound = 0,
#    upper_bound = 1,
#    base_name = "$(nw)_branch_z",
#    start = 0)
#end

#function variable_xe_s(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
#    _PMs.var(pm, nw)[:xe_d_s] = JuMP.@variable(pm.model,
#    [(l,i,j) in _PMs.ref(pm, nw, :arcs_damaged)],
#    base_name = "$(nw)_branch_xe_damaged_s",
#    binary = true,
#    start = 0)
#    _PMs.var(pm, nw)[:xe_n_s] = JuMP.@variable(pm.model,
#    [(l,i,j) in _PMs.ref(pm, nw, :arcs_new)],
#    base_name = "$(nw)_branch_xe_new_s",
#    binary = true,
#    start = 0)
#end

#function variable_te_s(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
#    _PMs.var(pm, nw)[:te_d_s] = JuMP.@variable(pm.model,
#    [(l,i,j) in _PMs.ref(pm, nw, :arcs_damaged)],
#    base_name = "$(nw)_branch_te_damaged_s",
#    binary = true,
#    start = 0)
#    _PMs.var(pm, nw)[:te_n_s] = JuMP.@variable(pm.model,
#    [(l,i,j) in _PMs.ref(pm, nw, :arcs_new)],
#    base_name = "$(nw)_branch_te_new_s",
#    binary = true,
#    start = 0)
#end

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
