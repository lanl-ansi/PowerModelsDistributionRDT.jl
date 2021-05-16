
function variable_branch_be(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
    variable_branch_be_p(pm)
    variable_branch_be_q(pm)
end

function variable_branch_be_p(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
    _PMs.var(pm, nw)[:be_p] = JuMP.@variable(pm.model,
    [(l,i,j) in _PMs.ref(pm, nw, :arcs_bal)],
    base_name = "$(nw)_branch_be_p",
    binary = true,
    start = 0)
end

function variable_branch_be_q(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
    _PMs.var(pm, nw)[:be_q] = JuMP.@variable(pm.model,
    [(l,i,j) in _PMs.ref(pm, nw, :arcs_bal)],
    base_name = "$(nw)_branch_be_q",
    binary = true,
    start = 0)
end

function variable_xe(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
    _PMs.var(pm, nw)[:xe_d] = JuMP.@variable(pm.model,
    [(l,i,j) in pm.ref[:arcs_damaged_all]],
    base_name = "$(nw)_branch_xe_damaged",
    binary = true,
    start = 0)
    _PMs.var(pm, nw)[:xe_n] = JuMP.@variable(pm.model,
    [(l,i,j) in pm.ref[:arcs_new_all]],
    base_name = "$(nw)_branch_xe_new",
    binary = true,
    start = 0)
    # println(pm.model)
    # println(stop)
end

function variable_te(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
    _PMs.var(pm, nw)[:te_d] = JuMP.@variable(pm.model,
    [(l,i,j) in pm.ref[:arcs_damaged_all]],
    base_name = "$(nw)_branch_te_damaged",
    binary = true,
    start = 0)
    _PMs.var(pm, nw)[:te_n] = JuMP.@variable(pm.model,
    [(l,i,j) in pm.ref[:arcs_new_all]],
    base_name = "$(nw)_branch_te_new",
    binary = true,
    start = 0)
end

function variable_he(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
    _PMs.var(pm, nw)[:he_d] = JuMP.@variable(pm.model,
    [(l,i,j) in pm.ref[:arcs_damaged_all]],
    base_name = "$(nw)_branch_he_damaged",
    binary = true,
    start = 0)
    _PMs.var(pm, nw)[:he_n] = JuMP.@variable(pm.model,
    [(l,i,j) in pm.ref[:arcs_new_all]],
    base_name = "$(nw)_branch_he_new",
    binary = true,
    start = 0)
end

function variable_z_branch(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
    _PMs.var(pm, nw)[:z_branch] = JuMP.@variable(pm.model,
    [(l,i,j) in _PMs.ref(pm, nw, :arcs)],
    lower_bound = 0,
    upper_bound = 1,
    base_name = "$(nw)_branch_z",
    start = 0)
end

function variable_xe_s(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
    _PMs.var(pm, nw)[:xe_d_s] = JuMP.@variable(pm.model,
    [(l,i,j) in _PMs.ref(pm, nw, :arcs_damaged)],
    base_name = "$(nw)_branch_xe_damaged_s",
    binary = true,
    start = 0)
    _PMs.var(pm, nw)[:xe_n_s] = JuMP.@variable(pm.model,
    [(l,i,j) in _PMs.ref(pm, nw, :arcs_new)],
    base_name = "$(nw)_branch_xe_new_s",
    binary = true,
    start = 0)
end

function variable_te_s(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
    _PMs.var(pm, nw)[:te_d_s] = JuMP.@variable(pm.model,
    [(l,i,j) in _PMs.ref(pm, nw, :arcs_damaged)],
    base_name = "$(nw)_branch_te_damaged_s",
    binary = true,
    start = 0)
    _PMs.var(pm, nw)[:te_n_s] = JuMP.@variable(pm.model,
    [(l,i,j) in _PMs.ref(pm, nw, :arcs_new)],
    base_name = "$(nw)_branch_te_new_s",
    binary = true,
    start = 0)
end

function variable_he_s(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
    _PMs.var(pm, nw)[:he_d_s] = JuMP.@variable(pm.model,
    [(l,i,j) in _PMs.ref(pm, nw, :arcs_damaged)],
    base_name = "$(nw)_branch_he_damaged_s",
    binary = true,
    start = 0)
    _PMs.var(pm, nw)[:he_n_s] = JuMP.@variable(pm.model,
    [(l,i,j) in _PMs.ref(pm, nw, :arcs_new)],
    base_name = "$(nw)_branch_he_new_s",
    binary = true,
    start = 0)
end

function variable_ye_s(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
    _PMs.var(pm, nw)[:ye_s] = JuMP.@variable(pm.model,
    [(l,i,j) in _PMs.ref(pm, nw, :arcs)],
    base_name = "$(nw)_branch_ye_s",
    binary = true,
    start = 0)
end
