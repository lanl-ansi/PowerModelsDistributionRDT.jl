
function constraint_branch_be(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw) 
    constraint_branch_be_p(pm, i, nw)
    constraint_branch_be_q(pm, i, nw)
end


function constraint_branch_be_p(pm::_PMs.AbstractPowerModel, i::Int, nw::Int) # check trans arcs and branch arches
    arcs = _PMs.ref(pm, nw, :arcs_bal)[i]
    arcs in _PMs.ref(pm, nw, :arcs_trans) ? trans = true : trans = false
    constraint_branch_be_p(pm, nw, arcs, trans)
end 

function constraint_branch_be_q(pm::_PMs.AbstractPowerModel, i::Int, nw::Int) 
    arcs = _PMs.ref(pm, nw, :arcs_bal)[i] 
    arcs in _PMs.ref(pm, nw, :arcs_trans) ? trans = true : trans = false
    constraint_branch_be_q(pm, nw, arcs, trans)
end 

function constraint_balance_flow(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw) 
    constraint_balance_p_flow(pm, i, nw)
    constraint_balance_q_flow(pm, i, nw)
end

function constraint_balance_p_flow(pm::_PMs.AbstractPowerModel, i::Int, nw::Int)
    arcs = _PMs.ref(pm, nw, :arcs_bal)[i]
    arcs in _PMs.ref(pm, nw, :arcs_trans) ? trans = true : trans = false
    constraint_balance_p_flow(pm, nw, arcs, trans)
end

function constraint_balance_q_flow(pm::_PMs.AbstractPowerModel, i::Int, nw::Int)
    arcs = _PMs.ref(pm, nw, :arcs_bal)[i]
    arcs in _PMs.ref(pm, nw, :arcs_trans) ? trans = true : trans = false
    constraint_balance_q_flow(pm, nw, arcs, trans)
end


function constraint_critical_load(pm::_PMs.AbstractPowerModel, nw::Int)
    haskey(_PMs.ref(pm, nw), :critical_level) ? limit = _PMs.ref(pm, nw)[:critical_level] : limit = 0.0
    bus = zeros(Int, 0)
    for i in _PMs.ids(pm, :load)
        _PMs.ref(pm, nw, :load, i)["weight"] == 100 ? append!(bus, i) : nothing 
    end 
    constraint_critical_load(pm, nw, bus, limit)
end

function constraint_non_critical_load(pm::_PMs.AbstractPowerModel, nw::Int)
    haskey(_PMs.ref(pm, nw), :demand_level) ? limit = _PMs.ref(pm, nw)[:demand_level] : limit = 1.0
    constraint_non_critical_load(pm, nw, limit)
end

function constraint_mc_vm_vuf(pm::_PMs.AbstractPowerModel, bus_id::Int; nw::Int=pm.cnw)
    bus = _PMs.ref(pm, nw, :bus_bal)[bus_id]
    vufmax = _PMs.ref(pm, nw, :bus)[bus]["vm_vuf_max"]
    constraint_mc_vm_vuf(pm, nw, bus, vufmax)
end

function constraint_mc_vm_vuf_mod(pm::_PMs.AbstractPowerModel, bus_id::Int; nw::Int=pm.cnw)
    bus = _PMs.ref(pm, nw, :bus_bal)[bus_id]
    vufmax = _PMs.ref(pm, nw, :bus)[bus]["vm_vuf_max"]
    constraint_mc_vm_vuf(pm, nw, bus, vufmax)
end

function constraint_switch(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
    constraint_switch_damaged(pm, nw)
    constraint_switch_new(pm, nw)
end

function constraint_switch_damaged(pm::_PMs.AbstractPowerModel, nw::Int)
    for (l,i,j) in _PMs.ref(pm, nw, :arcs_damaged) 
        xe = _PMs.var(pm, nw, :xe_d_s, (l,i,j))
        te = _PMs.var(pm, nw, :te_d_s, (l,i,j))
        constraint_switch(pm, xe, te)
    end
end

function constraint_switch_new(pm::_PMs.AbstractPowerModel, nw::Int)
    for (l,i,j) in _PMs.ref(pm, nw, :arcs_new)
        xe = _PMs.var(pm, nw, :xe_n_s, (l,i,j))
        te = _PMs.var(pm, nw, :te_n_s, (l,i,j))
        constraint_switch(pm, xe, te)
    end
end

function constraint_harden(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
    constraint_harden_damaged(pm, nw)
    constraint_harden_new(pm, nw)
end

function constraint_harden_damaged(pm::_PMs.AbstractPowerModel, nw::Int)
    for (l,i,j) in _PMs.ref(pm, nw, :arcs_damaged) 
        xe = _PMs.var(pm, nw, :xe_d_s, (l,i,j))
        he = _PMs.var(pm, nw, :he_d_s, (l,i,j))
        constraint_harden(pm, xe, he)
    end
end

function constraint_harden_new(pm::_PMs.AbstractPowerModel, nw::Int)
    for (l,i,j) in _PMs.ref(pm, nw, :arcs_new) 
        xe = _PMs.var(pm, nw, :xe_n_s, (l,i,j))
        he = _PMs.var(pm, nw, :he_n_s, (l,i,j))
        constraint_harden(pm, xe, he)
    end
end

function constraint_active_line(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
    constraint_active_line_existing(pm, nw)
    constraint_active_line_new(pm, nw)
    constraint_active_line_active(pm, nw)
end

function constraint_active_line_existing(pm::_PMs.AbstractPowerModel, nw::Int)
    for (l,i,j) in _PMs.ref(pm, nw, :arcs_damaged) 
        constraint_activation_damage(pm, nw, (l,i,j))
    end
end

function constraint_active_line_new(pm::_PMs.AbstractPowerModel, nw::Int)
    for (l,i,j) in _PMs.ref(pm, nw, :arcs_new)
        constraint_activation_new(pm, nw, (l,i,j))
    end
end

function constraint_active_line_active(pm::_PMs.AbstractPowerModel, nw::Int)
    for (l,i,j) in _PMs.ref(pm, nw, :arcs)
        (l,i,j) in _PMs.ref(pm, nw, :arcs_damaged) || (l,i,j) in _PMs.ref(pm, nw, :arcs_new) ? nothing : constraint_activation_active(pm, nw, (l,i,j))
    end
end

function constraint_mc_ohms_yt_from(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PMs.calc_branch_y(branch)
    tr, ti = _PMs.calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]

    constraint_mc_ohms_yt_from(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end


"ohms constraint for branches on the to-side"
function constraint_mc_ohms_yt_to(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PMs.calc_branch_y(branch)
    tr, ti = _PMs.calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    constraint_mc_ohms_yt_to(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end

function constraint_variable(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
    constraint_variable_he(pm, nw)
    constraint_variable_te(pm, nw)
    constraint_variable_xe(pm, nw)
end

function constraint_variable_he(pm::_PMs.AbstractPowerModel, nw::Int=pm.cnw)
    for (l,i,j) in _PMs.ref(pm, nw, :arcs_damaged)
        constraint_variable_he_d(pm, nw, (l,i,j))
    end
    for (l,i,j) in _PMs.ref(pm, nw, :arcs_new)
        constraint_variable_he_n(pm, nw, (l,i,j))
    end
end

function constraint_variable_te(pm::_PMs.AbstractPowerModel, nw::Int=pm.cnw)
    for (l,i,j) in _PMs.ref(pm, nw, :arcs_damaged)
        constraint_variable_te_d(pm, nw, (l,i,j))
    end
    for (l,i,j) in _PMs.ref(pm, nw, :arcs_new)
        constraint_variable_te_n(pm, nw, (l,i,j))
    end
end

function constraint_variable_xe(pm::_PMs.AbstractPowerModel, nw::Int=pm.cnw)
    for (l,i,j) in _PMs.ref(pm, nw, :arcs_damaged)
        constraint_variable_xe_d(pm, nw, (l,i,j))
    end
    for (l,i,j) in _PMs.ref(pm, nw, :arcs_new)
        constraint_variable_xe_n(pm, nw, (l,i,j))
    end
end

function constraint_mc_voltage_angle_difference(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    pair = (f_bus, t_bus)

    constraint_mc_voltage_angle_difference(pm, nw, f_idx, branch["angmin"], branch["angmax"])
end

function constraint_cycle_function(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    arcs = (i, f_bus, t_bus)

    constraint_cycle_function(pm, nw, arcs)
end

function constraint_cycle_elimination(pm, nw::Int, tour::Int)
    tours = _PMs.ref(pm, nw, :arc_tour, tour)
    
    constraint_cycle_elimination(pm, nw::Int, tours)
end