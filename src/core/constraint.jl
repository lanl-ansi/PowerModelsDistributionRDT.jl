
function constraint_ue(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, base_nw::Int, gen::Int)
    ue = _PMD.var(pm, base_nw, :ue, gen)
    ue_s = _PMD.var(pm, nw, :ue_s, gen)

    JuMP.@constraint(pm.model, ue >= ue_s)
end

function constraint_xe(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, base_nw::Int, branch::Int)
    xe = _PMD.var(pm, base_nw, :xe, branch)
    xe_s = _PMD.var(pm, nw, :xe_s, branch)

    JuMP.@constraint(pm.model, xe >= xe_s)
end

function constraint_te(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, base_nw::Int, switch::Int)
    te = _PMD.var(pm, base_nw, :te, switch)
    te_s = _PMD.var(pm, nw, :switch_inline_ne_state, switch)

    JuMP.@constraint(pm.model, te >= 1-te_s)
end

function constraint_he(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, base_nw::Int, branch::Int)
    he = _PMD.var(pm, base_nw, :he, branch)
    he_s = _PMD.var(pm, nw, :he_s, branch)

    JuMP.@constraint(pm.model, he >= he_s)
end



#function constraint_branch_be_p(pm::_PMs.AbstractPowerModel, nw::Int, arcs, trans) # check trans arcs and branch arches
#    be = _PMs.var(pm, nw, :be_p, (arcs))
#    for cnd in _PMs.conductor_ids(pm)
#        if trans
#            haskey(_PMs.ref(pm, nw, :transformer, arcs[1]), "rate_a") ? t = _PMs.ref(pm, nw, :transformer, arcs[1], "rate_a")[cnd] : t = 1
#            p = _PMs.var(pm, nw, :pt, (arcs))[cnd]
#        else
#            haskey(_PMs.ref(pm, nw, :branch, arcs[1]), "rate_a") ? t = _PMs.ref(pm, nw, :transformer, arcs[1], "rate_a")[cnd] : t = 1
#            p = _PMs.var(pm, nw, :p, (arcs))[cnd]
#        end
#        JuMP.@constraint(pm.model, p <= t * (1 - be))
#        JuMP.@constraint(pm.model, p >= -t * be)
#    end
#end

#function constraint_branch_be_q(pm::_PMs.AbstractPowerModel, nw::Int, arcs, trans)
#    be = _PMs.var(pm, nw, :be_q, (arcs))
#    for cnd in _PMs.conductor_ids(pm)
#        if trans
#            haskey(_PMs.ref(pm, nw, :transformer, arcs[1]), "rate_a") ? t = _PMs.ref(pm, nw, :transformer, arcs[1], "rate_a")[cnd] : t = 1
#            q = _PMs.var(pm, nw, :qt, (arcs))[cnd]
#        else
#            haskey(_PMs.ref(pm, nw, :branch, arcs[1]), "rate_a") ? t = _PMs.ref(pm, nw, :transformer, arcs[1], "rate_a")[cnd] : t = 1
#            q = _PMs.var(pm, nw, :q, (arcs))[cnd]
#        end
#        JuMP.@constraint(pm.model, q <= t * (1 - be))
#        JuMP.@constraint(pm.model, q >= -t * be)
#    end
#end

#function constraint_balance_p_flow(pm::_PMs.AbstractPowerModel, nw::Int, arcs, trans)
#    be = _PMs.var(pm, nw, :be_p, (arcs))
#    if trans
#        limit = _PMs.ref(pm, nw, :transformer, arcs[1])["pq_imbalance"]
#        p = [ _PMs.var(pm, nw, :pt, arcs)[c] for c in _PMs.conductor_ids(pm)]
#    else
#        limit = _PMs.ref(pm, nw, :branch, arcs[1])["pq_imbalance"]
#        p = [ _PMs.var(pm, nw, :p, arcs)[c] for c in _PMs.conductor_ids(pm)]
#    end
#    lb_beta = 1 - limit
#    ub_beta = 1 + limit
#    sum_p = JuMP.@variable(pm.model, upper_bound = 10, lower_bound = -10, start = 0)
#    JuMP.@constraint(pm.model, sum(p[c] for c in _PMs.conductor_ids(pm)) == sum_p)
#    p_be = JuMP.@variable(pm.model, upper_bound = 10, lower_bound = -10, start = 0)
#    _INs.relaxation_product(pm.model, be, sum_p, p_be)
#    for cnd in _PMs.conductor_ids(pm)
#        JuMP.@constraint(pm.model, ub_beta/3 * sum(p[c] for c in _PMs.conductor_ids(pm)) + 1/3*(lb_beta - ub_beta) * p_be >= p[cnd])
#        JuMP.@constraint(pm.model, lb_beta/3 * sum(p[c] for c in _PMs.conductor_ids(pm)) + 1/3*(ub_beta - lb_beta) * p_be <= p[cnd])
#    end
#end

#function constraint_balance_q_flow(pm::_PMs.AbstractPowerModel, nw::Int, arcs, trans)
#    be = _PMs.var(pm, nw, :be_q, (arcs))
#    if trans
#        limit = _PMs.ref(pm, nw, :transformer, arcs[1])["pq_imbalance"]
#        q = [ _PMs.var(pm, nw, :qt, arcs)[c] for c in _PMs.conductor_ids(pm)]
#    else
#        limit = _PMs.ref(pm, nw, :branch, arcs[1])["pq_imbalance"]
#        q = [ _PMs.var(pm, nw, :q, arcs)[c] for c in _PMs.conductor_ids(pm)]
#    end
#    lb_beta = 1 - limit
#    ub_beta = 1 + limit
#    sum_q = JuMP.@variable(pm.model, upper_bound = 10, lower_bound = -10, start = 0)
#    JuMP.@constraint(pm.model, sum(q[c] for c in _PMs.conductor_ids(pm)) == sum_q)
#    q_be = JuMP.@variable(pm.model, upper_bound = 10, lower_bound = -10, start = 0)
#    _INs.relaxation_product(pm.model, be, sum_q, q_be)
#    for cnd in _PMs.conductor_ids(pm)
#        JuMP.@constraint(pm.model, ub_beta/3 * sum(q[c] for c in _PMs.conductor_ids(pm)) + 1/3*(lb_beta - ub_beta) * q_be >= q[cnd])
#        JuMP.@constraint(pm.model, lb_beta/3 * sum(q[c] for c in _PMs.conductor_ids(pm)) + 1/3*(ub_beta - lb_beta) * q_be <= q[cnd])
#    end
#end

#function constraint_critical_load(pm::_PMs.AbstractPowerModel, nw::Int, bus, limit)
#    z_demand = _PMs.var(pm, nw)[:z_demand]
#    pd = sum(sum(_PMs.ref(pm, nw, :load, i, "pd")) for i in bus)
#    qd = sum(sum(_PMs.ref(pm, nw, :load, i, "qd")) for i in bus)
#    JuMP.@constraint(pm.model, limit * pd <= sum(sum(z_demand[i] * sum(_PMs.ref(pm, nw, :load, i, "pd")) for i in bus)))
#    JuMP.@constraint(pm.model, limit * qd <= sum(sum(z_demand[i] * sum(_PMs.ref(pm, nw, :load, i, "qd")) for i in bus)))
#end

#function constraint_non_critical_load(pm::_PMs.AbstractPowerModel, nw::Int, limit)
#    z_demand = _PMs.var(pm, nw)[:z_demand]
#    pd = sum(sum(_PMs.ref(pm, nw, :load, i, "pd")) for i in _PMs.ids(pm, :load))
#    qd = sum(sum(_PMs.ref(pm, nw, :load, i, "qd")) for i in _PMs.ids(pm, :load))
#    JuMP.@constraint(pm.model, limit * pd <= sum(sum(z_demand[i] * sum(_PMs.ref(pm, nw, :load, i, "pd")) for i in _PMs.ids(pm, :load))))
#    JuMP.@constraint(pm.model, limit * qd <= sum(sum(z_demand[i] * sum(_PMs.ref(pm, nw, :load, i, "qd")) for i in _PMs.ids(pm, :load))))
#end

#function constraint_mc_vm_vuf(pm::_PMs.AbstractPowerModel, nw::Int, bus, vufmax)
#    (vma, vmb, vmc) = [_PMs.var(pm, nw, :vm)[bus][i] for i in 1:3]
#    JuMP.@constraint(pm.model, vma - 1/3 * sum([_PMs.var(pm, nw, :vm)[bus][i] for i in 1:3]) <= vufmax/3 * sum([_PMs.var(pm, nw, :vm)[bus][i] for i in 1:3]))
#    JuMP.@constraint(pm.model, vmb - 1/3 * sum([_PMs.var(pm, nw, :vm)[bus][i] for i in 1:3]) <= vufmax/3 * sum([_PMs.var(pm, nw, :vm)[bus][i] for i in 1:3]))
#    JuMP.@constraint(pm.model, vmc - 1/3 * sum([_PMs.var(pm, nw, :vm)[bus][i] for i in 1:3]) <= vufmax/3 * sum([_PMs.var(pm, nw, :vm)[bus][i] for i in 1:3]))
#end

#function constraint_mc_vm_vuf_lin(pm::_PMs.AbstractPowerModel, nw::Int, bus, vufmax)
#    (wa, wb, wc) = [_PMs.var(pm, nw, :w)[bus][i] for i in 1:3]
#    JuMP.@constraint(pm.model, wa - 1/3 * sum([_PMs.var(pm, nw, :w)[bus][i] for i in 1:3]) <= vufmax/6 * sum([_PMs.var(pm, nw, :w)[bus][i] for i in 1:3]))
#    JuMP.@constraint(pm.model, wb - 1/3 * sum([_PMs.var(pm, nw, :w)[bus][i] for i in 1:3]) <= vufmax/6 * sum([_PMs.var(pm, nw, :w)[bus][i] for i in 1:3]))
#    JuMP.@constraint(pm.model, wc - 1/3 * sum([_PMs.var(pm, nw, :w)[bus][i] for i in 1:3]) <= vufmax/6 * sum([_PMs.var(pm, nw, :w)[bus][i] for i in 1:3]))
#end

#function constraint_switch(pm::_PMs.AbstractPowerModel, xe, te)
#    JuMP.@constraint(pm.model, xe >= te)
#end

#function constraint_harden(pm::_PMs.AbstractPowerModel, xe, he)
#    JuMP.@constraint(pm.model, xe == he)
#end

#function constraint_mc_ohms_yt_from(pm::_PMs.AbstractPowerModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
#    p_fr  = _PMs.var(pm, n,  :p, f_idx)
#    q_fr  = _PMs.var(pm, n,  :q, f_idx)
#    vm_fr = _PMs.var(pm, n, :vm, f_bus)
#    vm_to = _PMs.var(pm, n, :vm, t_bus)
#    va_fr = _PMs.var(pm, n, :va, f_bus)
#    va_to = _PMs.var(pm, n, :va, t_bus)
#    z_branch = _PMs.var(pm, n, :z_branch, f_idx)

#    cnds = _PMD.conductor_ids(pm; nw=n)
#    for c in cnds
#        JuMP.@NLconstraint(pm.model, p_fr[c] ==
#                                        ((g[c,c]+g_fr[c,c])*vm_fr[c]^2
#                                        +sum( (g[c,d]+g_fr[c,d]) * vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d])
#                                             +(b[c,d]+b_fr[c,d]) * vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d])
#                                             for d in cnds if d != c)
#                                        +sum(-g[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d])
#                                             -b[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d])
#                                             for d in cnds)
#                                        ) * z_branch
#                                    )
#        JuMP.@NLconstraint(pm.model, q_fr[c] == (-(b[c,c]+b_fr[c,c])*vm_fr[c]^2
#                                        -sum( (b[c,d]+b_fr[c,d])*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d])
#                                             -(g[c,d]+g_fr[c,d])*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d])
#                                             for d in cnds if d != c)
#                                        -sum(-b[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d])
#                                             +g[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d])
#                                             for d in cnds)
#                                         * z_branch)
#                                    )
#    end
#end

#function constraint_mc_ohms_yt_to(pm::_PMs.AbstractPowerModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
#    constraint_mc_ohms_yt_from(pm, n, t_bus, f_bus, t_idx, f_idx, g, b, g_to, b_to, tr, ti, tm)
#end

#function constraint_activation_damage(pm::_PMs.AbstractPowerModel, nw::Int, arcs)
#    z = _PMs.var(pm, nw, :z_branch, (arcs))
#    xe = _PMs.var(pm, nw, :xe_d_s, (arcs))
#    te = _PMs.var(pm, nw, :te_d_s, (arcs))
#    JuMP.@constraint(pm.model, z == xe - te)
#end

#function constraint_activation_new(pm::_PMs.AbstractPowerModel, nw::Int, arcs)
#    z = _PMs.var(pm, nw, :z_branch, (arcs))
#    xe = _PMs.var(pm, nw, :xe_n_s, (arcs))
#    te = _PMs.var(pm, nw, :te_n_s, (arcs))
#    JuMP.@constraint(pm.model, z == xe - te)
#end

#function constraint_activation_active(pm::_PMs.AbstractPowerModel, nw::Int, arcs)
#    z = _PMs.var(pm, nw, :z_branch, (arcs))
#    JuMP.@constraint(pm.model, z == 1)
#end

#function constraint_variable_he_d(pm::_PMs.AbstractPowerModel, nw::Int, arcs)
#    he_s = _PMs.var(pm, nw, :he_d_s, (arcs))
#    he = _PMs.var(pm, pm.cnw, :he_d, (arcs))
#    JuMP.@constraint(pm.model, he >= he_s)
#end

#function constraint_variable_he_n(pm::_PMs.AbstractPowerModel, nw::Int, arcs)
#    he_s = _PMs.var(pm, nw, :he_n_s, (arcs))
#    he = _PMs.var(pm, pm.cnw, :he_n, (arcs))
#    JuMP.@constraint(pm.model, he >= he_s)
#end

#function constraint_variable_te_d(pm::_PMs.AbstractPowerModel, nw::Int, arcs)
#    te_s = _PMs.var(pm, nw, :te_d_s, (arcs))
#    te = _PMs.var(pm, pm.cnw, :te_d, (arcs))
#    JuMP.@constraint(pm.model, te >= te_s)
#end

#function constraint_variable_te_n(pm::_PMs.AbstractPowerModel, nw::Int, arcs)
#    te_s = _PMs.var(pm, nw, :te_n_s, (arcs))
#    te = _PMs.var(pm, pm.cnw, :te_n, (arcs))
#    JuMP.@constraint(pm.model, te >= te_s)
#end

#function constraint_variable_xe_d(pm::_PMs.AbstractPowerModel, nw::Int, arcs)
#    xe_s = _PMs.var(pm, nw, :xe_d_s, (arcs))
#    xe = _PMs.var(pm, pm.cnw, :xe_d, (arcs))
#    JuMP.@constraint(pm.model, xe >= xe_s)
#end

#function constraint_variable_xe_n(pm::_PMs.AbstractPowerModel, nw::Int, arcs)
#    xe_s = _PMs.var(pm, nw, :xe_n_s, (arcs))
#    xe = _PMs.var(pm, pm.cnw, :xe_n, (arcs))
#    JuMP.@constraint(pm.model, xe >= xe_s)
#end

#function constraint_mc_voltage_angle_difference(pm::_PMs.AbstractACPModel, nw::Int, f_idx, angmin, angmax)
#    i, f_bus, t_bus = f_idx

#    va_fr = _PMs.var(pm, nw, :va, f_bus)
#    va_to = _PMs.var(pm, nw, :va, t_bus)
#    z_branch = _PMs.var(pm, nw, :z_branch, (f_idx))
#
#    cnds = _PMD.conductor_ids(pm; nw=nw)

#    for c in cnds
#        JuMP.@NLconstraint(pm.model, z_branch*(va_fr[c] - va_to[c]) <= angmax[c])
#        JuMP.@NLconstraint(pm.model, z_branch*(va_fr[c] - va_to[c]) >= angmin[c])
#    end
#end

# add constraints for cycle elimination
#function constraint_cycle_elimination(pm, nw::Int, tours)
#    tour_sum = 0
#    for arcs in tours
#        ye = _PMs.var(pm, nw, :ye_s, arcs)
#        tour_sum += ye;
#    end
#    JuMP.@constraint(pm.model, tour_sum <= length(tours) - 1);
#end

#function constraint_cycle_function(pm, nw::Int, arcs)
#    z = _PMs.var(pm,  nw, :z_branch, arcs)
#    ye = _PMs.var(pm, nw, :ye_s, arcs)
#    JuMP.@constraint(pm.model, z <= ye)
#end
