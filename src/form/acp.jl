""
function constraint_mc_power_balance_shed_ne(pm::_PMD.AbstractUnbalancedACPModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool},
                                             bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
                                             bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}},
                                             bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}},
                                             bus_shunts::Vector{Tuple{Int,Vector{Int}}}, bus_arcs_ne::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
                                             bus_arcs_sw_ne::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans_ne::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
                                             bus_gens_ne::Vector{Tuple{Int,Vector{Int}}})
    vm       = _PMD.var(pm, nw, :vm, i)
    va       = _PMD.var(pm, nw, :va, i)
    p        = get(_PMD.var(pm, nw),    :p, Dict());      _PMD._check_var_keys(p, bus_arcs, "active power", "branch")
    q        = get(_PMD.var(pm, nw),    :q, Dict());      _PMD._check_var_keys(q, bus_arcs, "reactive power", "branch")
    p_ne     = get(_PMD.var(pm, nw),    :p_ne, Dict());   _PMD._check_var_keys(p_ne, bus_arcs_ne, "active power", "branch_ne")
    q_ne     = get(_PMD.var(pm, nw),    :q_ne, Dict());   _PMD._check_var_keys(q_ne, bus_arcs_ne, "reactive power", "branch_ne")
    pg       = get(_PMD.var(pm, nw),    :pg, Dict());     _PMD._check_var_keys(pg, bus_gens, "active power", "generator")
    qg       = get(_PMD.var(pm, nw),    :qg, Dict());     _PMD._check_var_keys(qg, bus_gens, "reactive power", "generator")
    pg_ne    = get(_PMD.var(pm, nw),    :pg_ne, Dict());  _PMD._check_var_keys(pg, bus_gens_ne, "active power", "generator_ne")
    qg_ne    = get(_PMD.var(pm, nw),    :qg_ne, Dict());  _PMD._check_var_keys(qg, bus_gens_ne, "reactive power", "generator_ne")
    ps       = get(_PMD.var(pm, nw),    :ps, Dict());     _PMD._check_var_keys(ps, bus_storage, "active power", "storage")
    qs       = get(_PMD.var(pm, nw),    :qs, Dict());     _PMD._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw      = get(_PMD.var(pm, nw),    :psw, Dict());    _PMD._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw      = get(_PMD.var(pm, nw),    :qsw, Dict());    _PMD._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    psw_ne   = get(_PMD.var(pm, nw),    :psw_inline_ne, Dict()); _PMD._check_var_keys(psw, bus_arcs_sw_ne, "active power", "switch_ne")
    qsw_ne   = get(_PMD.var(pm, nw),    :qsw_inline_ne, Dict()); _PMD._check_var_keys(qsw, bus_arcs_sw_ne, "reactive power", "switch_ne")
    pt       = get(_PMD.var(pm, nw),    :pt, Dict());     _PMD._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt       = get(_PMD.var(pm, nw),    :qt, Dict());     _PMD._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")
    pt_ne    = get(_PMD.var(pm, nw),    :pt_ne, Dict());  _PMD._check_var_keys(pt, bus_arcs_trans_ne, "active power", "transformer_ne")
    qt_ne    = get(_PMD.var(pm, nw),    :qt_ne, Dict());  _PMD._check_var_keys(qt, bus_arcs_trans_ne, "reactive power", "transformer_ne")

    z_demand = _PMD.var(pm, nw, :z_demand)
    z_gen = haskey(_PMD.var(pm, nw), :z_gen) ? _PMD.var(pm, nw, :z_gen) : Dict(i => 1.0 for i in _PMD.ids(pm, nw, :gen))
    z_gen_ne = haskey(_PMD.var(pm, nw), :z_gen_ne) ? _PMD.var(pm, nw, :z_gen_ne) : Dict(i => 1.0 for i in _PMD.ids(pm, nw, :gen_ne))
    z_storage = haskey(_PMD.var(pm, nw), :z_storage) ? _PMD.var(pm, nw, :z_storage) : Dict(i => 1.0 for i in _PMD.ids(pm, nw, :storage))
    z_shunt  = haskey(_PMD.var(pm, nw), :z_shunt) ? _PMD.var(pm, nw, :z_shunt) : Dict(i => 1.0 for i in _PMD.ids(pm, nw, :shunt))

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx,t) in ungrounded_terminals
        cp = JuMP.@NLconstraint(pm.model,
              sum(p[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(psw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(pt[a_t][t] for (a_t, conns) in bus_arcs_trans if t in conns)
            - sum(pg[g][t]*z_gen[g] for (g, conns) in bus_gens if t in conns)
            + sum(ps[s][t]*z_storage[s] for (s, conns) in bus_storage if t in conns)
            + sum(_PMD.ref(pm, nw, :load, d, "pd")[findfirst(isequal(t), conns)]*z_demand[d] for (d, conns) in bus_loads if t in conns)
            + sum(z_shunt[s] *
                (_PMD.ref(pm, nw, :shunt, s)["gs"][findfirst(isequal(t), conns), findfirst(isequal(t), conns)] * vm[t]^2
                +sum( _PMD.ref(pm, nw, :shunt, s)["gs"][findfirst(isequal(t), conns), findfirst(isequal(u), conns)] * vm[t]*vm[u] * cos(va[t]-va[u])
                     + _PMD.ref(pm, nw, :shunt, s)["bs"][findfirst(isequal(t), conns), findfirst(isequal(u), conns)] * vm[t]*vm[u] * sin(va[t]-va[u])
                for (jdx, u) in ungrounded_terminals if idx != jdx ) )
            for (s, conns) in bus_shunts if t in conns )
            + sum(p_ne[a][t] for (a, conns) in bus_arcs_ne if t in conns)
            + sum(psw_ne[a_sw][t] for (a_sw, conns) in bus_arcs_sw_ne if t in conns)
            + sum(pt_ne[a_t][t] for (a_t, conns) in bus_arcs_trans_ne if t in conns)
            - sum(pg_ne[g][t]*z_gen_ne[g] for (g, conns) in bus_gens_ne if t in conns)
            ==
            0.0
        )
        push!(cstr_p, cp)

        cq = JuMP.@NLconstraint(pm.model,
              sum(q[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(qsw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(qt[a_t][t] for (a_t, conns) in bus_arcs_trans if t in conns)
            - sum(qg[g][t]*z_gen[g] for (g, conns) in bus_gens if t in conns)
            + sum(qs[s][t]*z_storage[s] for (s, conns) in bus_storage if t in conns)
            + sum(_PMD.ref(pm, nw, :load, l, "qd")[findfirst(isequal(t), conns)]*z_demand[l] for (l, conns) in bus_loads if t in conns)
            + sum(z_shunt[sh] *
                (-_PMD.ref(pm, nw, :shunt, sh)["bs"][findfirst(isequal(t), conns), findfirst(isequal(t), conns)] * vm[t]^2
                 -sum( _PMD.ref(pm, nw, :shunt, sh)["bs"][findfirst(isequal(t), conns), findfirst(isequal(u), conns)] * vm[t]*vm[u] * cos(va[t]-va[u])
                      -_PMD.ref(pm, nw, :shunt, sh)["gs"][findfirst(isequal(t), conns), findfirst(isequal(u), conns)] * vm[t]*vm[u] * sin(va[t]-va[u])
                for (jdx, u) in ungrounded_terminals if idx != jdx ) )
            for (sh, conns) in bus_shunts if t in conns )
            + sum(q_ne[a][t] for (a, conns) in bus_arcs_ne if t in conns)
            + sum(qsw_ne[a_sw][t] for (a_sw, conns) in bus_arcs_sw_ne if t in conns)
            + sum(qt_ne[a_t][t] for (a_t, conns) in bus_arcs_trans_ne if t in conns)
            - sum(qg_ne[g][t]*z_gen_ne[g] for (g, conns) in bus_gens_ne if t in conns)
            ==
            0.0
        )
        push!(cstr_q, cq)
    end

    _PMD.con(pm, nw, :lam_kcl_r)[i] = cstr_p
    _PMD.con(pm, nw, :lam_kcl_i)[i] = cstr_q

    if _IM.report_duals(pm)
        _PMD.sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        _PMD.sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end



@doc raw"""
    constraint_mc_ampacity_from_damaged(pm::AbstractUnbalancedACPModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on branches from-side for damaged branches

math```
p_{fr}^2 + q_{fr}^2 \leq vm_{fr}^2 i_{max}^2 *he_s
```
"""
function constraint_mc_ampacity_from_damaged(pm::_PMD.AbstractUnbalancedACPModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    p_fr  = [_PMD.var(pm, nw, :p, f_idx)[c] for c in f_connections]
    q_fr  = [_PMD.var(pm, nw, :q, f_idx)[c] for c in f_connections]
    vm_fr = [_PMD.var(pm, nw, :vm, f_idx[2])[c] for c in f_connections]
    he_s  = _PMD.var(pm, nw, :he_s, f_idx[1])

    # TODO: maybe introduce an auxillary varaible v_sqr = vm_fr[idx]^2, and do exact McCormick on v_sqr * he_s (and use @constraint)
    _PMD.con(pm, nw, :mu_cm_branch)[f_idx] = mu_cm_fr = [JuMP.@NLconstraint(pm.model, p_fr[idx]^2 + q_fr[idx]^2 <= vm_fr[idx]^2 * c_rating[idx]^2 * he_s) for idx in f_connections]

    if _IM.report_duals(pm)
        _PMD.sol(pm, nw, :branch, f_idx[1])[:mu_cm_fr] = mu_cm_fr
    end

    nothing
end


@doc raw"""
    constraint_mc_ampacity_to(pm::AbstractUnbalancedACPModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on branches to-side for damaged branches

math```
p_{to}^2 + q_{to}^2 \leq vm_{to}^2 i_{max}^2 * he_s
```
"""
function constraint_mc_ampacity_to_damaged(pm::_PMD.AbstractUnbalancedACPModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    p_to = [_PMD.var(pm, nw, :p, t_idx)[c] for c in t_connections]
    q_to = [_PMD.var(pm, nw, :q, t_idx)[c] for c in t_connections]
    vm_to = [_PMD.var(pm, nw, :vm, t_idx[2])[c] for c in t_connections]
    he_s  = _PMD.var(pm, nw, :he_s, t_idx[1])

    # TODO: maybe introduce an auxillary varaible v_sqr = vm_to[idx]^2, and do exact McCormick on v_sqr * he_s (and use @constraint)
    _PMD.con(pm, nw, :mu_cm_branch)[t_idx] = mu_cm_to = [JuMP.@NLconstraint(pm.model, p_to[idx]^2 + q_to[idx]^2 <= vm_to[idx]^2 * c_rating[idx]^2 * he_s) for idx in t_connections]

    if _IM.report_duals(pm)
        _PMD.sol(pm, nw, :branch, t_idx[1])[:mu_cm_to] = mu_cm_to
    end

    nothing
end



@doc raw"""
    constraint_mc_ampacity_from_ne(pm::AbstractUnbalancedACPModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on branches from-side for ne branches

math```
p_{fr}^2 + q_{fr}^2 \leq vm_{fr}^2 i_{max}^2 *xe_s
```
"""
function constraint_mc_ampacity_from_ne(pm::_PMD.AbstractUnbalancedACPModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    p_fr  = [_PMD.var(pm, nw, :p_ne, f_idx)[c] for c in f_connections]
    q_fr  = [_PMD.var(pm, nw, :q_ne, f_idx)[c] for c in f_connections]
    vm_fr = [_PMD.var(pm, nw, :vm, f_idx[2])[c] for c in f_connections]
    xe_s  = _PMD.var(pm, nw, :xe_s, f_idx[1])

    # TODO: maybe introduce an auxillary varaible v_sqr = vm_fr[idx]^2, and do exact McCormick on v_sqr * xe_s (and use @constraint)
    _PMD.con(pm, nw, :mu_cm_branch)[f_idx] = mu_cm_fr = [JuMP.@NLconstraint(pm.model, p_fr[idx]^2 + q_fr[idx]^2 <= vm_fr[idx]^2 * c_rating[idx]^2 * xe_s) for idx in f_connections]


    if _IM.report_duals(pm)
        _PMD.sol(pm, nw, :branch_ne, f_idx[1])[:mu_cm_fr_ne] = mu_cm_fr
    end

    nothing
end


@doc raw"""
    constraint_mc_ampacity_to_ne(pm::AbstractUnbalancedACPModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on branches to-side for ne branches

math```
p_{to}^2 + q_{to}^2 \leq vm_{to}^2 i_{max}^2 * xe_s
```
"""
function constraint_mc_ampacity_to_ne(pm::_PMD.AbstractUnbalancedACPModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    p_to = [_PMD.var(pm, nw, :p_ne, t_idx)[c] for c in t_connections]
    q_to = [_PMD.var(pm, nw, :q_ne, t_idx)[c] for c in t_connections]
    vm_to = [_PMD.var(pm, nw, :vm, t_idx[2])[c] for c in t_connections]
    xe_s  = _PMD.var(pm, nw, :xe_s, t_idx[1])

    # TODO: maybe introduce an auxillary varaible v_sqr = vm_to[idx]^2, and do exact McCormick on v_sqr * xe_s (and use @constraint)
    _PMD.con(pm, nw, :mu_cm_branch)[t_idx] = mu_cm_to = [JuMP.@NLconstraint(pm.model, p_to[idx]^2 + q_to[idx]^2 <= vm_to[idx]^2 * c_rating[idx]^2 * xe_s) for idx in t_connections]

    if _IM.report_duals(pm)
        _PMD.sol(pm, nw, :branch_ne, t_idx[1])[:mu_cm_to_ne] = mu_cm_to
    end

    nothing
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form) for damaged lines

```
p_fr ==    he * g[c,c] * vm_fr[c]^2 +
            sum( g[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) +
                 b[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in conductor_ids(pm) if d != c) +
            sum(-g[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                -b[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in conductor_ids(pm))
            + g_fr[c,c] * vm_fr[c]^2 +
            sum( g_fr[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) +
                 b_fr[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in conductor_ids(pm) if d != c)
            )
q_fr == he * -b[c,c] *vm_fr[c]^2 -
            sum( b[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) -
                 g[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in conductor_ids(pm) if d != c) -
            sum(-b[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                 g[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in conductor_ids(pm))
            -b_fr[c,c] *vm_fr[c]^2 -
            sum( b_fr[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) -
                 g_fr[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in conductor_ids(pm) if d != c)
            )
```
"""
function constraint_mc_ohms_yt_from_damaged(pm::_PMD.AbstractUnbalancedACPModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real}, vad_min::Vector{<:Real}, vad_max::Vector{<:Real})
    p_fr  = _PMD.var(pm, nw,  :p, f_idx)
    q_fr  = _PMD.var(pm, nw,  :q, f_idx)
    vm_fr = _PMD.var(pm, nw, :vm, f_bus)
    vm_to = _PMD.var(pm, nw, :vm, t_bus)
    va_fr = _PMD.var(pm, nw, :va, f_bus)
    va_to = _PMD.var(pm, nw, :va, t_bus)
    he_s  = _PMD.var(pm, nw, :he_s, f_idx[1])

    ohms_yt_p = JuMP.ConstraintRef[]
    ohms_yt_q = JuMP.ConstraintRef[]
    for (idx, (fc,tc)) in enumerate(zip(f_connections,t_connections))
        push!(ohms_yt_p, JuMP.@NLconstraint(pm.model, p_fr[fc] == he_s * ((G[idx,idx]+G_fr[idx,idx])*vm_fr[fc]^2
            +sum( (G[idx,jdx]+G_fr[idx,jdx]) * vm_fr[fc]*vm_fr[fd]*cos(va_fr[fc]-va_fr[fd])
                 +(B[idx,jdx]+B_fr[idx,jdx]) * vm_fr[fc]*vm_fr[fd]*sin(va_fr[fc]-va_fr[fd])
                for (jdx, (fd,td)) in enumerate(zip(f_connections,t_connections)) if idx != jdx)
            +sum( -G[idx,jdx]*vm_fr[fc]*vm_to[td]*cos(va_fr[fc]-va_to[td])
                  -B[idx,jdx]*vm_fr[fc]*vm_to[td]*sin(va_fr[fc]-va_to[td])
                for (jdx, (fd,td)) in enumerate(zip(f_connections,t_connections))))
            )
        )

        push!(ohms_yt_q, JuMP.@NLconstraint(pm.model, q_fr[fc] == he_s * (-(B[idx,idx]+B_fr[idx,idx])*vm_fr[fc]^2
            -sum( (B[idx,jdx]+B_fr[idx,jdx])*vm_fr[fc]*vm_fr[fd]*cos(va_fr[fc]-va_fr[fd])
                 -(G[idx,jdx]+G_fr[idx,jdx])*vm_fr[fc]*vm_fr[fd]*sin(va_fr[fc]-va_fr[fd])
                for (jdx, (fd,td)) in enumerate(zip(f_connections,t_connections)) if idx != jdx)
            -sum(-B[idx,jdx]*vm_fr[fc]*vm_to[td]*cos(va_fr[fc]-va_to[td])
                 +G[idx,jdx]*vm_fr[fc]*vm_to[td]*sin(va_fr[fc]-va_to[td])
                for (jdx, (fd,td)) in enumerate(zip(f_connections,t_connections))))
            )
        )
    end
    _PMD.con(pm, nw, :ohms_yt)[f_idx] = [ohms_yt_p, ohms_yt_q]
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form) for damaged lines

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_mc_ohms_yt_to_damaged(pm::_PMD.AbstractUnbalancedACPModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_to::Matrix{<:Real}, B_to::Matrix{<:Real}, vad_min::Vector{<:Real}, vad_max::Vector{<:Real})
    constraint_mc_ohms_yt_from_damaged(pm, nw, t_bus, f_bus, t_idx, f_idx, t_connections, f_connections, G, B, G_to, B_to, vad_min, vad_max)
end



"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form) for ne lines

```
p_fr ==    he * g[c,c] * vm_fr[c]^2 +
            sum( g[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) +
                 b[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in conductor_ids(pm) if d != c) +
            sum(-g[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                -b[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in conductor_ids(pm))
            + g_fr[c,c] * vm_fr[c]^2 +
            sum( g_fr[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) +
                 b_fr[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in conductor_ids(pm) if d != c)
            )
q_fr == he * -b[c,c] *vm_fr[c]^2 -
            sum( b[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) -
                 g[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in conductor_ids(pm) if d != c) -
            sum(-b[c,d]*vm_fr[c]*vm_to[d]*cos(va_fr[c]-va_to[d]) +
                 g[c,d]*vm_fr[c]*vm_to[d]*sin(va_fr[c]-va_to[d]) for d in conductor_ids(pm))
            -b_fr[c,c] *vm_fr[c]^2 -
            sum( b_fr[c,d]*vm_fr[c]*vm_fr[d]*cos(va_fr[c]-va_fr[d]) -
                 g_fr[c,d]*vm_fr[c]*vm_fr[d]*sin(va_fr[c]-va_fr[d]) for d in conductor_ids(pm) if d != c)
            )
```
"""
function constraint_mc_ohms_yt_from_ne(pm::_PMD.AbstractUnbalancedACPModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real}, vad_min::Vector{<:Real}, vad_max::Vector{<:Real})
    p_fr  = _PMD.var(pm, nw,  :p_ne, f_idx)
    q_fr  = _PMD.var(pm, nw,  :q_ne, f_idx)
    vm_fr = _PMD.var(pm, nw, :vm, f_bus)
    vm_to = _PMD.var(pm, nw, :vm, t_bus)
    va_fr = _PMD.var(pm, nw, :va, f_bus)
    va_to = _PMD.var(pm, nw, :va, t_bus)
    xe_s  = _PMD.var(pm, nw, :xe_s, f_idx[1])

    ohms_yt_p = JuMP.ConstraintRef[]
    ohms_yt_q = JuMP.ConstraintRef[]
    for (idx, (fc,tc)) in enumerate(zip(f_connections,t_connections))
        push!(ohms_yt_p, JuMP.@NLconstraint(pm.model, p_fr[fc] == xe_s * ((G[idx,idx]+G_fr[idx,idx])*vm_fr[fc]^2
            +sum( (G[idx,jdx]+G_fr[idx,jdx]) * vm_fr[fc]*vm_fr[fd]*cos(va_fr[fc]-va_fr[fd])
                 +(B[idx,jdx]+B_fr[idx,jdx]) * vm_fr[fc]*vm_fr[fd]*sin(va_fr[fc]-va_fr[fd])
                for (jdx, (fd,td)) in enumerate(zip(f_connections,t_connections)) if idx != jdx)
            +sum( -G[idx,jdx]*vm_fr[fc]*vm_to[td]*cos(va_fr[fc]-va_to[td])
                  -B[idx,jdx]*vm_fr[fc]*vm_to[td]*sin(va_fr[fc]-va_to[td])
                for (jdx, (fd,td)) in enumerate(zip(f_connections,t_connections))))
            )
        )

        push!(ohms_yt_q, JuMP.@NLconstraint(pm.model, q_fr[fc] == xe_s * (-(B[idx,idx]+B_fr[idx,idx])*vm_fr[fc]^2
            -sum( (B[idx,jdx]+B_fr[idx,jdx])*vm_fr[fc]*vm_fr[fd]*cos(va_fr[fc]-va_fr[fd])
                 -(G[idx,jdx]+G_fr[idx,jdx])*vm_fr[fc]*vm_fr[fd]*sin(va_fr[fc]-va_fr[fd])
                for (jdx, (fd,td)) in enumerate(zip(f_connections,t_connections)) if idx != jdx)
            -sum(-B[idx,jdx]*vm_fr[fc]*vm_to[td]*cos(va_fr[fc]-va_to[td])
                 +G[idx,jdx]*vm_fr[fc]*vm_to[td]*sin(va_fr[fc]-va_to[td])
                for (jdx, (fd,td)) in enumerate(zip(f_connections,t_connections))))
            )
        )
    end
    _PMD.con(pm, nw, :ohms_yt)[f_idx] = [ohms_yt_p, ohms_yt_q]
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form) for ne lines

```
p[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_mc_ohms_yt_to_ne(pm::_PMD.AbstractUnbalancedACPModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_to::Matrix{<:Real}, B_to::Matrix{<:Real}, vad_min::Vector{<:Real}, vad_max::Vector{<:Real})
    constraint_mc_ohms_yt_from_ne(pm, nw, t_bus, f_bus, t_idx, f_idx, t_connections, f_connections, G, B, G_to, B_to, vad_min, vad_max)
end



@doc raw"""
    constraint_mc_switch_state_voltage_open_closed(pm::PMD.AbstractUnbalancedACPModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_connections::Vector{Int}, t_connections::Vector{Int})

Linear switch power on/off constraint for ACPU form.

```math
\begin{align}
& |V^{fr}_{i,c}| - |V^{to}_{i,c}| \leq \left ( v^u_{i,c} - v^l_{i,c} \right ) \left ( 1 - z^{sw}_i \right )\ \forall i \in S,\forall c \in C \\
& |V^{fr}_{i,c}| - |V^{to}_{i,c}| \geq -\left ( v^u_{i,c} - v^l_{i,c} \right ) \left ( 1 - z^{sw}_i \right )\ \forall i \in S,\forall c \in C \\

\end{align}
```
"""
function constraint_mc_switch_inline_ne_voltage_open_close(pm::_PMD.AbstractUnbalancedACPModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_connections::Vector{Int}, t_connections::Vector{Int})
    vm_fr = _PMD.var(pm, nw, :vm, f_bus)
    vm_to = _PMD.var(pm, nw, :vm, t_bus)
    va_fr = _PMD.var(pm, nw, :va, f_bus)
    va_to = _PMD.var(pm, nw, :va, t_bus)

    f_bus = _PMD.ref(pm, nw, :bus, f_bus)
    t_bus = _PMD.ref(pm, nw, :bus, t_bus)

    f_vmin = f_bus["vmin"][[findfirst(isequal(c), f_bus["terminals"]) for c in f_connections]]
    t_vmin = t_bus["vmin"][[findfirst(isequal(c), t_bus["terminals"]) for c in t_connections]]

    f_vmax = f_bus["vmax"][[findfirst(isequal(c), f_bus["terminals"]) for c in f_connections]]
    t_vmax = t_bus["vmax"][[findfirst(isequal(c), t_bus["terminals"]) for c in t_connections]]

    vmin = max.(fill(0.0, length(f_vmax)), f_vmin, t_vmin)
    vmax = min.(fill(2.0, length(f_vmax)), f_vmax, t_vmax)

    angmin = get(_PMD.ref(pm, nw, :switch_inline_ne, i), "angmin", deg2rad.(fill(-5.0, length(f_connections))))
    angmax = get(_PMD.ref(pm, nw, :switch_inline_ne, i), "angmax", deg2rad.(fill( 5.0, length(f_connections))))

    state = _PMD.var(pm, nw, :switch_inline_ne_state, i)

    for (idx, (fc, tc)) in enumerate(zip(f_connections, t_connections))
        JuMP.@constraint(pm.model, vm_fr[fc] - vm_to[tc] <=  (vmax[idx]-vmin[idx]) * (1-state))
        JuMP.@constraint(pm.model, vm_fr[fc] - vm_to[tc] >= -(vmax[idx]-vmin[idx]) * (1-state))

        JuMP.@constraint(pm.model, va_fr[fc] - va_to[tc] <=  (angmax[idx]-angmin[idx]) * (1-state))
        JuMP.@constraint(pm.model, va_fr[fc] - va_to[tc] >= -(angmax[idx]-angmin[idx]) * (1-state))

        # Indicator constraint version, for reference
        # JuMP.@constraint(pm.model, state => {vm_fr[fc] == vm_to[tc]})
        # JuMP.@constraint(pm.model, state => {va_fr[fc] == va_to[tc]})
    end
end

@doc raw"""
    constraint_mc_switch_ampacity(pm::AbstractUnbalancedACPModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on switches

math```
p_{fr}^2 + q_{fr}^2 \leq vm_{fr}^2 i_{max}^2
```
"""
function constraint_mc_switch_inline_ne_ampacity(pm::_PMD.AbstractUnbalancedACPModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    psw_fr = [_PMD.var(pm, nw, :psw_inline_ne, f_idx)[c] for c in f_connections]
    qsw_fr = [_PMD.var(pm, nw, :qsw_inline_ne, f_idx)[c] for c in f_connections]
    vm_fr = [_PMD.var(pm, nw, :vm, f_idx[2])[c] for c in f_connections]

    _PMD.con(pm, nw, :mu_cm_switch)[f_idx] = mu_cm_fr = [JuMP.@constraint(pm.model, psw_fr[idx]^2 + qsw_fr[idx]^2 .<= vm_fr[idx]^2 * c_rating[idx]^2) for idx in findall(c_rating .< Inf)]

    if _IM.report_duals(pm)
        _PMD.sol(pm, nw, :switch_inline_ne, f_idx[1])[:mu_cm_fr] = mu_cm_fr
    end

    nothing
end
