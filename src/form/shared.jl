
function variable_he_s(pm::_PMD.AbstractUnbalancedWModels; nw::Int=_PMD.nw_id_default, relax::Bool=false, report::Bool=true)
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

     terminals = Dict(i => bus["terminals"] for (i,bus) in _PMD.ref(pm, nw, :bus))
     branch = _PMD.ref(pm, nw, :branch)

     he_s_w_fr = _PMD.var(pm, nw)[:he_s_w_fr] = Dict(i => JuMP.@variable(pm.model,
                 [t in terminals[branch[i]["f_bus"]]], base_name="$(nw)_he_s_w_fr_$(i)",
                 start=_PMD.comp_start_value(_PMD.ref(pm, nw, :branch, i), "he_start", i, 0.0)
             ) for i in _PMD.ref(pm, nw, :branch_harden)
         )

    he_s_w_to = _PMD.var(pm, nw)[:he_s_w_to] = Dict(i => JuMP.@variable(pm.model,
                [t in terminals[branch[i]["t_bus"]]], base_name="$(nw)_he_s_w_to_$(i)",
                start=_PMD.comp_start_value(_PMD.ref(pm, nw, :branch, i), "he_start", i, 0.0)
            ) for i in _PMD.ref(pm, nw, :branch_harden)
         )

    report && _IM.sol_component_value(pm, _PMD.pmd_it_sym, nw, :branch, :he_s, _PMD.ref(pm, nw, :branch_harden), he_s)
end


function variable_xe_s(pm::_PMD.AbstractUnbalancedWModels; nw::Int=_PMD.nw_id_default, relax::Bool=false, report::Bool=true)
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

     terminals = Dict(i => bus["terminals"] for (i,bus) in _PMD.ref(pm, nw, :bus))
     branch = _PMD.ref(pm, nw, :branch_ne)

     xe_s_w_fr = _PMD.var(pm, nw)[:xe_s_w_fr] = Dict(i => JuMP.@variable(pm.model,
                 [t in terminals[branch[i]["f_bus"]]], base_name="$(nw)_xe_s_w_fr_$(i)",
                 start=_PMD.comp_start_value(_PMD.ref(pm, nw, :branch_ne, i), "xe_start", i, 0.0)
             ) for i in _PMD.ids(pm, nw, :branch_ne)
         )

     xe_s_w_to = _PMD.var(pm, nw)[:xe_s_w_to] = Dict(i => JuMP.@variable(pm.model,
                [t in terminals[branch[i]["t_bus"]]], base_name="$(nw)_xe_s_w_to_$(i)",
                start=_PMD.comp_start_value(_PMD.ref(pm, nw, :branch_ne, i), "xe_start", i, 0.0)
            ) for i in _PMD.ids(pm, nw, :branch_ne)
         )


    report && _IM.sol_component_value(pm, _PMD.pmd_it_sym, nw, :branch_ne, :xe_s, _PMD.ids(pm, nw, :branch_ne), xe_s)
end


"KCL for load shed problem with transformers (AbstractWForms)"
function constraint_mc_power_balance_shed_ne(pm::_PMD.AbstractUnbalancedWModels, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool},
                                          bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
                                          bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}},
                                          bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}},
                                          bus_arcs_ne::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
                                          bus_arcs_sw_ne::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans_ne::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
                                          bus_gens_ne::Vector{Tuple{Int,Vector{Int}}})

    w        = _PMD.var(pm, nw, :w, i)
    p        = get(_PMD.var(pm, nw),   :p, Dict()); _PMD._check_var_keys(p, bus_arcs, "active power", "branch")
    q        = get(_PMD.var(pm, nw),   :q, Dict()); _PMD._check_var_keys(q, bus_arcs, "reactive power", "branch")
    p_ne     = get(_PMD.var(pm, nw),   :p_ne, Dict()); _PMD._check_var_keys(p_ne, bus_arcs_ne, "active power", "branch_ne")
    q_ne     = get(_PMD.var(pm, nw),   :q_ne, Dict()); _PMD._check_var_keys(q_ne, bus_arcs_ne, "reactive power", "branch_ne")
    pg       = get(_PMD.var(pm, nw),   :pg, Dict()); _PMD._check_var_keys(pg, bus_gens, "active power", "generator")
    qg       = get(_PMD.var(pm, nw),   :qg, Dict()); _PMD._check_var_keys(qg, bus_gens, "reactive power", "generator")
    pg_ne    = get(_PMD.var(pm, nw),   :pg_ne, Dict()); _PMD._check_var_keys(pg_ne, bus_gens_ne, "active power", "generator_ne")
    qg_ne    = get(_PMD.var(pm, nw),   :qg_ne, Dict()); _PMD._check_var_keys(qg_ne, bus_gens_ne, "reactive power", "generator_ne")
    ps       = get(_PMD.var(pm, nw),   :ps, Dict()); _PMD._check_var_keys(ps, bus_storage, "active power", "storage")
    qs       = get(_PMD.var(pm, nw),   :qs, Dict()); _PMD._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw      = get(_PMD.var(pm, nw),   :psw, Dict()); _PMD._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw      = get(_PMD.var(pm, nw),   :qsw, Dict()); _PMD._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    psw_ne   = get(_PMD.var(pm, nw),   :psw_inline_ne, Dict()); _PMD._check_var_keys(psw_ne, bus_arcs_sw_ne, "active power", "switch_inline_ne")
    qsw_ne   = get(_PMD.var(pm, nw),   :qsw_inline_ne, Dict()); _PMD._check_var_keys(qsw_ne, bus_arcs_sw_ne, "reactive power", "switch_inline_ne")
    pt       = get(_PMD.var(pm, nw),   :pt, Dict()); _PMD._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt       = get(_PMD.var(pm, nw),   :qt, Dict()); _PMD._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")
    pt_ne    = get(_PMD.var(pm, nw),   :pt_ne, Dict()); _PMD._check_var_keys(pt_ne, bus_arcs_trans_ne, "active power", "transformer_ne")
    qt_ne    = get(_PMD.var(pm, nw),   :qt_ne, Dict()); _PMD._check_var_keys(qt_ne, bus_arcs_trans_ne, "reactive power", "transformer_ne")
    z_demand = _PMD.var(pm, nw, :z_demand)
    z_shunt  = _PMD.var(pm, nw, :z_shunt)

    Gt, Bt = _PMD._build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx, t) in ungrounded_terminals
        cp = JuMP.@constraint(pm.model,
              sum(p[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(psw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(pt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
            + sum(p_ne[a][t] for (a, conns) in bus_arcs_ne if t in conns)
            + sum(psw_ne[a_sw][t] for (a_sw, conns) in bus_arcs_sw_ne if t in conns)
            + sum(pt_ne[a_trans][t] for (a_trans, conns) in bus_arcs_trans_ne if t in conns)
            ==
            sum(pg[g][t] for (g, conns) in bus_gens if t in conns)
            + sum(pg_ne[g][t] for (g, conns) in bus_gens_ne if t in conns)
            - sum(ps[s][t] for (s, conns) in bus_storage if t in conns)
            - sum(_PMD.ref(pm, nw, :load, l, "pd")[findfirst(isequal(t), conns)] * z_demand[l] for (l, conns) in bus_loads if t in conns)
            - sum(z_shunt[sh] *(w[t] * LinearAlgebra.diag(Gt')[idx]) for (sh, conns) in bus_shunts if t in conns)
        )
        push!(cstr_p, cp)
        cq = JuMP.@constraint(pm.model,
              sum(q[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(qsw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(qt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
            + sum(q_ne[a][t] for (a, conns) in bus_arcs_ne if t in conns)
            + sum(qsw_ne[a_sw][t] for (a_sw, conns) in bus_arcs_sw_ne if t in conns)
            + sum(qt_ne[a_trans][t] for (a_trans, conns) in bus_arcs_trans_ne if t in conns)
            ==
            sum(qg[g][t] for (g, conns) in bus_gens if t in conns)
            + sum(qg_ne[g][t] for (g, conns) in bus_gens_ne if t in conns)
            - sum(qs[s][t] for (s, conns) in bus_storage if t in conns)
            - sum(_PMD.ref(pm, nw, :load, l, "qd")[findfirst(isequal(t), conns)]*z_demand[l] for (l, conns) in bus_loads if t in conns)
            - sum(z_shunt[sh] * (-w[t] * LinearAlgebra.diag(Bt')[idx]) for (sh, conns) in bus_shunts if t in conns)
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
    constraint_mc_ampacity_from_damaged(pm::AbstractUnbalancedWModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on branches from-side that are damaged

math```
p_{fr}^2 + q_{fr}^2 \leq w_{fr} i_{max}^2 * he_s
```
"""
function constraint_mc_ampacity_from_damaged(pm::_PMD.AbstractUnbalancedWModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    p_fr       = [_PMD.var(pm, nw, :p, f_idx)[c] for c in f_connections]
    q_fr       = [_PMD.var(pm, nw, :q, f_idx)[c] for c in f_connections]
    w_fr       = [_PMD.var(pm, nw, :w, f_idx[2])[c] for c in f_connections]
    he_s       = _PMD.var(pm, nw, :he_s, f_idx[1])
    he_s_w_fr  = [_PMD.var(pm, nw, :he_s_w_fr, f_idx[1])[c] for c in f_connections]

    for c in f_connections
       _IM.relaxation_product(pm.model, he_s, w_fr[c], he_s_w_fr[c])
    end

#    _PMD.con(pm, nw, :mu_cm_branch)[f_idx] = mu_cm_fr = [JuMP.@constraint(pm.model, p_fr[idx]^2 + q_fr[idx]^2 .<= w_fr[idx] * c_rating[idx]^2 * he_s) for idx in findall(c_rating .< Inf)]
    _PMD.con(pm, nw, :mu_cm_branch)[f_idx] = mu_cm_fr = [JuMP.@constraint(pm.model, p_fr[idx]^2 + q_fr[idx]^2 .<= he_s_w_fr[idx] * c_rating[idx]^2) for idx in f_connections]


    if _IM.report_duals(pm)
        _PMD.sol(pm, nw, :branch, f_idx[1])[:mu_cm_fr] = mu_cm_fr
    end

    nothing
end


@doc raw"""
    constraint_mc_ampacity_to(pm::AbstractUnbalancedWModels, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on branches to-side that are damaged

math```
p_{to}^2 + q_{to}^2 \leq w_{to} i_{max}^2
```
"""
function constraint_mc_ampacity_to_damaged(pm::_PMD.AbstractUnbalancedWModels, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    p_to       = [_PMD.var(pm, nw, :p, t_idx)[c] for c in t_connections]
    q_to       = [_PMD.var(pm, nw, :q, t_idx)[c] for c in t_connections]
    w_to       = [_PMD.var(pm, nw, :w, t_idx[2])[c] for c in t_connections]
    he_s       = _PMD.var(pm, nw, :he_s, t_idx[1])
    he_s_w_to  = [_PMD.var(pm, nw, :he_s_w_to, t_idx[1])[c] for c in t_connections]

    for c in t_connections
       _IM.relaxation_product(pm.model, he_s, w_to[c], he_s_w_to[c])
    end

#    _PMD.con(pm, nw, :mu_cm_branch)[t_idx] = mu_cm_to = [JuMP.@constraint(pm.model, p_to[idx]^2 + q_to[idx]^2 .<= w_to[idx] * c_rating[idx]^2 * he_s) for idx in findall(c_rating .< Inf)]
    _PMD.con(pm, nw, :mu_cm_branch)[t_idx] = mu_cm_to = [JuMP.@constraint(pm.model, p_to[idx]^2 + q_to[idx]^2 .<= he_s_w_to[idx] * c_rating[idx]^2) for idx in t_connections]

    if _IM.report_duals(pm)
        _PMD.sol(pm, nw, :branch, t_idx[1])[:mu_cm_to] = mu_cm_to
    end

    nothing
end



@doc raw"""
    constraint_mc_ampacity_from_ne(pm::AbstractUnbalancedWModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on branches from-side that are ne

math```
p_{fr}^2 + q_{fr}^2 \leq w_{fr} i_{max}^2 * xe_s
```
"""
function constraint_mc_ampacity_from_ne(pm::_PMD.AbstractUnbalancedWModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    p_fr       = [_PMD.var(pm, nw, :p_ne, f_idx)[c] for c in f_connections]
    q_fr       = [_PMD.var(pm, nw, :q_ne, f_idx)[c] for c in f_connections]
    w_fr       = [_PMD.var(pm, nw, :w, f_idx[2])[c] for c in f_connections]
    xe_s       = _PMD.var(pm, nw, :xe_s, f_idx[1])
    xe_s_w_fr  = [_PMD.var(pm, nw, :xe_s_w_fr, f_idx[1])[c] for c in f_connections]

    for c in f_connections
       _IM.relaxation_product(pm.model, xe_s, w_fr[c], xe_s_w_fr[c])
    end

#    _PMD.con(pm, nw, :mu_cm_branch)[f_idx] = mu_cm_fr = [JuMP.@constraint(pm.model, p_fr[idx]^2 + q_fr[idx]^2 .<= w_fr[idx] * c_rating[idx]^2 * he_s) for idx in findall(c_rating .< Inf)]
    _PMD.con(pm, nw, :mu_cm_branch_ne)[f_idx] = mu_cm_fr = [JuMP.@constraint(pm.model, p_fr[idx]^2 + q_fr[idx]^2 .<= xe_s_w_fr[idx] * c_rating[idx]^2) for idx in f_connections]


    if _IM.report_duals(pm)
        _PMD.sol(pm, nw, :branch_ne, f_idx[1])[:mu_cm_fr_ne] = mu_cm_fr
    end

    nothing
end


@doc raw"""
    constraint_mc_ampacity_to_ne(pm::AbstractUnbalancedWModels, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on branches to-side that are ned
math```
p_{to}^2 + q_{to}^2 \leq w_{to} i_{max}^2 # xe_s
```
"""
function constraint_mc_ampacity_to_ne(pm::_PMD.AbstractUnbalancedWModels, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    p_to       = [_PMD.var(pm, nw, :p_ne, t_idx)[c] for c in t_connections]
    q_to       = [_PMD.var(pm, nw, :q_ne, t_idx)[c] for c in t_connections]
    w_to       = [_PMD.var(pm, nw, :w, t_idx[2])[c] for c in t_connections]
    xe_s       = _PMD.var(pm, nw, :xe_s, t_idx[1])
    xe_s_w_to  = [_PMD.var(pm, nw, :xe_s_w_to, t_idx[1])[c] for c in t_connections]

    for c in t_connections
       _IM.relaxation_product(pm.model, xe_s, w_to[c], xe_s_w_to[c])
    end

#    _PMD.con(pm, nw, :mu_cm_branch)[t_idx] = mu_cm_to = [JuMP.@constraint(pm.model, p_to[idx]^2 + q_to[idx]^2 .<= w_to[idx] * c_rating[idx]^2 * he_s) for idx in findall(c_rating .< Inf)]
    _PMD.con(pm, nw, :mu_cm_branch_ne)[t_idx] = mu_cm_to = [JuMP.@constraint(pm.model, p_to[idx]^2 + q_to[idx]^2 .<= xe_s_w_to[idx] * c_rating[idx]^2) for idx in t_connections]

    if _IM.report_duals(pm)
        _PMD.sol(pm, nw, :branch_ne, t_idx[1])[:mu_cm_to_ne] = mu_cm_to
    end

    nothing
end

""
function constraint_mc_voltage_angle_difference_damaged(pm::_PMD.AbstractUnbalancedPolarModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, angmin::Vector{<:Real}, angmax::Vector{<:Real})
    i, f_bus, t_bus = f_idx

    va_fr = [_PMD.var(pm, nw, :va, f_bus)[fc] for fc in f_connections]
    va_to = [_PMD.var(pm, nw, :va, t_bus)[tc] for tc in t_connections]
    he_s  = _PMD.var(pm, nw, :he_s, i)

    #TODO: A bit lazy, but this is how PowerModels.jl does on_off on phase angle difference constraints because va is unbounded.  In practice, if he_s is binary, we should do McCormick
    # on the lhs, and this constraint otherwise.  If we have an analytical result on the largest phase angle difference across the network, we could also make this a big-M style on-off constraint.
    JuMP.@constraint(pm.model, he_s * (va_fr .- va_to) .<= he_s * angmax)
    JuMP.@constraint(pm.model, he_s * (va_fr .- va_to) .>= he_s * angmin)
end


""
function constraint_mc_voltage_angle_difference_ne(pm::_PMD.AbstractUnbalancedPolarModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, angmin::Vector{<:Real}, angmax::Vector{<:Real})
    i, f_bus, t_bus = f_idx

    va_fr = [_PMD.var(pm, nw, :va, f_bus)[fc] for fc in f_connections]
    va_to = [_PMD.var(pm, nw, :va, t_bus)[tc] for tc in t_connections]
    xe_s  = _PMD.var(pm, nw, :xe_s, i)

    #TODO: A bit lazy, but this is how PowerModels.jl does on_off on phase angle difference constraints because va is unbounded.  In practice, if he_s is binary, we should do McCormick
    # on the lhs, and this constraint otherwise.  If we have an analytical result on the largest phase angle difference across the network, we could also make this a big-M style on-off constraint.
    JuMP.@constraint(pm.model, xe_s * (va_fr .- va_to) .<= xe_s * angmax)
    JuMP.@constraint(pm.model, xe_s * (va_fr .- va_to) .>= xe_s * angmin)
end


@doc raw"""
    constraint_mc_switch_inline_ne_ampacity(pm::AbstractUnbalancedWModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on switches from-side

math```
p_{fr}^2 + q_{fr}^2 \leq w_{fr} i_{max}^2
```
"""
function constraint_mc_switch_inline_ne_ampacity(pm::_PMD.AbstractUnbalancedWModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    psw_fr = [_PMD.var(pm, nw, :psw_inline_ne, f_idx)[c] for c in f_connections]
    qsw_fr = [_PMD.var(pm, nw, :qsw_inline_ne, f_idx)[c] for c in f_connections]
    w_fr = [_PMD.var(pm, nw, :w, f_idx[2])[c] for c in f_connections]

    _PMD.con(pm, nw, :mu_cm_switch_inline_ne)[f_idx] = mu_cm_fr = [JuMP.@constraint(pm.model, psw_fr[idx]^2 + qsw_fr[idx]^2 .<= w_fr[idx] * c_rating[idx]^2) for idx in findall(c_rating .< Inf)]

    if _IM.report_duals(pm)
        _PMD.sol(pm, nw, :switch_inline_ne, f_idx[1])[:mu_cm_fr] = mu_cm_fr
    end

    nothing
end
