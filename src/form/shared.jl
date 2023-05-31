
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
    psw_ne   = get(_PMD.var(pm, nw),   :psw_ne, Dict()); _PMD._check_var_keys(psw_ne, bus_arcs_sw_ne, "active power", "switch_inline_ne")
    qsw_ne   = get(_PMD.var(pm, nw),   :qsw_ne, Dict()); _PMD._check_var_keys(qsw_ne, bus_arcs_sw_ne, "reactive power", "switch_inline_ne")
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
