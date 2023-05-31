
function constraint_mc_power_ne_balance_shed(pm::_PMD.AbstractUnbalancedACRModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool},
                                             bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
                                             bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}},
                                             bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}},
                                             bus_shunts::Vector{Tuple{Int,Vector{Int}}}, bus_arcs_ne::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
                                             bus_arcs_sw_ne::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans_ne::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
                                             bus_gens_ne::Vector{Tuple{Int,Vector{Int}}})
    vr     = _PMD.var(pm, nw, :vr, i)
    vi     = _PMD.var(pm, nw, :vi, i)
    p      = get(_PMD.var(pm, nw), :p,      Dict()); _PMD._check_var_keys(p,   bus_arcs,          "active power",   "branch")
    q      = get(_PMD.var(pm, nw), :q,      Dict()); _PMD._check_var_keys(q,   bus_arcs,          "reactive power", "branch")
    p_ne   = get(_PMD.var(pm, nw), :p_ne,   Dict()); _PMD._check_var_keys(p,   bus_arcs_ne,       "active power",   "branch_ne")
    q_ne   = get(_PMD.var(pm, nw), :q_ne,   Dict()); _PMD._check_var_keys(q,   bus_arcs_ne,       "reactive power", "branch_ne")
    pg     = get(_PMD.var(pm, nw), :pg,     Dict()); _PMD._check_var_keys(pg,  bus_gens,          "active power",   "generator")
    qg     = get(_PMD.var(pm, nw), :qg,     Dict()); _PMD._check_var_keys(qg,  bus_gens,          "reactive power", "generator")
    pg_ne  = get(_PMD.var(pm, nw), :pg_ne,  Dict()); _PMD._check_var_keys(pg,  bus_gens_ne,       "active power",   "generator_ne")
    qg_ne  = get(_PMD.var(pm, nw), :qg_ne,  Dict()); _PMD._check_var_keys(qg,  bus_gens_ne,       "reactive power", "generator_ne")
    ps     = get(_PMD.var(pm, nw), :ps,     Dict()); _PMD._check_var_keys(ps,  bus_storage,       "active power",   "storage")
    qs     = get(_PMD.var(pm, nw), :qs,     Dict()); _PMD._check_var_keys(qs,  bus_storage,       "reactive power", "storage")
    psw    = get(_PMD.var(pm, nw), :psw,    Dict()); _PMD._check_var_keys(psw, bus_arcs_sw,       "active power",   "switch")
    qsw    = get(_PMD.var(pm, nw), :qsw,    Dict()); _PMD._check_var_keys(qsw, bus_arcs_sw,       "reactive power", "switch")
    psw_ne = get(_PMD.var(pm, nw), :psw_ne, Dict()); _PMD._check_var_keys(psw, bus_arcs_sw_ne,    "active power",   "switch_inline_ne")
    qsw_ne = get(_PMD.var(pm, nw), :qsw_ne, Dict()); _PMD._check_var_keys(qsw, bus_arcs_sw_ne,    "reactive power", "switch_inline_ne")
    pt     = get(_PMD.var(pm, nw), :pt,     Dict()); _PMD._check_var_keys(pt,  bus_arcs_trans,    "active power",   "transformer")
    qt     = get(_PMD.var(pm, nw), :qt,     Dict()); _PMD._check_var_keys(qt,  bus_arcs_trans,    "reactive power", "transformer")
    pt_ne  = get(_PMD.var(pm, nw), :pt_ne,  Dict()); _PMD._check_var_keys(pt,  bus_arcs_trans_ne, "active power",   "transformer_ne")
    qt_ne  = get(_PMD.var(pm, nw), :qt_ne,  Dict()); _PMD._check_var_keys(qt,  bus_arcs_trans_ne, "reactive power", "transformer_ne")


    zd = _PMD.var(pm, nw, :z_demand)
    z_shunt  = _PMD.var(pm, nw, :z_shunt)  # TODO add support for z_shunt in power balance shed
    zg = haskey(_PMD.var(pm, nw), :z_gen) ? _PMD.var(pm, nw, :z_gen) : Dict(i => 1.0 for i in _PMD.ids(pm, nw, :gen))
    zg_ne = haskey(_PMD.var(pm, nw), :z_gen_ne) ? _PMD.var(pm, nw, :z_gen_ne) : Dict(i => 1.0 for i in _PMD.ids(pm, nw, :gen_ne))
    zs = haskey(_PMD.var(pm, nw), :z_storage) ? _PMD.var(pm, nw, :z_storage) : Dict(i => 1.0 for i in _PMD.ids(pm, nw, :storage))

    Gt, Bt = _PMD._build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    # pd/qd can be NLexpressions, so cannot be vectorized
    for (idx, t) in ungrounded_terminals
        cp = _PMD.@smart_constraint(pm.model, [p, pg, ps, psw, pt],
              sum(p[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
            + sum(pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
            - sum(pg[g][t]*zg[g] for (g, conns) in bus_gens if t in conns)
            + sum(ps[s][t]*zs[s] for (s, conns) in bus_storage if t in conns)
            + sum(_PMD.ref(pm, nw, :load, d, "pd")[findfirst(isequal(t), conns)]*zd[d] for (d, conns) in bus_loads if t in conns)
            + (+vr[t] * sum(Gt[idx,jdx]*vr[u]-Bt[idx,jdx]*vi[u] for (jdx,u) in ungrounded_terminals)
               +vi[t] * sum(Gt[idx,jdx]*vi[u]+Bt[idx,jdx]*vr[u] for (jdx,u) in ungrounded_terminals))
            + sum( p_ne[a][t] for (a, conns) in bus_arcs_ne if t in conns)
            + sum(psw_ne[a][t] for (a, conns) in bus_arcs_sw_ne if t in conns)
            + sum(pt_ne[a][t] for (a, conns) in bus_arcs_trans_ne if t in conns)
            - sum(pg_ne[g][t]*zg_ne[g] for (g, conns) in bus_gens_ne if t in conns)
            ==
            0.0
        )
        push!(cstr_p, cp)

        cq = _PMD.@smart_constraint(pm.model, [q, qg, qs, qsw, qt],
              sum(q[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
            + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
            - sum(qg[g][t]*zg[g] for (g, conns) in bus_gens if t in conns)
            + sum(qs[s][t]*zs[s] for (s, conns) in bus_storage if t in conns)
            + sum(_PMD.ref(pm, nw, :load, d, "qd")[findfirst(isequal(t), conns)]*zd[d] for (d, conns) in bus_loads if t in conns)
            + (-vr[t] * sum(Gt[idx,jdx]*vi[u]+Bt[idx,jdx]*vr[u] for (jdx,u) in ungrounded_terminals)
               +vi[t] * sum(Gt[idx,jdx]*vr[u]-Bt[idx,jdx]*vi[u] for (jdx,u) in ungrounded_terminals))
            + sum(q_ne[a][t] for (a, conns) in bus_arcs_ne if t in conns)
            + sum(qsw_ne[a][t] for (a, conns) in bus_arcs_sw_ne if t in conns)
            + sum(qt_ne[a][t] for (a, conns) in bus_arcs_trans_ne if t in conns)
            - sum(qg_ne[g][t]*zg_ne[g] for (g, conns) in bus_gens_ne if t in conns)

            ==
            0.0
        )
        push!(cstr_q, cq)
    end

    con(pm, nw, :lam_kcl_r)[i] = cstr_p
    con(pm, nw, :lam_kcl_i)[i] = cstr_q

    if _INs.report_duals(pm)
        _PMD.sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        _PMD.sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end