"power balance constraint with line shunts and transformers for load shed problem, DCP formulation"
function constraint_mc_power_balance_shed_ne(pm::_PMD.AbstractUnbalancedDCPModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
                                          bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
                                          bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}},
                                          bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}}, bus_arcs_ne::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
                                          bus_arcs_sw_ne::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans_ne::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
                                          bus_gens_ne::Vector{Tuple{Int,Vector{Int}}})

    p       = get(_PMD.var(pm, nw), :p,      Dict()); _PMD._check_var_keys(p, bus_arcs, "active power", "branch")
    pg      = get(_PMD.var(pm, nw), :pg_bus, Dict()); _PMD._check_var_keys(pg, bus_gens, "active power", "generator")
    ps      = get(_PMD.var(pm, nw), :ps,     Dict()); _PMD._check_var_keys(ps, bus_storage, "active power", "storage")
    psw     = get(_PMD.var(pm, nw), :psw,    Dict()); _PMD._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    pt      = get(_PMD.var(pm, nw), :pt,     Dict()); _PMD._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    p_ne    = get(_PMD.var(pm, nw), :p_ne,   Dict()); _PMD._check_var_keys(p, bus_arcs_ne, "active power", "branch_ne")
    psw_ne  = get(_PMD.var(pm, nw), :psw,    Dict()); _PMD._check_var_keys(psw, bus_arcs_sw_ne, "active power", "switch_inline_ne")
    pt_ne   = get(_PMD.var(pm, nw), :pt_ne,  Dict()); _PMD._check_var_keys(pt, bus_arcs_trans_ne, "active power", "transformer_ne")
    pg_ne   = get(_PMD.var(pm, nw), :pg_ne, Dict()); _PMD._check_var_keys(pg_ne, bus_gens_ne, "active power", "generator_ne")


    z_demand = _PMD.var(pm, nw, :z_demand)
    z_gen = haskey(_PMD.var(pm, nw), :z_gen) ? _PMD.var(pm, nw, :z_gen) : Dict(i => 1.0 for i in _PMD.ids(pm, nw, :gen))
    z_gen_ne = haskey(_PMD.var(pm, nw), :z_gen_ne) ? _PMD.var(pm, nw, :z_gen_ne) : Dict(i => 1.0 for i in _PMD.ids(pm, nw, :gen_ne))
    z_storage = haskey(_PMD.var(pm, nw), :z_storage) ? _PMD.var(pm, nw, :z_storage) : Dict(i => 1.0 for i in _PMD.ids(pm, nw, :storage))
    z_shunt  = haskey(_PMD.var(pm, nw), :z_shunt) ? _PMD.var(pm, nw, :z_shunt) : Dict(i => 1.0 for i in _PMD.ids(pm, nw, :shunt))

    Gt, Bt = _PMD._build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx,t) in ungrounded_terminals
        cp = JuMP.@constraint(pm.model,
              sum(p[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(psw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
            + sum(pt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
            + sum(ref(pm, nw, :load, d, "pd")[findfirst(isequal(t), conns)]*z_demand[d] for (d, conns) in bus_loads if t in conns)
            - sum(pg[g][t]*z_gen[g] for (g, conns) in bus_gens if t in conns)
            - sum(ps[s][t]*z_storage[s] for (s, conns) in bus_storage if t in conns)
            + sum(LinearAlgebra.diag(Gt)[idx]*z_shunt[sh] for (sh, conns) in bus_shunts if t in conns)
            + sum(p_ne[a][t] for (a, conns) in bus_arcs_ne if t in conns)
            + sum(psw_ne[a_sw][t] for (a_sw, conns) in bus_arcs_sw_ne if t in conns)
            + sum(pt_ne[a_trans][t] for (a_trans, conns) in bus_arcs_trans_ne if t in conns)
            - sum(pg_ne[g][t]*z_gen_ne[g] for (g, conns) in bus_gens_ne if t in conns)
            == 0
        )
        push!(cstr_p, cp)
    end

    con(pm, nw, :lam_kcl_r)[i] = cstr_p
    con(pm, nw, :lam_kcl_i)[i] = []

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
    end
end
