
function variable_xe(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    if relax
        xe = var(pm, nw)[:xe] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :branch_ne)],
            base_name = "$(nw)_xe",
            lower_bound = 0,
            upper_bound = 1,
            start = _PMD.comp_start_value(ref(pm, nw, :branch_ne, i), "xe_start", i, 0.0)
        )
    else
        xe = var(pm, nw)[:xe] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :branch_ne)],
            base_name = "$(nw)_xe",
            binary = true,
            start = _PMD.comp_start_value(ref(pm, nw, :branch_ne, i), "xe_start", i, 0.0)
        )
    end
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :branch_ne, :xe, ids(pm, nw, :branch_ne), xe)
end

function variable_ue(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    if relax
        ue = var(pm, nw)[:ue] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :gen_ne)],
            base_name = "$(nw)_ue",
            lower_bound = 0,
            upper_bound = 1,
            start = _PMD.comp_start_value(ref(pm, nw, :gen_ne, i), "ue_start", i, 0.0)
        )
    else
        ue = var(pm, nw)[:ue] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :gen_ne)],
            base_name = "$(nw)_ue",
            binary = true,
            start = _PMD.comp_start_value(ref(pm, nw, :gen_ne, i), "ue_start", i, 0.0)
        )
    end
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen_ne, :ue, ids(pm, nw, :gen_ne), ue)
end

function variable_he(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    if relax
        he = var(pm, nw)[:he] = JuMP.@variable(pm.model,
            [i in ref(pm, nw, :branch_harden)],
            base_name = "$(nw)_he",
            lower_bound = 0,
            upper_bound = 1,
            start = _PMD.comp_start_value(ref(pm, nw, :branch, i), "he_start", i, 0.0)
        )
    else
        he = var(pm, nw)[:he] = JuMP.@variable(pm.model,
            [i in ref(pm, nw, :branch_harden)],
            base_name = "$(nw)_he",
            binary = true,
            start = _PMD.comp_start_value(ref(pm, nw, :branch, i), "he_start", i, 0.0)
        )
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :branch, :he, ref(pm, nw, :branch_harden), he)
end

function variable_te(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    if relax
        te = var(pm, nw)[:te] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :switch_inline_ne)],
            base_name = "$(nw)_te",
            lower_bound = 0,
            upper_bound = 1,
            start = _PMD.comp_start_value(ref(pm, nw, :switch_inline_ne, i), "te_start", i, 0.0)
        )
    else
        te = var(pm, nw)[:te] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :switch_inline_ne)],
            base_name = "$(nw)_te",
            binary = true,
            start = _PMD.comp_start_value(ref(pm, nw, :switch_inline_ne, i), "te_start", i, 0.0)
        )
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :switch_inline_ne, :te, ids(pm, nw, :switch_inline_ne), te)
end

function variable_xe_s(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    if relax
        xe_s = var(pm, nw)[:xe_s] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :branch_ne)],
            base_name = "$(nw)_xe_s_ne",
            lower_bound = 0,
            upper_bound = 1,
            start = _PMD.comp_start_value(ref(pm, nw, :branch_ne, i), "xe_start", i, 0.0)
        )
    else
        xe_s = var(pm, nw)[:xe_s] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :branch_ne)],
            base_name = "$(nw)_xe_s_ne",
            binary = true,
            start = _PMD.comp_start_value(ref(pm, nw, :branch_ne, i), "xe_start", i, 0.0)
        )
    end
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :branch_ne, :xe_s, ids(pm, nw, :branch_ne), xe_s)
end

function variable_ze_s(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    if relax
        ze_s = var(pm, nw)[:ze_s] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :branch)],
            base_name = "$(nw)_ze_s",
            lower_bound = 0,
            upper_bound = 1,
            start = _PMD.comp_start_value(ref(pm, nw, :branch, i), "ze_start", i, 0.0)
        )
        ze_s_xfr = var(pm, nw)[:ze_s_xfr] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :transformer)],
            base_name = "$(nw)_ze_s_xfr",
            lower_bound = 0,
            upper_bound = 1,
            start = _PMD.comp_start_value(ref(pm, nw, :transformer, i), "ze_start", i, 0.0)
        )
    else
        ze_s = var(pm, nw)[:ze_s] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :branch)],
            base_name = "$(nw)_ze_s",
            binary = true,
            start = _PMD.comp_start_value(ref(pm, nw, :branch, i), "ze_start", i, 0.0)
        )
        ze_s_xfr = var(pm, nw)[:ze_s_xfr] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :transformer)],
            base_name = "$(nw)_ze_s_xfr",
            binary = true,
            start = _PMD.comp_start_value(ref(pm, nw, :transformer, i), "ze_start", i, 0.0)
        )
    end
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :branch, :ze_s, ids(pm, nw, :branch), ze_s)
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :transformer, :ze_s_xfr, ids(pm, nw, :transformer), ze_s_xfr)
end


function variable_he_s(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    if relax
        he_s = var(pm, nw)[:he_s] = JuMP.@variable(pm.model,
            [i in ref(pm, nw, :branch_harden)],
            base_name = "$(nw)_he_s",
            lower_bound = 0,
            upper_bound = 1,
            start = _PMD.comp_start_value(ref(pm, nw, :branch, i), "he_start", i, 0.0)
        )
    else
        he_s = var(pm, nw)[:he_s] = JuMP.@variable(pm.model,
            [i in ref(pm, nw, :branch_harden)],
            base_name = "$(nw)_he_s",
            binary = true,
            start = _PMD.comp_start_value(ref(pm, nw, :branch, i), "he_start", i, 0.0)
        )
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :branch, :he_s, ref(pm, nw, :branch_harden), he_s)
end

function variable_ue_s(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    if relax
        ue_s = var(pm, nw)[:ue_s] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :gen_ne)],
            base_name = "$(nw)_ue_s",
            lower_bound = 0,
            upper_bound = 1,
            start = _PMD.comp_start_value(ref(pm, nw, :gen_ne, i), "ue_start", i, 0.0)
        )
    else
        ue_s = var(pm, nw)[:ue_s] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :gen_ne)],
            base_name = "$(nw)_ue_s",
            binary = true,
            start = _PMD.comp_start_value(ref(pm, nw, :gen_ne, i), "ue_start", i, 0.0)
        )
    end
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen_ne, :ue_s, ids(pm, nw, :gen_ne), ue_s)
end


"switch_inline_ne state (open/close) variables"
function variable_mc_switch_inline_ne_state(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, report::Bool=true, relax::Bool=false)
    if relax
        state = var(pm, nw)[:switch_inline_ne_state] = JuMP.@variable(
            pm.model,
            [l in ids(pm, nw, :switch_inline_ne)],
            base_name = "$(nw)_switch_inline_ne_state_$(l)",
            lower_bound = 0,
            upper_bound = 1,
            start = _PMD.comp_start_value(ref(pm, nw, :switch_inline_ne, l), "state_start", 0)
        )
    else
        state = var(pm, nw)[:switch_inline_ne_state] = JuMP.@variable(
            pm.model,
            [l in ids(pm, nw, :switch_inline_ne)],
            base_name = "$(nw)_switch_inline_ne_state_$(l)",
            binary = true,
            start = _PMD.comp_start_value(ref(pm, nw, :switch_inline_ne, l), "state_start", 0)
        )
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :switch_inline_ne, :switch_inline_ne_state, ids(pm, nw, :switch_inline_ne), state)
end


"branch_ne flow variables"
function variable_mc_branch_ne_power(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_branch_ne_power_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_branch_ne_power_imaginary(pm; nw=nw, bounded=bounded, report=report)
end

"variable: `p[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_branch_ne_power_real(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l, i, j) => connections for (bus, entry) in ref(pm, nw, :bus_arcs_conns_branch_ne) for ((l, i, j), connections) in entry)
    p = var(pm, nw)[:p_ne] = Dict((l, i, j) => JuMP.@variable(pm.model,
        [c in connections[(l, i, j)]], base_name = "$(nw)_p_ne_$((l,i,j))",
        start = _PMD.comp_start_value(ref(pm, nw, :branch_ne, l), "p_start", c, 0.0)
    ) for (l, i, j) in ref(pm, nw, :arcs_branch_ne)
    )

    if bounded
        for (l, i, j) in ref(pm, nw, :arcs_branch_ne)
            smax = _PMD._calc_branch_power_max(ref(pm, nw, :branch_ne, l), ref(pm, nw, :bus, i))
            for (idx, c) in enumerate(connections[(l, i, j)])
                JuMP.set_upper_bound(p[(l, i, j)][c], smax[idx])
                JuMP.set_lower_bound(p[(l, i, j)][c], -smax[idx])
            end
        end
    end

    for (l, branch) in ref(pm, nw, :branch_ne)
        if haskey(branch, "pf_start")
            f_idx = (l, branch["f_bus"], branch["t_bus"])
            for (idx, c) in enumerate(connections[f_idx])
                JuMP.set_start_value(p[f_idx][c], branch["pf_start"][idx])
            end
        end
        if haskey(branch, "pt_start")
            t_idx = (l, branch["t_bus"], branch["f_bus"])
            for (idx, c) in enumerate(connections[t_idx])
                JuMP.set_start_value(p[t_idx][c], branch["pt_start"][idx])
            end
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :branch_ne, :pf_ne, :pt_ne, ref(pm, nw, :arcs_branch_ne_from), ref(pm, nw, :arcs_branch_ne_to), p)
end


"variable: `q[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_branch_ne_power_imaginary(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l, i, j) => connections for (bus, entry) in ref(pm, nw, :bus_arcs_conns_branch_ne) for ((l, i, j), connections) in entry)
    q = var(pm, nw)[:q_ne] = Dict((l, i, j) => JuMP.@variable(pm.model,
        [c in connections[(l, i, j)]], base_name = "$(nw)_q_ne_$((l,i,j))",
        start = _PMD.comp_start_value(ref(pm, nw, :branch_ne, l), "q_start", c, 0.0)
    ) for (l, i, j) in ref(pm, nw, :arcs_branch_ne)
    )

    if bounded
        for (l, i, j) in ref(pm, nw, :arcs_branch_ne)
            smax = _PMD._calc_branch_power_max(ref(pm, nw, :branch_ne, l), ref(pm, nw, :bus, i))
            for (idx, c) in enumerate(connections[(l, i, j)])
                JuMP.set_upper_bound(q[(l, i, j)][c], smax[idx])
                JuMP.set_lower_bound(q[(l, i, j)][c], -smax[idx])
            end
        end
    end

    for (l, branch) in ref(pm, nw, :branch_ne)
        if haskey(branch, "qf_start")
            f_idx = (l, branch["f_bus"], branch["t_bus"])
            for (idx, c) in enumerate(connections[f_idx])
                JuMP.set_start_value(q[f_idx][c], branch["qf_start"][idx])
            end
        end
        if haskey(branch, "qt_start")
            t_idx = (l, branch["t_bus"], branch["f_bus"])
            for (idx, c) in enumerate(connections[t_idx])
                JuMP.set_start_value(q[t_idx][c], branch["qt_start"][idx])
            end
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :branch_ne, :qf_ne, :qt_ne, ref(pm, nw, :arcs_branch_ne_from), ref(pm, nw, :arcs_branch_ne_to), q)
end



function variable_mc_switch_inline_ne_power(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_switch_inline_ne_power_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_switch_inline_ne_power_imaginary(pm; nw=nw, bounded=bounded, report=report)
end


""
function variable_mc_switch_inline_ne_power_real(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l, i, j) => connections for (bus, entry) in ref(pm, nw, :bus_arcs_conns_switch_inline_ne) for ((l, i, j), connections) in entry)
    psw = Dict((l, i, j) => JuMP.@variable(pm.model,
        [c in connections[(l, i, j)]], base_name = "$(nw)_psw_inline_ne_$((l,i,j))",
        start = _PMD.comp_start_value(ref(pm, nw, :switch_inline_ne, l), "psw_start", c, 0.0)
    ) for (l, i, j) in ref(pm, nw, :arcs_switch_inline_ne)
    )

    if bounded
        for (l, i, j) in ref(pm, nw, :arcs_switch_inline_ne)
            smax = _PMD._calc_branch_power_max(ref(pm, nw, :switch_inline_ne, l), ref(pm, nw, :bus, i))
            for (idx, c) in enumerate(connections[(l, i, j)])
                JuMP.set_upper_bound(psw[(l, i, j)][c], smax[idx])
                JuMP.set_lower_bound(psw[(l, i, j)][c], -smax[idx])
            end
        end
    end

    # this explicit type erasure is necessary
    psw_expr = Dict{Any,Any}((l, i, j) => psw[(l, i, j)] for (l, i, j) in ref(pm, nw, :arcs_switch_inline_ne_from))
    psw_expr = merge(psw_expr, Dict((l, j, i) => -1.0 .* psw[(l, i, j)] for (l, i, j) in ref(pm, nw, :arcs_switch_inline_ne_from)))

    # This is needed to get around error: "unexpected affine expression in nlconstraint"
    psw_auxes = Dict{Any,Any}(
        (l, i, j) => JuMP.@variable(
            pm.model, [c in connections[(l, i, j)]],
            base_name = "$(nw)_psw_inline_ne_aux_$((l,i,j))"
        ) for (l, i, j) in ref(pm, nw, :arcs_switch_inline_ne)
    )
    for ((l, i, j), psw_aux) in psw_auxes
        for (idx, c) in enumerate(connections[(l, i, j)])
            JuMP.@constraint(pm.model, psw_expr[(l, i, j)][c] == psw_aux[c])
        end
    end

    var(pm, nw)[:psw_inline_ne] = psw_auxes

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :switch_inline_ne, :pf_ne, :pt_ne, ref(pm, nw, :arcs_switch_inline_ne_from), ref(pm, nw, :arcs_switch_inline_ne_to), psw_expr)
end


""
function variable_mc_switch_inline_ne_power_imaginary(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l, i, j) => connections for (bus, entry) in ref(pm, nw, :bus_arcs_conns_switch_inline_ne) for ((l, i, j), connections) in entry)
    qsw = Dict((l, i, j) => JuMP.@variable(pm.model,
        [c in connections[(l, i, j)]], base_name = "$(nw)_qsw_inline_ne_$((l,i,j))",
        start = _PMD.comp_start_value(ref(pm, nw, :switch_inline_ne, l), "qsw_start", c, 0.0)
    ) for (l, i, j) in ref(pm, nw, :arcs_switch_inline_ne)
    )

    if bounded
        for (l, i, j) in ref(pm, nw, :arcs_switch_inline_ne)
            smax = _PMD._calc_branch_power_max(ref(pm, nw, :switch_inline_ne, l), ref(pm, nw, :bus, i))
            for (idx, c) in enumerate(connections[(l, i, j)])
                JuMP.set_upper_bound(qsw[(l, i, j)][c], smax[idx])
                JuMP.set_lower_bound(qsw[(l, i, j)][c], -smax[idx])
            end
        end
    end

    # this explicit type erasure is necessary
    qsw_expr = Dict{Any,Any}((l, i, j) => qsw[(l, i, j)] for (l, i, j) in ref(pm, nw, :arcs_switch_inline_ne_from))
    qsw_expr = merge(qsw_expr, Dict((l, j, i) => -1.0 * qsw[(l, i, j)] for (l, i, j) in ref(pm, nw, :arcs_switch_inline_ne_from)))

    # This is needed to get around error: "unexpected affine expression in nlconstraint"
    qsw_auxes = Dict{Any,Any}(
        (l, i, j) => JuMP.@variable(
            pm.model, [c in connections[(l, i, j)]],
            base_name = "$(nw)_qsw_aux_inline_ne_$((l,i,j))"
        ) for (l, i, j) in ref(pm, nw, :arcs_switch_inline_ne)
    )
    for ((l, i, j), qsw_aux) in qsw_auxes
        for (idx, c) in enumerate(connections[(l, i, j)])
            JuMP.@constraint(pm.model, qsw_expr[(l, i, j)][c] == qsw_aux[c])
        end
    end

    var(pm, nw)[:qsw_inline_ne] = qsw_auxes

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :switch_inline_ne, :qf_ne, :qt_ne, ref(pm, nw, :arcs_switch_inline_ne_from), ref(pm, nw, :arcs_switch_inline_ne_to), qsw_expr)
end


"Creates variables for both `active` and `reactive` power flow at each transformer."
function variable_mc_transformer_ne_power(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_transformer_ne_power_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_transformer_ne_power_imaginary(pm; nw=nw, bounded=bounded, report=report)
end


"Create variables for the active power flowing into all transformer windings."
function variable_mc_transformer_ne_power_real(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l, i, j) => connections for (bus, entry) in ref(pm, nw, :bus_arcs_conns_transformer_ne) for ((l, i, j), connections) in entry)
    pt = var(pm, nw)[:pt_ne] = Dict((l, i, j) => JuMP.@variable(pm.model,
        [c in connections[(l, i, j)]],
        base_name = "$(nw)_pt_ne_$((l,i,j))",
        start = 0.0,
    ) for (l, i, j) in ref(pm, nw, :arcs_transformer_ne)
    )

    if bounded
        for arc in ref(pm, nw, :arcs_transformer_ne_from)
            (l, i, j) = arc
            rate_a_fr, rate_a_to = _PMD._calc_transformer_power_ub_frto(ref(pm, nw, :transformer_ne, l), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))
            for (idx, (fc, tc)) in enumerate(zip(connections[(l, i, j)], connections[(l, j, i)]))
                JuMP.set_lower_bound(pt[(l, i, j)][fc], -rate_a_fr[idx])
                JuMP.set_upper_bound(pt[(l, i, j)][fc], rate_a_fr[idx])
                JuMP.set_lower_bound(pt[(l, j, i)][tc], -rate_a_to[idx])
                JuMP.set_upper_bound(pt[(l, j, i)][tc], rate_a_to[idx])
            end
        end
    end

    for (l, transformer) in ref(pm, nw, :transformer_ne)
        if haskey(transformer, "pf_start")
            f_idx = (l, transformer["f_bus"], transformer["t_bus"])
            for (idx, c) in enumerate(connections[f_idx])
                JuMP.set_start_value(pt[f_idx][c], transformer["pf_start"][idx])
            end
        end
        if haskey(transformer, "pt_start")
            t_idx = (l, transformer["t_bus"], transformer["f_bus"])
            for (idx, c) in enumerate(connections[t_idx])
                JuMP.set_start_value(pt[t_idx][c], transformer["pt_start"][idx])
            end
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :transformer_ne, :pf_ne, :pt_ne, ref(pm, nw, :arcs_transformer_ne_from), ref(pm, nw, :arcs_transformer_ne_to), pt)
end


"Create variables for the reactive power flowing into all transformer windings."
function variable_mc_transformer_ne_power_imaginary(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l, i, j) => connections for (bus, entry) in ref(pm, nw, :bus_arcs_conns_transformer_ne) for ((l, i, j), connections) in entry)
    qt = var(pm, nw)[:qt_ne] = Dict((l, i, j) => JuMP.@variable(pm.model,
        [c in connections[(l, i, j)]], base_name = "$(nw)_qt_ne_$((l,i,j))",
        start = 0.0
    ) for (l, i, j) in ref(pm, nw, :arcs_transformer_ne)
    )

    if bounded
        for arc in ref(pm, nw, :arcs_transformer_ne_from)
            (l, i, j) = arc
            rate_a_fr, rate_a_to = _PMD._calc_transformer_power_ub_frto(ref(pm, nw, :transformer_ne, l), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))

            for (idx, (fc, tc)) in enumerate(zip(connections[(l, i, j)], connections[(l, j, i)]))
                JuMP.set_lower_bound(qt[(l, i, j)][fc], -rate_a_fr[idx])
                JuMP.set_upper_bound(qt[(l, i, j)][fc], rate_a_fr[idx])
                JuMP.set_lower_bound(qt[(l, j, i)][tc], -rate_a_to[idx])
                JuMP.set_upper_bound(qt[(l, j, i)][tc], rate_a_to[idx])
            end
        end
    end

    for (l, transformer) in ref(pm, nw, :transformer_ne)
        if haskey(transformer, "qf_start")
            f_idx = (l, transformer["f_bus"], transformer["t_bus"])
            for (idx, fc) in enumerate(connections[f_idx])
                JuMP.set_start_value(qt[f_idx][fc], transformer["qf_start"][idx])
            end
        end
        if haskey(transformer, "qt_start")
            t_idx = (l, transformer["t_bus"], transformer["f_bus"])
            for (idx, tc) in enumerate(connections[t_idx])
                JuMP.set_start_value(qt[t_idx][tc], transformer["qt_start"][idx])
            end
        end
    end

    report && _IM.sol_component_value_edge(pm, pmd_it_sym, nw, :transformer_ne, :qf_ne, :qt_ne, ref(pm, nw, :arcs_transformer_ne_from), ref(pm, nw, :arcs_transformer_ne_to), qt)
end

"Create variables for generator expansion status"
function variable_mc_gen_ne_indicator(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    if !relax
        z_gen_ne = var(pm, nw)[:z_gen_ne] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :gen_ne)], base_name = "$(nw)_z_gen_ne",
            binary = true,
            start = _PMD.comp_start_value(ref(pm, nw, :gen_ne, i), "z_gen_start", 1.0)
        )
    else
        z_gen_ne = var(pm, nw)[:z_gen_ne] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :gen_ne)], base_name = "$(nw)_z_gen_ne",
            lower_bound = 0,
            upper_bound = 1,
            start = _PMD.comp_start_value(ref(pm, nw, :gen_ne, i), "z_gen_start", 1.0)
        )
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen_ne, :gen_status, ids(pm, nw, :gen_ne), z_gen_ne)
end

""
function variable_mc_generator_ne_power_on_off(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_generator_ne_power_real_on_off(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_generator_ne_power_imaginary_on_off(pm; nw=nw, bounded=bounded, report=report)
end


""
function variable_mc_generator_ne_power_real_on_off(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => gen["connections"] for (i, gen) in ref(pm, nw, :gen_ne))
    pg = var(pm, nw)[:pg_ne] = Dict(i => JuMP.@variable(pm.model,
        [c in connections[i]], base_name = "$(nw)_pg_ne_$(i)",
        start = _PMD.comp_start_value(ref(pm, nw, :gen_ne, i), ["pg_start", "pg", "pmin"], c, 0.0)
    ) for i in ids(pm, nw, :gen_ne))

    if bounded
        for (i, gen) in ref(pm, nw, :gen_ne)
            if haskey(gen, "pmin")
                for (idx, c) in enumerate(connections[i])
                    JuMP.set_lower_bound(pg[i][c], min(gen["pmin"][idx], 0.0))
                end
            end

            if haskey(gen, "pmax")
                for (idx, c) in enumerate(connections[i])
                    JuMP.set_upper_bound(pg[i][c], max(gen["pmax"][idx], 0.0))
                end
            end
        end
    end

    var(pm, nw)[:pg_bus_ne] = Dict{Int,Any}()

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen_ne, :pg_ne, ids(pm, nw, :gen_ne), pg)
end


""
function variable_mc_generator_ne_power_imaginary_on_off(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => gen["connections"] for (i, gen) in ref(pm, nw, :gen_ne))
    qg = var(pm, nw)[:qg_ne] = Dict(i => JuMP.@variable(pm.model,
        [c in connections[i]], base_name = "$(nw)_qg_ne_$(i)",
        start = _PMD.comp_start_value(ref(pm, nw, :gen_ne, i), ["qg_start", "qg", "qmin"], c, 0.0)
    ) for i in ids(pm, nw, :gen_ne))

    if bounded
        for (i, gen) in ref(pm, nw, :gen_ne)
            if haskey(gen, "qmin")
                for (idx, c) in enumerate(connections[i])
                    JuMP.set_lower_bound(qg[i][c], min(gen["qmin"][idx], 0.0))
                end
            end

            if haskey(gen, "qmax")
                for (idx, c) in enumerate(connections[i])
                    JuMP.set_upper_bound(qg[i][c], max(gen["qmax"][idx], 0.0))
                end
            end
        end
    end

    var(pm, nw)[:qg_bus_ne] = Dict{Int,Any}()

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :qg_ne, ids(pm, nw, :gen_ne), qg)
end
