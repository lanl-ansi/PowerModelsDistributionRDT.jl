
function constraint_mc_power_ne_balance_shed(pm::_PMD.AbstractUnbalancedACRModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool},
                                             bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
                                             bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}},
                                             bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}},
                                             bus_shunts::Vector{Tuple{Int,Vector{Int}}}, bus_arcs_ne::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
                                             bus_arcs_sw_ne::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans_ne::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
                                             bus_gens_ne::Vector{Tuple{Int,Vector{Int}}})
    vr     = var(pm, nw, :vr, i)
    vi     = var(pm, nw, :vi, i)
    p      = get(var(pm, nw), :p,      Dict()); _PMD._check_var_keys(p,   bus_arcs,          "active power",   "branch")
    q      = get(var(pm, nw), :q,      Dict()); _PMD._check_var_keys(q,   bus_arcs,          "reactive power", "branch")
    p_ne   = get(var(pm, nw), :p_ne,   Dict()); _PMD._check_var_keys(p,   bus_arcs_ne,       "active power",   "branch_ne")
    q_ne   = get(var(pm, nw), :q_ne,   Dict()); _PMD._check_var_keys(q,   bus_arcs_ne,       "reactive power", "branch_ne")
    pg     = get(var(pm, nw), :pg,     Dict()); _PMD._check_var_keys(pg,  bus_gens,          "active power",   "generator")
    qg     = get(var(pm, nw), :qg,     Dict()); _PMD._check_var_keys(qg,  bus_gens,          "reactive power", "generator")
    pg_ne  = get(var(pm, nw), :pg_ne,  Dict()); _PMD._check_var_keys(pg,  bus_gens_ne,       "active power",   "generator_ne")
    qg_ne  = get(var(pm, nw), :qg_ne,  Dict()); _PMD._check_var_keys(qg,  bus_gens_ne,       "reactive power", "generator_ne")
    ps     = get(var(pm, nw), :ps,     Dict()); _PMD._check_var_keys(ps,  bus_storage,       "active power",   "storage")
    qs     = get(var(pm, nw), :qs,     Dict()); _PMD._check_var_keys(qs,  bus_storage,       "reactive power", "storage")
    psw    = get(var(pm, nw), :psw,    Dict()); _PMD._check_var_keys(psw, bus_arcs_sw,       "active power",   "switch")
    qsw    = get(var(pm, nw), :qsw,    Dict()); _PMD._check_var_keys(qsw, bus_arcs_sw,       "reactive power", "switch")
    psw_ne = get(var(pm, nw), :psw_inline_ne, Dict()); _PMD._check_var_keys(psw, bus_arcs_sw_ne,    "active power",   "switch_inline_ne")
    qsw_ne = get(var(pm, nw), :qsw_inline_ne, Dict()); _PMD._check_var_keys(qsw, bus_arcs_sw_ne,    "reactive power", "switch_inline_ne")
    pt     = get(var(pm, nw), :pt,     Dict()); _PMD._check_var_keys(pt,  bus_arcs_trans,    "active power",   "transformer")
    qt     = get(var(pm, nw), :qt,     Dict()); _PMD._check_var_keys(qt,  bus_arcs_trans,    "reactive power", "transformer")
    pt_ne  = get(var(pm, nw), :pt_ne,  Dict()); _PMD._check_var_keys(pt,  bus_arcs_trans_ne, "active power",   "transformer_ne")
    qt_ne  = get(var(pm, nw), :qt_ne,  Dict()); _PMD._check_var_keys(qt,  bus_arcs_trans_ne, "reactive power", "transformer_ne")


    zd = var(pm, nw, :z_demand)
    z_shunt  = var(pm, nw, :z_shunt)  # TODO add support for z_shunt in power balance shed
    zg = haskey(var(pm, nw), :z_gen) ? var(pm, nw, :z_gen) : Dict(i => 1.0 for i in ids(pm, nw, :gen))
    zg_ne = haskey(var(pm, nw), :z_gen_ne) ? var(pm, nw, :z_gen_ne) : Dict(i => 1.0 for i in ids(pm, nw, :gen_ne))
    zs = haskey(var(pm, nw), :z_storage) ? var(pm, nw, :z_storage) : Dict(i => 1.0 for i in ids(pm, nw, :storage))

    Gt, Bt = _PMD._build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    cstr_p = []
    cstr_q = []

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    # pd/qd can be NLexpressions, so cannot be vectorized
    for (idx, t) in ungrounded_terminals
        cp = JuMP.@constraint(pm.model,
              sum(p[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(psw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
            + sum(pt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
            - sum(pg[g][t]*zg[g] for (g, conns) in bus_gens if t in conns)
            + sum(ps[s][t]*zs[s] for (s, conns) in bus_storage if t in conns)
            + sum(ref(pm, nw, :load, d, "pd")[findfirst(isequal(t), conns)]*zd[d] for (d, conns) in bus_loads if t in conns)
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

        cq = JuMP.@constraint(pm.model,
              sum(q[a][t] for (a, conns) in bus_arcs if t in conns)
            + sum(qsw[a][t] for (a, conns) in bus_arcs_sw if t in conns)
            + sum( qt[a][t] for (a, conns) in bus_arcs_trans if t in conns)
            - sum(qg[g][t]*zg[g] for (g, conns) in bus_gens if t in conns)
            + sum(qs[s][t]*zs[s] for (s, conns) in bus_storage if t in conns)
            + sum(ref(pm, nw, :load, d, "qd")[findfirst(isequal(t), conns)]*zd[d] for (d, conns) in bus_loads if t in conns)
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

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end




@doc raw"""
    constraint_mc_ampacity_from_damaged(pm::AbstractUnbalancedRectangularModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on branches from-side on branches that are damaged

math```
p_{fr}^2 + q_{fr}^2 \leq (vr_{fr}^2 + vi_{fr}^2) i_{max}^2 * he_s
```
"""
function constraint_mc_ampacity_from_damaged(pm::_PMD.AbstractUnbalancedRectangularModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    p_fr = [var(pm, nw, :p, f_idx)[c] for c in f_connections]
    q_fr = [var(pm, nw, :q, f_idx)[c] for c in f_connections]
    vr_fr = [var(pm, nw, :vr, f_idx[2])[c] for c in f_connections]
    vi_fr = [var(pm, nw, :vi, f_idx[2])[c] for c in f_connections]
    he_s  = var(pm, nw, :he_s, f_idx[1])

    # TODO: maybe introduce an auxillary varaible v_sqr = vr_fr[idx]^2 + vi_fr[idx]^2, and do exact McCormick on v_sqr * he_s (and use @constraint)
    con(pm, nw, :mu_cm_branch)[f_idx] = mu_cm_fr = [JuMP.@constraint(pm.model, p_fr[idx]^2 + q_fr[idx]^2 .<= (vr_fr[idx]^2 + vi_fr[idx]^2) * c_rating[idx]^2 * he_s) for idx in f_connections]

    if _IM.report_duals(pm)
        sol(pm, nw, :branch, f_idx[1])[:mu_cm_fr] = mu_cm_fr
    end

    nothing
end


@doc raw"""
    constraint_mc_ampacity_to_damaged(pm::AbstractUnbalancedRectangularModels, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on branches to-side that are damaged

math```
p_{to}^2 + q_{to}^2 \leq (vr_{to}^2 + vi_{to}^2) i_{max}^2 * he_s
```
"""
function constraint_mc_ampacity_to_damaged(pm::_PMD.AbstractUnbalancedRectangularModels, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    p_to = [var(pm, nw, :p, t_idx)[c] for c in t_connections]
    q_to = [var(pm, nw, :q, t_idx)[c] for c in t_connections]
    vr_to = [var(pm, nw, :vr, t_idx[2])[c] for c in t_connections]
    vi_to = [var(pm, nw, :vi, t_idx[2])[c] for c in t_connections]
    he_s  = var(pm, nw, :he_s, t_idx[1])

    # TODO: maybe introduce an auxillary varaible v_sqr = vr_to[idx]^2 + vi_to[idx]^2, and do exact McCormick on v_sqr * he_s (and use @constraint)
    con(pm, nw, :mu_cm_branch)[t_idx] = mu_cm_to = [JuMP.@constraint(pm.model, p_to[idx]^2 + q_to[idx]^2 .<= (vr_to[idx]^2 + vi_to[idx]^2) * c_rating[idx]^2 * he_s) for idx in t_connections]

    if _IM.report_duals(pm)
        sol(pm, nw, :branch, t_idx[1])[:mu_cm_to] = mu_cm_to
    end

    nothing
end


@doc raw"""
    constraint_mc_ampacity_from_ne(pm::AbstractUnbalancedRectangularModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on branches from-side on branches that are ne

math```
p_{fr}^2 + q_{fr}^2 \leq (vr_{fr}^2 + vi_{fr}^2) i_{max}^2 * xe_s
```
"""
function constraint_mc_ampacity_from_ne(pm::_PMD.AbstractUnbalancedRectangularModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    p_fr = [var(pm, nw, :p_ne, f_idx)[c] for c in f_connections]
    q_fr = [var(pm, nw, :q_ne, f_idx)[c] for c in f_connections]
    vr_fr = [var(pm, nw, :vr, f_idx[2])[c] for c in f_connections]
    vi_fr = [var(pm, nw, :vi, f_idx[2])[c] for c in f_connections]
    xe_s  = var(pm, nw, :xe_s, f_idx[1])

    # TODO: maybe introduce an auxillary varaible v_sqr = vr_fr[idx]^2 + vi_fr[idx]^2, and do exact McCormick on v_sqr * he_s (and use @constraint)
    con(pm, nw, :mu_cm_branch_ne)[f_idx] = mu_cm_fr = [JuMP.@constraint(pm.model, p_fr[idx]^2 + q_fr[idx]^2 .<= (vr_fr[idx]^2 + vi_fr[idx]^2) * c_rating[idx]^2 * xe_s) for idx in f_connections]

    if _IM.report_duals(pm)
        sol(pm, nw, :branch_ne, f_idx[1])[:mu_cm_fr_ne] = mu_cm_fr
    end

    nothing
end


@doc raw"""
    constraint_mc_ampacity_to_ne(pm::AbstractUnbalancedRectangularModels, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on branches to-side that are ne

math```
p_{to}^2 + q_{to}^2 \leq (vr_{to}^2 + vi_{to}^2) i_{max}^2 * he_s
```
"""
function constraint_mc_ampacity_to_ne(pm::_PMD.AbstractUnbalancedRectangularModels, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    p_to = [var(pm, nw, :p_ne, t_idx)[c] for c in t_connections]
    q_to = [var(pm, nw, :q_ne, t_idx)[c] for c in t_connections]
    vr_to = [var(pm, nw, :vr, t_idx[2])[c] for c in t_connections]
    vi_to = [var(pm, nw, :vi, t_idx[2])[c] for c in t_connections]
    xe_s  = var(pm, nw, :xe_s, t_idx[1])

    # TODO: maybe introduce an auxillary varaible v_sqr = vr_to[idx]^2 + vi_to[idx]^2, and do exact McCormick on v_sqr * xe_s (and use @constraint)
    con(pm, nw, :mu_cm_branch_ne)[t_idx] = mu_cm_to = [JuMP.@constraint(pm.model, p_to[idx]^2 + q_to[idx]^2 .<= (vr_to[idx]^2 + vi_to[idx]^2) * c_rating[idx]^2 * xe_s) for idx in t_connections]

    if _IM.report_duals(pm)
        sol(pm, nw, :branch_ne, t_idx[1])[:mu_cm_to_ne] = mu_cm_to
    end

    nothing
end


""
function constraint_mc_voltage_angle_difference_damaged(pm::_PMD.AbstractUnbalancedACRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, angmin::Vector{<:Real}, angmax::Vector{<:Real})
    i, f_bus, t_bus = f_idx

    vr_fr = var(pm, nw, :vr, f_bus)
    vi_fr = var(pm, nw, :vi, f_bus)
    vr_to = var(pm, nw, :vr, t_bus)
    vi_to = var(pm, nw, :vi, t_bus)
    he_s  = var(pm, nw, :he_s, i)

    #TODO: A bit lazy, but this is how PowerModels.jl does on_off on phase angle difference constraints.  If we had the absolute maximum voltage angle difference, from bound on the v variables
    for (idx, (fc,tc)) in enumerate(zip(f_connections, t_connections))
        JuMP.@constraint(pm.model, he_s * (vi_fr[fc] * vr_to[tc] .- vr_fr[fc] * vi_to[tc]) <= he_s * tan(angmax[idx]) * (vr_fr[fc] * vr_to[tc] .+ vi_fr[fc] * vi_to[tc]))
        JuMP.@constraint(pm.model, he_s * (vi_fr[fc] * vr_to[tc] .- vr_fr[fc] * vi_to[tc]) >= he_s * tan(angmin[idx]) * (vr_fr[fc] * vr_to[tc] .+ vi_fr[fc] * vi_to[tc]))
    end
end


""
function constraint_mc_voltage_angle_difference_ne(pm::_PMD.AbstractUnbalancedACRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, angmin::Vector{<:Real}, angmax::Vector{<:Real})
    i, f_bus, t_bus = f_idx

    vr_fr = var(pm, nw, :vr, f_bus)
    vi_fr = var(pm, nw, :vi, f_bus)
    vr_to = var(pm, nw, :vr, t_bus)
    vi_to = var(pm, nw, :vi, t_bus)
    xe_s  = var(pm, nw, :xe_s, i)

    #TODO: A bit lazy, but this is how PowerModels.jl does on_off on phase angle difference constraints.  If we had the absolute maximum voltage angle difference, from bound on the v variables
    for (idx, (fc,tc)) in enumerate(zip(f_connections, t_connections))
        JuMP.@constraint(pm.model, xe_s * (vi_fr[fc] * vr_to[tc] .- vr_fr[fc] * vi_to[tc]) <= xe_s * tan(angmax[idx]) * (vr_fr[fc] * vr_to[tc] .+ vi_fr[fc] * vi_to[tc]))
        JuMP.@constraint(pm.model, xe_s * (vi_fr[fc] * vr_to[tc] .- vr_fr[fc] * vi_to[tc]) >= xe_s * tan(angmin[idx]) * (vr_fr[fc] * vr_to[tc] .+ vi_fr[fc] * vi_to[tc]))
    end
end

"""
Creates Ohms constraints for damaged lines

s_fr = he * v_fr.*conj(Y*(v_fr-v_to))
s_fr = he * (vr_fr+im*vi_fr).*(G-im*B)*([vr_fr-vr_to]-im*[vi_fr-vi_to])
s_fr = he * (vr_fr+im*vi_fr).*([G*vr_fr-G*vr_to-B*vi_fr+B*vi_to]-im*[G*vi_fr-G*vi_to+B*vr_fr-B*vr_to])
"""
function constraint_mc_ohms_yt_from_damaged(pm::_PMD.AbstractUnbalancedACRModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real}, vad_min::Vector{<:Real}, vad_max::Vector{<:Real})
    p_fr  = [var(pm, nw, :p, f_idx)[t] for t in f_connections]
    q_fr  = [var(pm, nw, :q, f_idx)[t] for t in f_connections]
    vr_fr = [var(pm, nw, :vr, f_bus)[t] for t in f_connections]
    vr_to = [var(pm, nw, :vr, t_bus)[t] for t in t_connections]
    vi_fr = [var(pm, nw, :vi, f_bus)[t] for t in f_connections]
    vi_to = [var(pm, nw, :vi, t_bus)[t] for t in t_connections]
    he_s  = var(pm, nw, :he_s, f_idx[1])

    con(pm, nw, :ohms_yt)[f_idx] = [
        JuMP.@constraint(pm.model,
            p_fr .==  he_s * (vr_fr.*(G*vr_fr-G*vr_to-B*vi_fr+B*vi_to)
                     +vi_fr.*(G*vi_fr-G*vi_to+B*vr_fr-B*vr_to)
                     # shunt
                     +vr_fr.*(G_fr*vr_fr-B_fr*vi_fr)
                     +vi_fr.*(G_fr*vi_fr+B_fr*vr_fr))
        ),
        JuMP.@constraint(pm.model,
            q_fr .== he_s * (-vr_fr.*(G*vi_fr-G*vi_to+B*vr_fr-B*vr_to)
                     +vi_fr.*(G*vr_fr-G*vr_to-B*vi_fr+B*vi_to)
                     # shunt
                     -vr_fr.*(G_fr*vi_fr+B_fr*vr_fr)
                     +vi_fr.*(G_fr*vr_fr-B_fr*vi_fr))
        )
    ]
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form) for damaged lines

```
p[t_idx] ==  he * (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == he * -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_mc_ohms_yt_to_damaged(pm::_PMD.AbstractUnbalancedACRModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix, B::Matrix, G_to::Matrix, B_to::Matrix, vad_min::Vector{<:Real}, vad_max::Vector{<:Real})
    constraint_mc_ohms_yt_from_damaged(pm, nw, t_bus, f_bus, t_idx, f_idx, t_connections, f_connections, G, B, G_to, B_to, vad_min, vad_max)
end


"""
Creates Ohms constraints for ne lines

s_fr = he * v_fr.*conj(Y*(v_fr-v_to))
s_fr = he * (vr_fr+im*vi_fr).*(G-im*B)*([vr_fr-vr_to]-im*[vi_fr-vi_to])
s_fr = he * (vr_fr+im*vi_fr).*([G*vr_fr-G*vr_to-B*vi_fr+B*vi_to]-im*[G*vi_fr-G*vi_to+B*vr_fr-B*vr_to])
"""
function constraint_mc_ohms_yt_from_ne(pm::_PMD.AbstractUnbalancedACRModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix{<:Real}, B::Matrix{<:Real}, G_fr::Matrix{<:Real}, B_fr::Matrix{<:Real}, vad_min::Vector{<:Real}, vad_max::Vector{<:Real})
    p_fr  = [var(pm, nw, :p_ne, f_idx)[t] for t in f_connections]
    q_fr  = [var(pm, nw, :q_ne, f_idx)[t] for t in f_connections]
    vr_fr = [var(pm, nw, :vr, f_bus)[t] for t in f_connections]
    vr_to = [var(pm, nw, :vr, t_bus)[t] for t in t_connections]
    vi_fr = [var(pm, nw, :vi, f_bus)[t] for t in f_connections]
    vi_to = [var(pm, nw, :vi, t_bus)[t] for t in t_connections]
    xe_s  = var(pm, nw, :xe_s, f_idx[1])

    con(pm, nw, :ohms_yt)[f_idx] = [
        JuMP.@constraint(pm.model,
            p_fr .==  xe_s * (vr_fr.*(G*vr_fr-G*vr_to-B*vi_fr+B*vi_to)
                     +vi_fr.*(G*vi_fr-G*vi_to+B*vr_fr-B*vr_to)
                     # shunt
                     +vr_fr.*(G_fr*vr_fr-B_fr*vi_fr)
                     +vi_fr.*(G_fr*vi_fr+B_fr*vr_fr))
        ),
        JuMP.@constraint(pm.model,
            q_fr .== xe_s * (-vr_fr.*(G*vi_fr-G*vi_to+B*vr_fr-B*vr_to)
                     +vi_fr.*(G*vr_fr-G*vr_to-B*vi_fr+B*vi_to)
                     # shunt
                     -vr_fr.*(G_fr*vi_fr+B_fr*vr_fr)
                     +vi_fr.*(G_fr*vr_fr-B_fr*vi_fr))
        )
    ]
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form) for damaged lines

```
p[t_idx] ==  he * (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
q[t_idx] == he * -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))
```
"""
function constraint_mc_ohms_yt_to_ne(pm::_PMD.AbstractUnbalancedACRModel, nw::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, G::Matrix, B::Matrix, G_to::Matrix, B_to::Matrix, vad_min::Vector{<:Real}, vad_max::Vector{<:Real})
    constraint_mc_ohms_yt_from_ne(pm, nw, t_bus, f_bus, t_idx, f_idx, t_connections, f_connections, G, B, G_to, B_to, vad_min, vad_max)
end


@doc raw"""
    constraint_mc_switch_voltage_open_close(pm::PMD.AbstractUnbalancedACRModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_connections::Vector{Int}, t_connections::Vector{Int})

nonlinear switch power on/off constraint for ac-rect form

```math
\begin{align}
& \\
&
\end{align}
```
"""
function constraint_mc_switch_inline_ne_voltage_open_close(pm::_PMD.AbstractUnbalancedACRModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_connections::Vector{Int}, t_connections::Vector{Int})
    vr_fr = var(pm, nw, :vr, f_bus)
    vr_to = var(pm, nw, :vr, t_bus)
    vi_fr = var(pm, nw, :vi, f_bus)
    vi_to = var(pm, nw, :vi, t_bus)

    f_bus = ref(pm, nw, :bus, f_bus)
    t_bus = ref(pm, nw, :bus, t_bus)

    f_vmin = f_bus["vmin"][[findfirst(isequal(c), f_bus["terminals"]) for c in f_connections]]
    t_vmin = t_bus["vmin"][[findfirst(isequal(c), t_bus["terminals"]) for c in t_connections]]

    f_vmax = f_bus["vmax"][[findfirst(isequal(c), f_bus["terminals"]) for c in f_connections]]
    t_vmax = t_bus["vmax"][[findfirst(isequal(c), t_bus["terminals"]) for c in t_connections]]

    vmin = max.(fill(0.0, length(f_vmax)), f_vmin, t_vmin)
    vmax = min.(fill(2.0, length(f_vmax)), f_vmax, t_vmax)

    state = var(pm, nw, :switch_inline_ne_state, i)

    for (idx, (fc, tc)) in enumerate(zip(f_connections, t_connections))
        JuMP.@constraint(pm.model, (vr_fr[fc]^2 + vi_fr[fc]^2) - (vr_to[tc]^2 + vi_to[tc]^2) <=  (vmax[idx]^2-vmin[idx]^2) * (1-state))
        JuMP.@constraint(pm.model, (vr_fr[fc]^2 + vi_fr[fc]^2) - (vr_to[tc]^2 + vi_to[tc]^2) >= -(vmax[idx]^2-vmin[idx]^2) * (1-state))
    end
end

@doc raw"""
    constraint_mc_switch_inline_ne_ampacity(pm::AbstractUnbalancedRectangularModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing

ACP current limit constraint on switches

math```
p_{fr}^2 + q_{fr}^2 \leq (vr_{fr}^2 + vi_{fr}^2) i_{max}^2
```
"""
function constraint_mc_switch_inline_ne_ampacity(pm::_PMD.AbstractUnbalancedRectangularModels, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, c_rating::Vector{<:Real})::Nothing
    psw_fr = [var(pm, nw, :psw_inline_ne, f_idx)[c] for c in f_connections]
    qsw_fr = [var(pm, nw, :qsw_inline_ne, f_idx)[c] for c in f_connections]
    vr_fr = [var(pm, nw, :vr, f_idx[2])[c] for c in f_connections]
    vi_fr = [var(pm, nw, :vi, f_idx[2])[c] for c in f_connections]

    con(pm, nw, :mu_cm_switch_inline_ne)[f_idx] = mu_cm_fr = [JuMP.@constraint(pm.model, psw_fr[idx]^2 + qsw_fr[idx]^2 .<= (vr_fr[idx]^2 + vi_fr[idx]^2) * c_rating[idx]^2) for idx in findall(c_rating .< Inf)]

    if _IM.report_duals(pm)
        sol(pm, nw, :switch_inline_ne, f_idx[1])[:mu_cm_fr] = mu_cm_fr
    end

    nothing
end
