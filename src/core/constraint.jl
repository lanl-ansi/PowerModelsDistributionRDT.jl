
function constraint_ue(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, base_nw::Int, gen::Int)
    ue = _PMD.var(pm, base_nw, :ue, gen)
    ue_s = _PMD.var(pm, nw, :ue_s, gen)

    JuMP.@constraint(pm.model, ue == ue_s)
end

function constraint_xe(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, base_nw::Int, branch::Int)
    xe = _PMD.var(pm, base_nw, :xe, branch)
    xe_s = _PMD.var(pm, nw, :xe_s, branch)

    JuMP.@constraint(pm.model, xe == xe_s)
end

"""
    Adding an inline switch to existing line is modeled by adding the switch explcitly as a seperate sequential edge that is a lossless closed edge when the switch is not built
    This constraint forces the open state when this happens
"""
function constraint_te(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, base_nw::Int, switch::Int)
    te = _PMD.var(pm, base_nw, :te, switch)
    te_s = _PMD.var(pm, nw, :switch_inline_ne_state, switch)

    JuMP.@constraint(pm.model, te >= 1-te_s)
end

function constraint_he(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, base_nw::Int, branch::Int)
    he = _PMD.var(pm, base_nw, :he, branch)
    he_s = _PMD.var(pm, nw, :he_s, branch)

    JuMP.@constraint(pm.model, he == he_s)
end


"""
    constraint_mc_thermal_limit_from_damaged(pm::AbstractUnbalancedPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})::Nothing

Generic thermal limit constraint for branches (from-side) that are damaged
"""
function constraint_mc_thermal_limit_from_damaged(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})::Nothing
    p_fr = [_PMD.var(pm, nw, :p, f_idx)[c] for c in f_connections]
    q_fr = [_PMD.var(pm, nw, :q, f_idx)[c] for c in f_connections]
    he_s = _PMD.var(pm, nw, :he_s, f_idx[1])

    _PMD.con(pm, nw, :mu_sm_branch)[f_idx] = mu_sm_fr = [JuMP.@constraint(pm.model, p_fr[idx]^2 + q_fr[idx]^2 <= rate_a[idx]^2 * he_s) for idx in f_connections]

    if _IM.report_duals(pm)
        _PMD.sol(pm, nw, :branch, f_idx[1])[:mu_sm_fr] = mu_sm_fr
    end
    nothing
end


"""
    constraint_mc_thermal_limit_to_damaged(pm::AbstractUnbalancedPowerModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})::Nothing

Generic thermal limit constraint for branches (to-side) that are damaged
"""
function constraint_mc_thermal_limit_to_damaged(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})::Nothing
    p_to = [_PMD.var(pm, nw, :p, t_idx)[c] for c in t_connections]
    q_to = [_PMD.var(pm, nw, :q, t_idx)[c] for c in t_connections]
    he_s = _PMD.var(pm, nw, :he_s, t_idx[1])

    _PMD.con(pm, nw, :mu_sm_branch)[t_idx] = mu_sm_to = [JuMP.@constraint(pm.model, p_to[idx]^2 + q_to[idx]^2 <= rate_a[idx]^2 * he_s) for idx in t_connections]

    if _IM.report_duals(pm)
        _PMD.sol(pm, nw, :branch, t_idx[1])[:mu_sm_to] = mu_sm_to
    end
    nothing
end


"""
    constraint_mc_thermal_limit_from_ne(pm::AbstractUnbalancedPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})::Nothing

Generic thermal limit constraint for branches (from-side) that are expansions
"""
function constraint_mc_thermal_limit_from_ne(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})::Nothing
    p_fr = [_PMD.var(pm, nw, :p_ne, f_idx)[c] for c in f_connections]
    q_fr = [_PMD.var(pm, nw, :q_ne, f_idx)[c] for c in f_connections]
    xe_s = _PMD.var(pm, nw, :xe_s, f_idx[1])

    _PMD.con(pm, nw, :mu_sm_branch_ne)[f_idx] = mu_sm_fr = [JuMP.@constraint(pm.model, p_fr[idx]^2 + q_fr[idx]^2 <= rate_a[idx]^2 * xe_s) for idx in f_connections]

    if _IM.report_duals(pm)
        _PMD.sol(pm, nw, :branch_ne, f_idx[1])[:mu_sm_fr_ne] = mu_sm_fr
    end
    nothing
end


"""
    constraint_mc_thermal_limit_to_damaged(pm::AbstractUnbalancedPowerModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})::Nothing

Generic thermal limit constraint for branches (to-side) that are expansions
"""
function constraint_mc_thermal_limit_to_ne(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})::Nothing
    p_to = [_PMD.var(pm, nw, :p_ne, t_idx)[c] for c in t_connections]
    q_to = [_PMD.var(pm, nw, :q_ne, t_idx)[c] for c in t_connections]
    xe_s = _PMD.var(pm, nw, :xe_s, t_idx[1])

    _PMD.con(pm, nw, :mu_sm_branch_ne)[t_idx] = mu_sm_to = [JuMP.@constraint(pm.model, p_to[idx]^2 + q_to[idx]^2 <= rate_a[idx]^2 * xe_s) for idx in t_connections]

    if _IM.report_duals(pm)
        _PMD.sol(pm, nw, :branch_ne, t_idx[1])[:mu_sm_to_ne] = mu_sm_to
    end
    nothing
end


function constraint_critical_load(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, loads::Set{Int64}, limit::Float64, total_pd::Vector{Float64}, total_qd::Vector{Float64})
    pd = _PMD.var(pm, nw, :pd)
    qd = _PMD.var(pm, nw, :qd)

    all_phase_pd = sum(total_pd)
    all_phase_qd = sum(total_qd)

    JuMP.@constraint(pm.model, limit * all_phase_pd <= sum(sum(pd[i]) for i in loads))
    JuMP.@constraint(pm.model, limit * all_phase_qd <= sum(sum(qd[i]) for i in loads))
end

function constraint_total_load(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, loads::Set{Int64}, limit::Float64, total_pd::Vector{Float64}, total_qd::Vector{Float64})
    pd = _PMD.var(pm, nw, :pd)
    qd = _PMD.var(pm, nw, :qd)

    all_phase_pd = sum(total_pd)
    all_phase_qd = sum(total_qd)

    JuMP.@constraint(pm.model, limit * all_phase_pd <= sum(sum(pd[i]) for i in loads))
    JuMP.@constraint(pm.model, limit * all_phase_qd <= sum(sum(qd[i]) for i in loads))
end



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


"on/off constraint for ne generators"
function constraint_mc_gen_power_ne(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, i::Int, connections::Vector{<:Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real})
    pg = _PMD.var(pm, nw, :pg_ne, i)
    qg = _PMD.var(pm, nw, :qg_ne, i)
    z = _PMD.var(pm, nw, :z_gen_ne, i)
    u = _PMD.var(pm, nw, :ue_s, i)

    mc = JuMP.@variable(pm.model,
            lower_bound = 0,
            upper_bound = 1
         )

    _IM.relaxation_product(pm.model, z, u, mc)

    for (idx, c) in enumerate(connections)
        if isfinite(pmax[idx])
            JuMP.@constraint(pm.model, pg[c] .<= pmax[idx].*mc)
        end

        if isfinite(pmin[idx])
            JuMP.@constraint(pm.model, pg[c] .>= pmin[idx].*mc)
        end

        if isfinite(qmax[idx])
            JuMP.@constraint(pm.model, qg[c] .<= qmax[idx].*mc)
        end

        if isfinite(qmin[idx])
            JuMP.@constraint(pm.model, qg[c] .>= qmin[idx].*mc)
        end
    end
    nothing
end


@doc raw"""
    constraint_radial_topology_ne(pm::AbstractUnbalancedPowerModel, nw::Int; relax::Bool=false)

Constraint to enforce a radial topology

See 10.1109/TSG.2020.2985087

```math
\begin{align}
\mathbf{\beta} \in \mathbf{\Omega} \\
\alpha_{ij} \leq \beta_{ij},\forall(i,j) \in L \\
\sum_{\substack{(j,i_r)\in L}}f^{k}_{ji_r} - \sum_{\substack{(i_r,j)\in L}}f^{k}_{i_rj}=-1,~\forall k \in N\setminus i_r \\
\sum_{\substack{(j,k)\in L}}f^{k}_{jk} - \sum_{\substack{(k,j)\in L}}f^k_{kj} = 1,~\forall k \in N\setminus i_r \\
\sum_{\substack{(j,i)\in L}}f^k_{ji}-\sum_{\substack{(i,j)\in L}}f^k_{ij}=0,~\forall k \in N\setminus i_r,\forall i \in N\setminus {i_r,k} \\
0 \leq f^k_{ij} \leq \lambda_{ij},0 \leq f^k_{ji} \leq \lambda_{ji},\forall k \in N\setminus i_r,\forall(i,j)\in L \\
\sum_{\substack{(i,j)\in L}}\left(\lambda_{ij} + \lambda_{ji} \right ) = \left | N \right | - 1 \\
\lambda_{ij} + \lambda_{ji} = \beta_{ij},\forall(i,j)\in L \\
\lambda_{ij},\lambda_{ji}\in\left \{ 0,1 \right \},\forall(i,j)\in L
\end{align}
```
"""
function constraint_radial_topology_ne(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int; relax::Bool=false)
    # doi: 10.1109/TSG.2020.2985087
    _PMD.var(pm, nw)[:f] = Dict{Tuple{Int,Int,Int},JuMP.VariableRef}()
    _PMD.var(pm, nw)[:lambda] = Dict{Tuple{Int,Int},JuMP.VariableRef}()
    _PMD.var(pm, nw)[:beta] = Dict{Tuple{Int,Int},JuMP.VariableRef}()
    _PMD.var(pm, nw)[:alpha] = Dict{Tuple{Int,Int},Union{JuMP.VariableRef,JuMP.AffExpr,Int}}()

    # "real" node and branch sets, based on load blocks
    N₀ = _PMD.ids(pm, nw, :blocks)
    L₀ = _PMD.ref(pm, nw, :block_pairs)

    # Add "virtual" iᵣ to N
    virtual_iᵣ = maximum(N₀)+1
    N = [N₀..., virtual_iᵣ]
    iᵣ = [virtual_iᵣ]

    # create a set L of all branches, including virtual branches between iᵣ and all other nodes in L₀
    L = [L₀..., [(virtual_iᵣ, n) for n in N₀]...]

    # create a set L′ that inlcudes the branch reverses
#    L′ = union(L, Set([(j,i) for (i,j) in L]))

    # create variables fᵏ and λ over all L, including virtual branches connected to iᵣ
#    for (i,j) in L′
#        for k in filter(kk->kk∉iᵣ,N)
#            var(pm, nw, :f)[(k, i, j)] = JuMP.@variable(pm.model, base_name="$(nw)_f_$((k,i,j))", start=(k,i,j) == (k,virtual_iᵣ,k) ? 1 : 0)
#        end
#        var(pm, nw, :lambda)[(i,j)] = JuMP.@variable(pm.model, base_name="$(nw)_lambda_$((i,j))", binary=!relax, lower_bound=0, upper_bound=1, start=(i,j) == (virtual_iᵣ,j) ? 1 : 0)

        # create variable β over only original set L₀
#        if (i,j) ∈ L₀
#            var(pm, nw, :beta)[(i,j)] = JuMP.@variable(pm.model, base_name="$(nw)_beta_$((i,j))", lower_bound=0, upper_bound=1)
#        end
#    end

    # create an aux varible α that maps to the switch states
#    switch_lookup = Dict{Tuple{Int,Int},Vector{Int}}((ref(pm, nw, :bus_block_map, sw["f_bus"]), ref(pm, nw, :bus_block_map, sw["t_bus"])) => Int[ss for (ss,ssw) in ref(pm, nw, :switch) if (ref(pm, nw, :bus_block_map, sw["f_bus"])==ref(pm, nw, :bus_block_map, ssw["f_bus"]) && ref(pm, nw, :bus_block_map, sw["t_bus"])==ref(pm, nw, :bus_block_map, ssw["t_bus"])) || (ref(pm, nw, :bus_block_map, sw["f_bus"])==ref(pm, nw, :bus_block_map, ssw["t_bus"]) && ref(pm, nw, :bus_block_map, sw["t_bus"])==ref(pm, nw, :bus_block_map, ssw["f_bus"]))] for (s,sw) in ref(pm, nw, :switch))
#    for ((i,j), switches) in switch_lookup
#        var(pm, nw, :alpha)[(i,j)] = JuMP.@expression(pm.model, sum(var(pm, nw, :switch_state, s) for s in switches))
#        JuMP.@constraint(pm.model, var(pm, nw, :alpha, (i,j)) <= 1)
#    end

#    f = var(pm, nw, :f)
#    λ = var(pm, nw, :lambda)
#    β = var(pm, nw, :beta)
#    α = var(pm, nw, :alpha)

    # Eq. (1) -> Eqs. (3-8)
#    for k in filter(kk->kk∉iᵣ,N)
#        # Eq. (3)
#        for _iᵣ in iᵣ
#            jiᵣ = filter(((j,i),)->i==_iᵣ&&i!=j,L)
#            iᵣj = filter(((i,j),)->i==_iᵣ&&i!=j,L)
#            if !(isempty(jiᵣ) && isempty(iᵣj))
#                c = JuMP.@constraint(
#                    pm.model,
#                    sum(f[(k,j,i)] for (j,i) in jiᵣ) -
#                    sum(f[(k,i,j)] for (i,j) in iᵣj)
#                    ==
#                    -1.0
#                )
#            end
#        end

        # Eq. (4)
#        jk = filter(((j,i),)->i==k&&i!=j,L′)
#        kj = filter(((i,j),)->i==k&&i!=j,L′)
#        if !(isempty(jk) && isempty(kj))
#            c = JuMP.@constraint(
#                pm.model,
#                sum(f[(k,j,k)] for (j,i) in jk) -
#                sum(f[(k,k,j)] for (i,j) in kj)
#                ==
#                1.0
#            )
#        end

        # Eq. (5)
#        for i in filter(kk->kk∉iᵣ&&kk!=k,N)
#            ji = filter(((j,ii),)->ii==i&&ii!=j,L′)
#            ij = filter(((ii,j),)->ii==i&&ii!=j,L′)
#            if !(isempty(ji) && isempty(ij))
#                c = JuMP.@constraint(
#                    pm.model,
#                    sum(f[(k,j,i)] for (j,ii) in ji) -
#                    sum(f[(k,i,j)] for (ii,j) in ij)
#                    ==
#                    0.0
#                )
#            end
#        end

        # Eq. (6)
#        for (i,j) in L
#            JuMP.@constraint(pm.model, f[(k,i,j)] >= 0)
#            JuMP.@constraint(pm.model, f[(k,i,j)] <= λ[(i,j)])
#            JuMP.@constraint(pm.model, f[(k,j,i)] >= 0)
#            JuMP.@constraint(pm.model, f[(k,j,i)] <= λ[(j,i)])
#        end
#    end

    # Eq. (7)
#    JuMP.@constraint(pm.model, sum((λ[(i,j)] + λ[(j,i)]) for (i,j) in L) == length(N) - 1)

    # Connect λ and β, map β back to α, over only real switches (L₀)
#    for (i,j) in L₀
        # Eq. (8)
#        JuMP.@constraint(pm.model, λ[(i,j)] + λ[(j,i)] == β[(i,j)])

        # Eq. (2)
#        JuMP.@constraint(pm.model, α[(i,j)] <= β[(i,j)])
#    end
end


@doc raw"""
    constraint_mc_switch_power_open_close(
        pm::AbstractUnbalancedPowerModel,
        nw::Int,
        i::Int,
        f_bus::Int,
        t_bus::Int,
        f_connections::Vector{Int},
        t_connections::Vector{Int}
    )

generic switch power open/closed constraint

```math
\begin{align}
& S^{sw}_{i,c} \leq S^{swu}_{i,c} z^{sw}_i\ \forall i \in S,\forall c \in C \\
& S^{sw}_{i,c} \geq -S^{swu}_{i,c} z^{sw}_i\ \forall i \in S,\forall c \in C
\end{align}
```
"""
function constraint_mc_switch_inline_ne_power_open_close(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_connections::Vector{Int}, t_connections::Vector{Int})
    psw = _PMD.var(pm, nw, :psw_inline_ne, (i, f_bus, t_bus))
    qsw = _PMD.var(pm, nw, :qsw_inline_ne, (i, f_bus, t_bus))

    state = _PMD.var(pm, nw, :switch_inline_ne_state, i)

    rating = min.(fill(1.0, length(f_connections)), _PMD._calc_branch_power_max_frto(_PMD.ref(pm, nw, :switch_inline_ne, i), _PMD.ref(pm, nw, :bus, f_bus), _PMD.ref(pm, nw, :bus, t_bus))...)

    for (idx, c) in enumerate(f_connections)
        JuMP.@constraint(pm.model, psw[c] <=  rating[idx] * state)
        JuMP.@constraint(pm.model, psw[c] >= -rating[idx] * state)
        JuMP.@constraint(pm.model, qsw[c] <=  rating[idx] * state)
        JuMP.@constraint(pm.model, qsw[c] >= -rating[idx] * state)

        # Indicator constraint version, for reference
        # JuMP.@constraint(pm.model, !state => {psw[c] == 0.0})
        # JuMP.@constraint(pm.model, !state => {qsw[c] == 0.0})
    end
end

"""
    constraint_mc_switch_inline_ne_ampacity(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for switch current limit constraint from-side
"""
function constraint_mc_switch_inline_ne_ampacity(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    switch = _PMD.ref(pm, nw, :switch_inline_ne, i)
    f_idx = (i, switch["f_bus"], switch["t_bus"])

    if !haskey(_PMD.con(pm, nw), :mu_cm_switch_inline_ne)
        _PMD.con(pm, nw)[:mu_cm_switch_inline_ne] = Dict{Tuple{Int,Int,Int}, Vector{JuMP.ConstraintRef}}()
    end

    if haskey(switch, "current_rating") && any(switch["current_rating"] .< Inf)
        constraint_mc_switch_inline_ne_ampacity(pm, nw, f_idx, switch["f_connections"], switch["current_rating"])
    end
    nothing
end
