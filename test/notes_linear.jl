
#######################################################################################################################################################################
function variable_mc_bus_voltage_indicator(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, relax::Bool=false, report::Bool=true)
    if !relax
        z_voltage = var(pm, nw)[:z_voltage] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :bus)], base_name="$(nw)_z_voltage",
            binary = true,
            start = comp_start_value(ref(pm, nw, :bus, i), "z_voltage_start", 1.0)
        )
    else
        z_voltage =var(pm, nw)[:z_voltage] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :bus)], base_name="$(nw)_z_voltage",
            lower_bound = 0,
            upper_bound = 1,
            start = comp_start_value(ref(pm, nw, :bus, i), "z_voltage_start", 1.0)
        )
    end

    report && _IM.sol_component_value(pm, nw, :bus, :status, ids(pm, nw, :bus), z_voltage)
end

#######################################################################################################################################################################
function variable_mc_bus_voltage_on_off(pm::LPUBFDiagModel; kwargs...) # need LPUBFDiagModel for this one 
    variable_mc_bus_voltage_magnitude_sqr_on_off(pm; kwargs...)
end
function variable_mc_bus_voltage_magnitude_sqr_on_off(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    w = var(pm, nw)[:w] = Dict(i => JuMP.@variable(pm.model,
        [c in 1:ncnds], base_name="$(nw)_w_$(i)",
        start = comp_start_value(ref(pm, nw, :bus, i), "w_start", c, 1.001)
    ) for i in ids(pm, nw, :bus))

    if bounded
        for (i, bus) in ref(pm, nw, :bus)
            set_lower_bound.(w[i], 0.0)

            if haskey(bus, "vmax")
                set_upper_bound.(w[i], bus["vmax"].^2)
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :bus, :w, ids(pm, nw, :bus), w)
end

#######################################################################################################################################################################
function variable_mc_branch_power(pm::LPUBFDiagModel; n_cond::Int=3, nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    @assert(n_cond == 3)
    variable_mc_branch_power_real(pm, nw=nw, bounded=bounded)
    variable_mc_branch_power_imaginary(pm, nw=nw, bounded=bounded)
end

"variable: `p[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_branch_power_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    p = var(pm, nw)[:p] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_p_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "p_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs)
            smax = _calc_branch_power_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i))
            set_upper_bound.(p[(l,i,j)],  smax)
            set_lower_bound.(p[(l,i,j)], -smax)
        end
    end

    for (l,branch) in ref(pm, nw, :branch)
        if haskey(branch, "pf_start")
            f_idx = (l, branch["f_bus"], branch["t_bus"])
            set_start_value(p[f_idx], branch["pf_start"])
        end
        if haskey(branch, "pt_start")
            t_idx = (l, branch["t_bus"], branch["f_bus"])
            set_start_value(p[t_idx], branch["pt_start"])
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :branch, :pf, :pt, ref(pm, nw, :arcs_from), ref(pm, nw, :arcs_to), p)
end

"variable: `q[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_mc_branch_power_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    q = var(pm, nw)[:q] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_q_$((l,i,j))",
            start = comp_start_value(ref(pm, nw, :branch, l), "q_start", c, 0.0)
        ) for (l,i,j) in ref(pm, nw, :arcs)
    )

    if bounded
        for (l,i,j) in ref(pm, nw, :arcs)
            smax = _calc_branch_power_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i))
            set_upper_bound.(q[(l,i,j)],  smax)
            set_lower_bound.(q[(l,i,j)], -smax)
        end
    end

    for (l,branch) in ref(pm, nw, :branch)
        if haskey(branch, "qf_start")
            f_idx = (l, branch["f_bus"], branch["t_bus"])
            set_start_value(q[f_idx], branch["qf_start"])
        end
        if haskey(branch, "qt_start")
            t_idx = (l, branch["t_bus"], branch["f_bus"])
            set_start_value(q[t_idx], branch["qt_start"])
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :branch, :qf, :qt, ref(pm, nw, :arcs_from), ref(pm, nw, :arcs_to), q)
end

#######################################################################################################################################################################
"Creates variables for both `active` and `reactive` power flow at each transformer."
function variable_mc_transformer_power(pm::_PM.AbstractPowerModel; kwargs...)
    variable_mc_transformer_power_real(pm; kwargs...)
    variable_mc_transformer_power_imaginary(pm; kwargs...)
end


"Create variables for the active power flowing into all transformer windings."
function variable_mc_transformer_power_real(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    pt = var(pm, nw)[:pt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_pt_$((l,i,j))",
        ) for (l,i,j) in ref(pm, nw, :arcs_trans)
    )

    if bounded
        for arc in ref(pm, nw, :arcs_from_trans)
            (t,i,j) = arc
            rate_a_fr, rate_a_to = _calc_transformer_power_ub_frto(ref(pm, nw, :transformer, t), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))
            set_lower_bound.(pt[(t,i,j)], -rate_a_fr)
            set_upper_bound.(pt[(t,i,j)],  rate_a_fr)
            set_lower_bound.(pt[(t,j,i)], -rate_a_to)
            set_upper_bound.(pt[(t,j,i)],  rate_a_to)
        end
    end

    for (l,transformer) in ref(pm, nw, :transformer)
        if haskey(transformer, "pf_start")
            f_idx = (l, transformer["f_bus"], transformer["t_bus"])
            set_start_value(pt[f_idx], branch["pf_start"])
        end
        if haskey(transformer, "pt_start")
            t_idx = (l, transformer["t_bus"], transformer["f_bus"])
            set_start_value(pt[t_idx], transformer["pt_start"])
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :transformer, :pf, :pt, ref(pm, nw, :arcs_from_trans), ref(pm, nw, :arcs_to_trans), pt)
end


"Create variables for the reactive power flowing into all transformer windings."
function variable_mc_transformer_power_imaginary(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    qt = var(pm, nw)[:qt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:ncnds], base_name="$(nw)_qt_$((l,i,j))",
            start = 0.0
        ) for (l,i,j) in ref(pm, nw, :arcs_trans)
    )

    if bounded
        for arc in ref(pm, nw, :arcs_from_trans)
            (t,i,j) = arc
            rate_a_fr, rate_a_to = _calc_transformer_power_ub_frto(ref(pm, nw, :transformer, t), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))

            set_lower_bound.(qt[(t,i,j)], -rate_a_fr)
            set_upper_bound.(qt[(t,i,j)],  rate_a_fr)
            set_lower_bound.(qt[(t,j,i)], -rate_a_to)
            set_upper_bound.(qt[(t,j,i)],  rate_a_to)
        end
    end

    for (l,transformer) in ref(pm, nw, :transformer)
        if haskey(transformer, "qf_start")
            f_idx = (l, transformer["f_bus"], transformer["t_bus"])
            set_start_value(qt[f_idx], branch["qf_start"])
        end
        if haskey(transformer, "qt_start")
            t_idx = (l, transformer["t_bus"], transformer["f_bus"])
            set_start_value(qt[t_idx], transformer["qt_start"])
        end
    end

    report && _IM.sol_component_value_edge(pm, nw, :transformer, :qf, :qt, ref(pm, nw, :arcs_from_trans), ref(pm, nw, :arcs_to_trans), qt)
end

#######################################################################################################################################################################
"Create variables for generator status"
function variable_mc_gen_indicator(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, relax::Bool=false, report::Bool=true)
    if !relax
        z_gen = var(pm, nw)[:z_gen] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :gen)], base_name="$(nw)_z_gen",
            binary = true,
            start = comp_start_value(ref(pm, nw, :gen, i), "z_gen_start", 1.0)
        )
    else
        z_gen = var(pm, nw)[:z_gen] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :gen)], base_name="$(nw)_z_gen",
            lower_bound = 0,
            upper_bound = 1,
            start = comp_start_value(ref(pm, nw, :gen, i), "z_gen_start", 1.0)
        )
    end

    report && _IM.sol_component_value(pm, nw, :gen, :gen_status, ids(pm, nw, :gen), z_gen)
end

#######################################################################################################################################################################
function variable_mc_gen_power_setpoint_on_off(pm::_PM.AbstractPowerModel; kwargs...)
    variable_mc_gen_power_setpoint_real_on_off(pm; kwargs...)
    variable_mc_gen_power_setpoint_imaginary_on_off(pm; kwargs...)
end

function variable_mc_gen_power_setpoint_real_on_off(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(conductor_ids(pm, nw))

    pg = var(pm, nw)[:pg] = Dict(i => JuMP.@variable(pm.model,
        [cnd in 1:ncnds], base_name="$(nw)_pg_$(i)",
        start = comp_start_value(ref(pm, nw, :gen, i), "pg_start", cnd, 0.0)
    ) for i in ids(pm, nw, :gen))

    if bounded
        for (i, gen) in ref(pm, nw, :gen)
            if haskey(gen, "pmin")
                set_lower_bound.(pg[i], gen["pmin"])
            end

            if haskey(gen, "pmax")
                set_upper_bound.(pg[i], gen["pmax"])
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :gen, :pg, ids(pm, nw, :gen), pg)
end


function variable_mc_gen_power_setpoint_imaginary_on_off(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    qg = var(pm, nw)[:qg] = Dict(i => JuMP.@variable(pm.model,
        [cnd in 1:ncnds], base_name="$(nw)_qg_$(i)",
        start = comp_start_value(ref(pm, nw, :gen, i), "qg_start", cnd, 0.0)
    ) for i in ids(pm, nw, :gen))

    if bounded
        for (i, gen) in ref(pm, nw, :gen)
            if haskey(gen, "qmin")
                set_lower_bound.(qg[i], gen["qmin"])
            end

            if haskey(gen, "qmax")
                set_upper_bound.(qg[i], gen["qmax"])
            end
        end
    end

    report && _IM.sol_component_value(pm, nw, :gen, :qg, ids(pm, nw, :gen), qg)
end
#######################################################################################################################################################################
"Create variables for demand status"
function variable_mc_load_indicator(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, relax::Bool=false, report::Bool=true)
    # this is not indexedon cnd; why used in start value?
    cnd = 1
    if relax
        z_demand = var(pm, nw)[:z_demand] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :load)], base_name="$(nw)_z_demand",
            lower_bound = 0,
            upper_bound = 1,
            start = comp_start_value(ref(pm, nw, :load, i), "z_demand_on_start", 1.0)
        )
    else
        z_demand = var(pm, nw)[:z_demand] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :load)], base_name="$(nw)_z_demand",
            binary = true,
            start = comp_start_value(ref(pm, nw, :load, i), "z_demand_on_start", 1.0)
        )
    end

    # expressions for pd and qd
    pd = var(pm, nw)[:pd] = Dict(i => var(pm, nw)[:z_demand][i].*ref(pm, nw, :load, i)["pd"]
     for i in ids(pm, nw, :load))
    qd = var(pm, nw)[:qd] = Dict(i => var(pm, nw)[:z_demand][i].*ref(pm, nw, :load, i)["qd"]
     for i in ids(pm, nw, :load))

    report && _IM.sol_component_value(pm, nw, :load, :status, ids(pm, nw, :load), z_demand)
    report && _IM.sol_component_value(pm, nw, :load, :pd, ids(pm, nw, :load), pd)
    report && _IM.sol_component_value(pm, nw, :load, :qd, ids(pm, nw, :load), qd)
end

#######################################################################################################################################################################
function variable_mc_shunt_indicator(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, relax=false, report::Bool=true)
    # this is not indexedon cnd; why used in start value?
    cnd = 1
    if relax
        z_shunt = var(pm, nw)[:z_shunt] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :shunt)], base_name="$(nw)_z_shunt",
            lower_bound = 0,
            upper_bound = 1,
            start = comp_start_value(ref(pm, nw, :shunt, i), "z_shunt_on_start", 1.0)
        )
    else
        z_shunt = var(pm, nw)[:z_shunt] = JuMP.@variable(pm.model,
            [i in ids(pm, nw, :shunt)], base_name="$(nw)_z_shunt",
            binary=true,
            start = comp_start_value(ref(pm, nw, :shunt, i), "z_shunt_on_start", 1.0)
        )
    end

    report && _IM.sol_component_value(pm, nw, :shunt, :status, ids(pm, nw, :shunt), z_shunt)
end

#######################################################################################################################################################################
""
function constraint_mc_model_voltage(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw)
    constraint_mc_model_voltage(pm, nw)
end

function constraint_mc_model_voltage(pm::_PM.AbstractPowerModel, n::Int)
end

#######################################################################################################################################################################
function constraint_mc_theta_ref(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    va_ref = ref(pm, nw, :bus, i, "va")
    constraint_mc_theta_ref(pm, nw, i, va_ref)
end

function constraint_mc_theta_ref(pm::LPUBFDiagModel, n::Int, i::Int, va_ref)
    ncnds = length(conductor_ids(pm))
    @assert(ncnds >= 2)
    w = var(pm, n, :w)[i]
    JuMP.@constraint(pm.model, w[2:ncnds]   .== w[1])
end

#######################################################################################################################################################################
function constraint_mc_bus_voltage_on_off(pm::_PM.AbstractPowerModel; nw::Int=pm.cnw, kwargs...)
    constraint_mc_bus_voltage_on_off(pm, nw; kwargs...)
end

"on/off bus voltage constraint for relaxed forms"
function constraint_mc_bus_voltage_on_off(pm::_PM.AbstractWModels, n::Int; kwargs...)
    for (i, bus) in ref(pm, n, :bus)
        constraint_mc_bus_voltage_magnitude_sqr_on_off(pm, i, nw=n)
    end
end

function constraint_mc_bus_voltage_magnitude_sqr_on_off(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    bus = ref(pm, nw, :bus, i)

    constraint_mc_bus_voltage_magnitude_sqr_on_off(pm, nw, i, bus["vmin"], bus["vmax"])
end

"on/off bus voltage magnitude squared constraint for relaxed formulations"
function constraint_mc_bus_voltage_magnitude_sqr_on_off(pm::_PM.AbstractPowerModel, n::Int, i::Int, vmin, vmax)
    w = var(pm, n, :w, i)
    z_voltage = var(pm, n, :z_voltage, i)

    for c in conductor_ids(pm, n)
        if isfinite(vmax[c])
            JuMP.@constraint(pm.model, w[c] <= vmax[c]^2*z_voltage)
        end

        if isfinite(vmin[c])
            JuMP.@constraint(pm.model, w[c] >= vmin[c]^2*z_voltage)
        end
    end
end

#######################################################################################################################################################################
function constraint_mc_gen_power_on_off(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    gen = ref(pm, nw, :gen, i)
    ncnds = length(conductor_ids(pm; nw=nw))

    pmin = get(gen, "pmin", fill(-Inf, ncnds))
    pmax = get(gen, "pmax", fill( Inf, ncnds))
    qmin = get(gen, "qmin", fill(-Inf, ncnds))
    qmax = get(gen, "qmax", fill( Inf, ncnds))

    constraint_mc_gen_power_on_off(pm, nw, i, pmin, pmax, qmin, qmax)
end

function constraint_mc_gen_power_on_off(pm::_PM.AbstractPowerModel, n::Int, i::Int, pmin, pmax, qmin, qmax)
    pg = var(pm, n, :pg, i)
    qg = var(pm, n, :qg, i)
    z = var(pm, n, :z_gen, i)

    for c in conductor_ids(pm, n)
        if isfinite(pmax[c])
            JuMP.@constraint(pm.model, pg[c] .<= pmax[c].*z)
        end

        if isfinite(pmin[c])
            JuMP.@constraint(pm.model, pg[c] .>= pmin[c].*z)
        end

        if isfinite(qmax[c])
            JuMP.@constraint(pm.model, qg[c] .<= qmax[c].*z)
        end

        if isfinite(qmin[c])
            JuMP.@constraint(pm.model, qg[c] .>= qmin[c].*z)
        end
    end
end

#######################################################################################################################################################################
"KCL for load shed problem with transformers"
function constraint_mc_shed_power_balance(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs, i)
    bus_arcs_sw = ref(pm, nw, :bus_arcs_sw, i)
    bus_arcs_trans = ref(pm, nw, :bus_arcs_trans, i)
    bus_gens = ref(pm, nw, :bus_gens, i)
    bus_storage = ref(pm, nw, :bus_storage, i)
    bus_loads = ref(pm, nw, :bus_loads, i)
    bus_shunts = ref(pm, nw, :bus_shunts, i)

    bus_pd = Dict(k => ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs = Dict(k => ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    if !haskey(con(pm, nw), :lam_kcl_r)
        con(pm, nw)[:lam_kcl_r] = Dict{Int,Array{JuMP.ConstraintRef}}()
    end

    if !haskey(con(pm, nw), :lam_kcl_i)
        con(pm, nw)[:lam_kcl_i] = Dict{Int,Array{JuMP.ConstraintRef}}()
    end

    constraint_mc_shed_power_balance(pm, nw, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
end

"KCL for load shed problem with transformers (AbstractWForms)"
function constraint_mc_shed_power_balance(pm::_PM.AbstractWModels, nw::Int, i, bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    w        = var(pm, nw, :w, i)
    p        = get(var(pm, nw),    :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    q        = get(var(pm, nw),    :q, Dict()); _PM._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg       = get(var(pm, nw),   :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    qg       = get(var(pm, nw),   :qg, Dict()); _PM._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps       = get(var(pm, nw),   :ps, Dict()); _PM._check_var_keys(ps, bus_storage, "active power", "storage")
    qs       = get(var(pm, nw),   :qs, Dict()); _PM._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw      = get(var(pm, nw),  :psw, Dict()); _PM._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw      = get(var(pm, nw),  :qsw, Dict()); _PM._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    pt       = get(var(pm, nw),   :pt, Dict()); _PM._check_var_keys(pt, bus_arcs_trans, "active power", "transformer")
    qt       = get(var(pm, nw),   :qt, Dict()); _PM._check_var_keys(qt, bus_arcs_trans, "reactive power", "transformer")
    z_demand = var(pm, nw, :z_demand)
    z_shunt  = var(pm, nw, :z_shunt)

    cnds = conductor_ids(pm; nw=nw)
    ncnds = length(cnds)

    Gt = isempty(bus_gs) ? fill(0.0, ncnds, ncnds) : sum(values(bus_gs))
    Bt = isempty(bus_bs) ? fill(0.0, ncnds, ncnds) : sum(values(bus_bs))

    bus_GsBs = [(n,bus_gs[n], bus_bs[n]) for n in keys(bus_gs)]

    cstr_p = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        + sum(pt[a_trans] for a_trans in bus_arcs_trans)
        .==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd .*z_demand[n] for (n,pd) in bus_pd)
        - sum(z_shunt[n].*(w.*diag(Gt')) for (n,Gs,Bs) in bus_GsBs)
    )
    cstr_q = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        + sum(qt[a_trans] for a_trans in bus_arcs_trans)
        .==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd.*z_demand[n] for (n,qd) in bus_qd)
        - sum(z_shunt[n].*(-w.*diag(Bt')) for (n,Gs,Bs) in bus_GsBs)
    )

    con(pm, nw, :lam_kcl_r)[i] = isa(cstr_p, Array) ? cstr_p : [cstr_p]
    con(pm, nw, :lam_kcl_i)[i] = isa(cstr_q, Array) ? cstr_q : [cstr_q]

    if _IM.report_duals(pm)
        sol(pm, nw, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, nw, :bus, i)[:lam_kcl_i] = cstr_q
    end
end

#######################################################################################################################################################################
"ohms constraint for branches on the from-side"
function constraint_mc_ohms_yt_from(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PM.calc_branch_y(branch)
    tr, ti = _PM.calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]

    constraint_mc_ohms_yt_from(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end

"delegate back to PowerModels"
function constraint_mc_ohms_yt_from(pm::_PM.AbstractWModels, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    _PM.constraint_ohms_yt_from(pm, n, c, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
end

#######################################################################################################################################################################
"ohms constraint for branches on the to-side"
function constraint_mc_ohms_yt_to(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PM.calc_branch_y(branch)
    tr, ti = _PM.calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    constraint_mc_ohms_yt_to(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end

function constraint_mc_ohms_yt_to(pm::_PM.AbstractWModels, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    _PM.constraint_ohms_yt_to(pm, n, c, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
end

#######################################################################################################################################################################
"This is duplicated at PowerModelsDistribution level to correctly handle the indexing of the shunts."
function constraint_mc_voltage_angle_difference(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    pair = (f_bus, t_bus)

    constraint_mc_voltage_angle_difference(pm, nw, f_idx, branch["angmin"], branch["angmax"])
end

function constraint_mc_voltage_angle_difference(pm::_PM.AbstractWModels, n::Int, f_idx, angmin, angmax)
    i, f_bus, t_bus = f_idx

    ncnds = length(conductor_ids(pm, n))

    w_fr = var(pm, n, :w, f_bus)
    w_to = var(pm, n, :w, t_bus)
    wr   = [var(pm, n, :wr)[(f_bus, t_bus, c, c)] for c in 1:ncnds]
    wi   = [var(pm, n, :wi)[(f_bus, t_bus, c, c)] for c in 1:ncnds]

    JuMP.@constraint(pm.model, wi .<= tan.(angmax).*wr)
    JuMP.@constraint(pm.model, wi .>= tan.(angmin).*wr)

    for c in 1:ncnds
        _PM.cut_complex_product_and_angle_difference(pm.model, w_fr[c], w_to[c], wr[c], wi[c], angmin[c], angmax[c])
    end
end

#######################################################################################################################################################################
"branch thermal constraints from"
function constraint_mc_thermal_limit_from(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    if haskey(branch, "rate_a")
        constraint_mc_thermal_limit_from(pm, nw, f_idx, branch["rate_a"])
    end
end

"Generic thermal limit constraint from-side"
function constraint_mc_thermal_limit_from(pm::_PM.AbstractPowerModel, n::Int, f_idx, rate_a)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)

    mu_sm_fr = JuMP.@constraint(pm.model, p_fr.^2 + q_fr.^2 .<= rate_a.^2)

    if _IM.report_duals(pm)
        sol(pm, n, :branch, f_idx[1])[:mu_sm_fr] = mu_sm_fr
    end
end

#######################################################################################################################################################################
"branch thermal constraints to"
function constraint_mc_thermal_limit_to(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)

    if haskey(branch, "rate_a")
        constraint_mc_thermal_limit_to(pm, nw, t_idx, branch["rate_a"])
    end
end

"Generic thermal limit constraint to-side"
function constraint_mc_thermal_limit_to(pm::_PM.AbstractPowerModel, n::Int, t_idx, rate_a)
    p_to = var(pm, n, :p, t_idx)
    q_to = var(pm, n, :q, t_idx)

    mu_sm_to = JuMP.@constraint(pm.model, p_to.^2 + q_to.^2 .<= rate_a.^2)

    if _IM.report_duals(pm)
        sol(pm, n, :branch, t_idx[1])[:mu_sm_to] = mu_sm_to
    end
end

#######################################################################################################################################################################

#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
#######################################################################################################################################################################