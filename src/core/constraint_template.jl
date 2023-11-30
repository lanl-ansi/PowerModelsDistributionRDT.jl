

function constraint_ue(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=_PMD.nw_id_default, base_nw::Int=_PMD.nw_id_default)
    constraint_ue(pm, nw, base_nw, i)
end

function constraint_xe(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=_PMD.nw_id_default, base_nw::Int=_PMD.nw_id_default)
    constraint_xe(pm, nw, base_nw, i)
end

function constraint_te(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=_PMD.nw_id_default, base_nw::Int=_PMD.nw_id_default)
    constraint_te(pm, nw, base_nw, i)
end

function constraint_he(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=_PMD.nw_id_default, base_nw::Int=_PMD.nw_id_default)
    constraint_he(pm, nw, base_nw, i)
end


#function constraint_branch_be(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
#    constraint_branch_be_p(pm, i, nw)
#    constraint_branch_be_q(pm, i, nw)
#end


#function constraint_branch_be_p(pm::_PMs.AbstractPowerModel, i::Int, nw::Int) # check trans arcs and branch arches
#    arcs = _PMs.ref(pm, nw, :arcs_bal)[i]
#    arcs in _PMs.ref(pm, nw, :arcs_trans) ? trans = true : trans = false
#    constraint_branch_be_p(pm, nw, arcs, trans)
#end

#function constraint_branch_be_q(pm::_PMs.AbstractPowerModel, i::Int, nw::Int)
#    arcs = _PMs.ref(pm, nw, :arcs_bal)[i]
#    arcs in _PMs.ref(pm, nw, :arcs_trans) ? trans = true : trans = false
#    constraint_branch_be_q(pm, nw, arcs, trans)
#end

#function constraint_balance_flow(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
#    constraint_balance_p_flow(pm, i, nw)
#    constraint_balance_q_flow(pm, i, nw)
#end

#function constraint_balance_p_flow(pm::_PMs.AbstractPowerModel, i::Int, nw::Int)
#    arcs = _PMs.ref(pm, nw, :arcs_bal)[i]
#    arcs in _PMs.ref(pm, nw, :arcs_trans) ? trans = true : trans = false
#    constraint_balance_p_flow(pm, nw, arcs, trans)
#end

#function constraint_balance_q_flow(pm::_PMs.AbstractPowerModel, i::Int, nw::Int)
#    arcs = _PMs.ref(pm, nw, :arcs_bal)[i]
#    arcs in _PMs.ref(pm, nw, :arcs_trans) ? trans = true : trans = false
#    constraint_balance_q_flow(pm, nw, arcs, trans)
#end


function constraint_critical_load(pm::_PMD.AbstractUnbalancedPowerModel; nw::Int=_PMD.nw_id_default)
    limit = pm.ref[:it][_PMD.pmd_it_sym][:critical_load_met]
    conductors = _PMD.ref(pm,nw,:conductors)

    total_pd = Vector{Float64}()
    total_qd = Vector{Float64}()

    for i in 1:conductors
        push!(total_pd, 0.0)
        push!(total_qd, 0.0)
    end

    loads = Set{Int64}()
    for (i,load) in _PMD.ref(pm, nw, :load)
        if (load["is_critical"] == true)
            push!(loads,i)
            total_pd = total_pd + load["pd_nom"]
            total_qd = total_qd + load["qd_nom"]
        end
    end
    constraint_critical_load(pm, nw, loads, limit, total_pd, total_qd)
end


function constraint_total_load(pm::_PMD.AbstractUnbalancedPowerModel; nw::Int=_PMD.nw_id_default)
    limit = pm.ref[:it][_PMD.pmd_it_sym][:total_load_met]
    conductors = _PMD.ref(pm,nw,:conductors)

    total_pd = Vector{Float64}()
    total_qd = Vector{Float64}()

    for i in 1:conductors
        push!(total_pd, 0.0)
        push!(total_qd, 0.0)
    end

    loads = Set{Int64}()
    for (i,load) in _PMD.ref(pm, nw, :load)
        push!(loads,i)
        total_pd = total_pd + load["pd_nom"]
        total_qd = total_qd + load["qd_nom"]
    end
    constraint_total_load(pm, nw, loads, limit, total_pd, total_qd)
end


""""
    constraint_mc_power_balance_shed(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
Template function for KCL constraints for load shed problem. It is based on PowerModelsDistribution.constraint_mc_power_balance_shed and
adds balance for inline switches, branch_ne, and switch_ne
"""
function constraint_mc_power_balance_shed_ne(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    bus = _PMD.ref(pm, nw, :bus, i)
    bus_arcs = _PMD.ref(pm, nw, :bus_arcs_conns_branch, i)
    bus_arcs_ne = _PMD.ref(pm, nw, :bus_arcs_conns_branch_ne, i)
    bus_arcs_sw = _PMD.ref(pm, nw, :bus_arcs_conns_switch, i)
    bus_arcs_sw_ne = _PMD.ref(pm, nw, :bus_arcs_conns_switch_inline_ne, i)
    bus_arcs_trans = _PMD.ref(pm, nw, :bus_arcs_conns_transformer, i)
    bus_arcs_trans_ne = _PMD.ref(pm, nw, :bus_arcs_conns_transformer_ne, i)
    bus_gens = _PMD.ref(pm, nw, :bus_conns_gen, i)
    bus_gens_ne = _PMD.ref(pm, nw, :bus_conns_gen_ne, i)
    bus_storage = _PMD.ref(pm, nw, :bus_conns_storage, i)
    bus_loads = _PMD.ref(pm, nw, :bus_conns_load, i)
    bus_shunts = _PMD.ref(pm, nw, :bus_conns_shunt, i)

    if !haskey(_PMD.con(pm, nw), :lam_kcl_r)
        _PMD.con(pm, nw)[:lam_kcl_r] = Dict{Int,Array{JuMP.ConstraintRef}}()
    end

    if !haskey(_PMD.con(pm, nw), :lam_kcl_i)
        _PMD.con(pm, nw)[:lam_kcl_i] = Dict{Int,Array{JuMP.ConstraintRef}}()
    end

    constraint_mc_power_balance_shed_ne(pm, nw, i, bus["terminals"], bus["grounded"], bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_shunts, bus_arcs_ne, bus_arcs_sw_ne, bus_arcs_trans_ne, bus_gens_ne)
    nothing
end


"""
    constraint_mc_thermal_limit_from_damaged(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for branch thermal constraints (from-side) for damaged lines
"""
function constraint_mc_thermal_limit_from_damaged(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = _PMD.ref(pm, nw, :branch, i)
    f_idx = (i, branch["f_bus"], branch["t_bus"])
    bus = _PMD.ref(pm,nw,:bus,branch["f_bus"])

    # using the same constraint store since a branch is either damaged or not damaged
    if !haskey(_PMD.con(pm, nw), :mu_sm_branch)
        _PMD.con(pm, nw)[:mu_sm_branch] = Dict{Tuple{Int,Int,Int}, Vector{JuMP.ConstraintRef}}()
    end

    rating = calc_branch_power_max(branch, bus, _PMD.ref(pm, nw, :total_real_load), _PMD.ref(pm, nw, :total_reactive_load))
    constraint_mc_thermal_limit_from_damaged(pm, nw, f_idx, branch["f_connections"], rating)
    nothing
end


"""
    constraint_mc_thermal_limit_to(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for branch thermal constraints (to-side) for damaged lines
"""
function constraint_mc_thermal_limit_to_damaged(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = _PMD.ref(pm, nw, :branch, i)
    t_idx = (i, branch["t_bus"], branch["f_bus"])
    bus = _PMD.ref(pm,nw,:bus,branch["t_bus"])

    # using the same constraint store since a branch is either damaged or not damaged
    if !haskey(_PMD.con(pm, nw), :mu_sm_branch)
        _PMD.con(pm, nw)[:mu_sm_branch] = Dict{Tuple{Int,Int,Int}, Vector{JuMP.ConstraintRef}}()
    end

    rating = calc_branch_power_max(branch, bus, _PMD.ref(pm, nw, :total_real_load), _PMD.ref(pm, nw, :total_reactive_load))
    constraint_mc_thermal_limit_to_damaged(pm, nw, t_idx, branch["t_connections"], rating)
    nothing
end


"""
    constraint_mc_thermal_limit_from_ne(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for branch thermal constraints (from-side) for expansion lines
"""
function constraint_mc_thermal_limit_from_ne(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = _PMD.ref(pm, nw, :branch_ne, i)
    f_idx = (i, branch["f_bus"], branch["t_bus"])
    bus = _PMD.ref(pm,nw,:bus,branch["f_bus"])

    # using the same constraint store since a branch is either damaged or not damaged
    if !haskey(_PMD.con(pm, nw), :mu_sm_branch_ne)
        _PMD.con(pm, nw)[:mu_sm_branch_ne] = Dict{Tuple{Int,Int,Int}, Vector{JuMP.ConstraintRef}}()
    end

    rating = calc_branch_power_max(branch, bus, _PMD.ref(pm, nw, :total_real_load), _PMD.ref(pm, nw, :total_reactive_load))
    constraint_mc_thermal_limit_from_ne(pm, nw, f_idx, branch["f_connections"], rating)
    nothing
end


"""
    constraint_mc_thermal_limit_to_ne(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for branch thermal constraints (to-side) for expansion lines
"""
function constraint_mc_thermal_limit_to_ne(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = _PMD.ref(pm, nw, :branch_ne, i)
    t_idx = (i, branch["t_bus"], branch["f_bus"])
    bus = _PMD.ref(pm,nw,:bus,branch["t_bus"])

    # using the same constraint store since a branch is either damaged or not damaged
    if !haskey(_PMD.con(pm, nw), :mu_sm_branch_ne)
        _PMD.con(pm, nw)[:mu_sm_branch_ne] = Dict{Tuple{Int,Int,Int}, Vector{JuMP.ConstraintRef}}()
    end

    rating = calc_branch_power_max(branch, bus, _PMD.ref(pm, nw, :total_real_load), _PMD.ref(pm, nw, :total_reactive_load))
    constraint_mc_thermal_limit_to_ne(pm, nw, t_idx, branch["t_connections"], rating)
    nothing
end


"""
    constraint_mc_ampacity_from_damaged(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for branch current limit constraint from-side for branches that are damaged
"""
function constraint_mc_ampacity_from_damaged(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = _PMD.ref(pm, nw, :branch, i)
    f_idx = (i, branch["f_bus"], branch["t_bus"])
    bus = _PMD.ref(pm,nw,:bus,branch["f_bus"])

    if !haskey(_PMD.con(pm, nw), :mu_cm_branch)
        _PMD.con(pm, nw)[:mu_cm_branch] = Dict{Tuple{Int,Int,Int}, Vector{JuMP.ConstraintRef}}()
    end

    c_rating = calc_branch_current_max(branch, bus, _PMD.ref(pm, nw, :total_real_load), _PMD.ref(pm, nw, :total_reactive_load))
    constraint_mc_ampacity_from_damaged(pm, nw, f_idx, branch["f_connections"], c_rating)
    nothing
end


"""
    constraint_mc_ampacity_to_damaged(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for branch current limit constraint to-side for branches that are damaged
"""
function constraint_mc_ampacity_to_damaged(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = _PMD.ref(pm, nw, :branch, i)
    t_idx = (i, branch["t_bus"], branch["f_bus"])
    bus = _PMD.ref(pm,nw,:bus,branch["t_bus"])

    if !haskey(_PMD.con(pm, nw), :mu_cm_branch)
        _PMD.con(pm, nw)[:mu_cm_branch] = Dict{Tuple{Int,Int,Int}, Vector{JuMP.ConstraintRef}}()
    end

    c_rating = calc_branch_current_max(branch, bus, _PMD.ref(pm, nw, :total_real_load), _PMD.ref(pm, nw, :total_reactive_load))
    constraint_mc_ampacity_to_damaged(pm, nw, t_idx, branch["t_connections"], c_rating)
    nothing
end



"""
    constraint_mc_ampacity_from_ne(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for branch current limit constraint from-side for branches that are expansion
"""
function constraint_mc_ampacity_from_ne(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = _PMD.ref(pm, nw, :branch_ne, i)
    f_idx = (i, branch["f_bus"], branch["t_bus"])
    bus = _PMD.ref(pm,nw,:bus,branch["f_bus"])

    if !haskey(_PMD.con(pm, nw), :mu_cm_branch_ne)
        _PMD.con(pm, nw)[:mu_cm_branch_ne] = Dict{Tuple{Int,Int,Int}, Vector{JuMP.ConstraintRef}}()
    end

    c_rating = calc_branch_current_max(branch, bus, _PMD.ref(pm, nw, :total_real_load), _PMD.ref(pm, nw, :total_reactive_load))
    constraint_mc_ampacity_from_ne(pm, nw, f_idx, branch["f_connections"], c_rating)
    nothing
end


"""
    constraint_mc_ampacity_to_ne(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for branch current limit constraint to-side for branches that are expansion
"""
function constraint_mc_ampacity_to_ne(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = _PMD.ref(pm, nw, :branch_ne, i)
    t_idx = (i, branch["t_bus"], branch["f_bus"])
    bus = _PMD.ref(pm,nw,:bus,branch["t_bus"])

    if !haskey(_PMD.con(pm, nw), :mu_cm_branch_ne)
        _PMD.con(pm, nw)[:mu_cm_branch] = Dict{Tuple{Int,Int,Int}, Vector{JuMP.ConstraintRef}}()
    end

    c_rating = calc_branch_current_max(branch, bus, _PMD.ref(pm, nw, :total_real_load), _PMD.ref(pm, nw, :total_reactive_load))
    constraint_mc_ampacity_to_ne(pm, nw, t_idx, branch["t_connections"], c_rating)
    nothing
end


"""
    constraint_mc_voltage_angle_difference_damaged(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for constraints of the voltage angle difference across branches that are damaged
"""
function constraint_mc_voltage_angle_difference_damaged(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = _PMD.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    constraint_mc_voltage_angle_difference_damaged(pm, nw, f_idx, branch["f_connections"], branch["t_connections"], branch["angmin"], branch["angmax"])
    nothing
end


"""
    constraint_mc_voltage_angle_difference_damaged(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for constraints of the voltage angle difference across branches that are expansion
"""
function constraint_mc_voltage_angle_difference_ne(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = _PMD.ref(pm, nw, :branch_ne, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    constraint_mc_voltage_angle_difference_ne(pm, nw, f_idx, branch["f_connections"], branch["t_connections"], branch["angmin"], branch["angmax"])
    nothing
end



# Branch constraints

"""
    constraint_mc_ohms_yt_from(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)

Template function for ohms constraint for branches on the from-side for damaged lines
"""
function constraint_mc_ohms_yt_from_damaged(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)
    branch = _PMD.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    vad_min = _PMD.ref(pm, nw, :off_angmin)
    vad_max = _PMD.ref(pm, nw, :off_angmax)

    if all(all(isapprox.(branch[k], 0.0)) for k in ["br_r", "br_x", "g_fr", "g_to", "b_fr", "b_to"])
        @debug "branch $(branch["source_id"]) being treated as superconducting (effective zero impedance)"
        if !haskey(_PMD.con(pm, nw), :branch_flow)
            _PMD.con(pm, nw)[:branch_flow] = Dict{Int,Vector{Vector{<:JuMP.ConstraintRef}}}()
        end
        constraint_mc_branch_flow_damaged(pm, nw, f_idx, t_idx, branch["f_connections"], branch["t_connections"])
    else
        if !haskey(_PMD.con(pm, nw), :ohms_yt_from)
            _PMD.con(pm, nw)[:ohms_yt] = Dict{Tuple{Int,Int,Int},Vector{Vector{<:JuMP.ConstraintRef}}}()
        end
        G, B = _PMD.calc_branch_y(branch)
        constraint_mc_ohms_yt_from_damaged(pm, nw, f_bus, t_bus, f_idx, t_idx, branch["f_connections"], branch["t_connections"], G, B, branch["g_fr"], branch["b_fr"], vad_min, vad_max)
    end
end


"""
    constraint_mc_ohms_yt_to(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)

Template function for ohms constraint for branches on the to-side for damaged lines
"""
function constraint_mc_ohms_yt_to_damaged(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)
    branch = _PMD.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    vad_min = _PMD.ref(pm, nw, :off_angmin)
    vad_max = _PMD.ref(pm, nw, :off_angmax)

    if all(all(isapprox.(branch[k], 0.0)) for k in ["br_r", "br_x", "g_fr", "g_to", "b_fr", "b_to"])
        @debug "branch $(branch["source_id"]) being treated as superconducting (effective zero impedance)"
        if !haskey(_PMD.con(pm, nw), :branch_flow)
            _PMD.con(pm, nw)[:branch_flow] = Dict{Int,Vector{Vector{<:JuMP.ConstraintRef}}}()
        end
        constraint_mc_branch_flow_damaged(pm, nw, f_idx, t_idx, branch["f_connections"], branch["t_connections"])
    else
        if !haskey(_PMD.con(pm, nw), :ohms_yt_to)
            _PMD.con(pm, nw)[:ohms_yt] = Dict{Tuple{Int,Int,Int},Vector{Vector{<:JuMP.ConstraintRef}}}()
        end
        G, B = _PMD.calc_branch_y(branch)
        constraint_mc_ohms_yt_to_damaged(pm, nw, f_bus, t_bus, f_idx, t_idx, branch["f_connections"], branch["t_connections"], G, B, branch["g_to"], branch["b_to"], vad_min, vad_max)
    end
end


"""
    constraint_mc_ohms_yt_from(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)

Template function for ohms constraint for branches on the from-side for ne lines
"""
function constraint_mc_ohms_yt_from_ne(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)
    branch = _PMD.ref(pm, nw, :branch_ne, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    vad_min = _PMD.ref(pm, nw, :off_angmin)
    vad_max = _PMD.ref(pm, nw, :off_angmax)

    if all(all(isapprox.(branch[k], 0.0)) for k in ["br_r", "br_x", "g_fr", "g_to", "b_fr", "b_to"])
        @debug "branch $(branch["source_id"]) being treated as superconducting (effective zero impedance)"
        if !haskey(_PMD.con(pm, nw), :branch_flow)
            _PMD.con(pm, nw)[:branch_flow] = Dict{Int,Vector{Vector{<:JuMP.ConstraintRef}}}()
        end
        constraint_mc_branch_flow_damaged(pm, nw, f_idx, t_idx, branch["f_connections"], branch["t_connections"])
    else
        if !haskey(_PMD.con(pm, nw), :ohms_yt_from)
            _PMD.con(pm, nw)[:ohms_yt] = Dict{Tuple{Int,Int,Int},Vector{Vector{<:JuMP.ConstraintRef}}}()
        end
        G, B = _PMD.calc_branch_y(branch)
        constraint_mc_ohms_yt_from_ne(pm, nw, f_bus, t_bus, f_idx, t_idx, branch["f_connections"], branch["t_connections"], G, B, branch["g_fr"], branch["b_fr"], vad_min, vad_max)
    end
end


"""
    constraint_mc_ohms_yt_to(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)

Template function for ohms constraint for branches on the to-side for ne lines
"""
function constraint_mc_ohms_yt_to_ne(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)
    branch = _PMD.ref(pm, nw, :branch_ne, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    vad_min = _PMD.ref(pm, nw, :off_angmin)
    vad_max = _PMD.ref(pm, nw, :off_angmax)

    if all(all(isapprox.(branch[k], 0.0)) for k in ["br_r", "br_x", "g_fr", "g_to", "b_fr", "b_to"])
        @debug "branch $(branch["source_id"]) being treated as superconducting (effective zero impedance)"
        if !haskey(_PMD.con(pm, nw), :branch_flow)
            _PMD.con(pm, nw)[:branch_flow] = Dict{Int,Vector{Vector{<:JuMP.ConstraintRef}}}()
        end
        constraint_mc_branch_flow_damaged(pm, nw, f_idx, t_idx, branch["f_connections"], branch["t_connections"])
    else
        if !haskey(_PMD.con(pm, nw), :ohms_yt_to)
            _PMD.con(pm, nw)[:ohms_yt] = Dict{Tuple{Int,Int,Int},Vector{Vector{<:JuMP.ConstraintRef}}}()
        end
        G, B = _PMD.calc_branch_y(branch)
        constraint_mc_ohms_yt_to_ne(pm, nw, f_bus, t_bus, f_idx, t_idx, branch["f_connections"], branch["t_connections"], G, B, branch["g_to"], branch["b_to"], vad_min, vad_max)
    end
end


"""
    constraint_mc_gen_power_ne(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for ne generator power on/off constraints (MLD problems)
"""
function constraint_mc_gen_power_ne(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    gen = _PMD.ref(pm, nw, :gen_ne, i)
    ncnds = length(gen["connections"])

    pmin = get(gen, "pmin", fill(-Inf, ncnds))
    pmax = get(gen, "pmax", fill( Inf, ncnds))
    qmin = get(gen, "qmin", fill(-Inf, ncnds))
    qmax = get(gen, "qmax", fill( Inf, ncnds))

    constraint_mc_gen_power_ne(pm, nw, i, gen["connections"], pmin, pmax, qmin, qmax)
    nothing
end


""""
    constraint_radial_topology(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, relax::Bool=false)

Template function radial topology constraint with expansion edges.
"""
function constraint_radial_topology_ne(pm::_PMD.AbstractUnbalancedPowerModel; nw::Int=nw_id_default, relax::Bool=false)
    constraint_radial_topology_ne(pm, nw; relax=relax)
end



#function constraint_mc_vm_vuf(pm::_PMs.AbstractPowerModel, bus_id::Int; nw::Int=pm.cnw)
#    bus = _PMs.ref(pm, nw, :bus_bal)[bus_id]
#    vufmax = _PMs.ref(pm, nw, :bus)[bus]["vm_vuf_max"]
#    constraint_mc_vm_vuf(pm, nw, bus, vufmax)
#end

#function constraint_mc_vm_vuf_mod(pm::_PMs.AbstractPowerModel, bus_id::Int; nw::Int=pm.cnw)
#    bus = _PMs.ref(pm, nw, :bus_bal)[bus_id]
#    vufmax = _PMs.ref(pm, nw, :bus)[bus]["vm_vuf_max"]
#    constraint_mc_vm_vuf(pm, nw, bus, vufmax)
#end

#function constraint_switch(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
#    constraint_switch_damaged(pm, nw)
#    constraint_switch_new(pm, nw)
#end

#function constraint_switch_damaged(pm::_PMs.AbstractPowerModel, nw::Int)
#    for (l,i,j) in _PMs.ref(pm, nw, :arcs_damaged)
#        xe = _PMs.var(pm, nw, :xe_d_s, (l,i,j))
#        te = _PMs.var(pm, nw, :te_d_s, (l,i,j))
#        constraint_switch(pm, xe, te)
#    end
#end

#function constraint_switch_new(pm::_PMs.AbstractPowerModel, nw::Int)
#    for (l,i,j) in _PMs.ref(pm, nw, :arcs_new)
#        xe = _PMs.var(pm, nw, :xe_n_s, (l,i,j))
#        te = _PMs.var(pm, nw, :te_n_s, (l,i,j))
#        constraint_switch(pm, xe, te)
#    end
#end

#function constraint_harden(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
#    constraint_harden_damaged(pm, nw)
#    constraint_harden_new(pm, nw)
#end

#function constraint_harden_damaged(pm::_PMs.AbstractPowerModel, nw::Int)
#    for (l,i,j) in _PMs.ref(pm, nw, :arcs_damaged)
#        xe = _PMs.var(pm, nw, :xe_d_s, (l,i,j))
#        he = _PMs.var(pm, nw, :he_d_s, (l,i,j))
#        constraint_harden(pm, xe, he)
#    end
#end

#function constraint_harden_new(pm::_PMs.AbstractPowerModel, nw::Int)
#    for (l,i,j) in _PMs.ref(pm, nw, :arcs_new)
#        xe = _PMs.var(pm, nw, :xe_n_s, (l,i,j))
#        he = _PMs.var(pm, nw, :he_n_s, (l,i,j))
#        constraint_harden(pm, xe, he)
#    end
#end

#function constraint_active_line(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
#    constraint_active_line_existing(pm, nw)
#    constraint_active_line_new(pm, nw)
#    constraint_active_line_active(pm, nw)
#end

#function constraint_active_line_existing(pm::_PMs.AbstractPowerModel, nw::Int)
#    for (l,i,j) in _PMs.ref(pm, nw, :arcs_damaged)
#        constraint_activation_damage(pm, nw, (l,i,j))
#    end
#end

#function constraint_active_line_new(pm::_PMs.AbstractPowerModel, nw::Int)
#    for (l,i,j) in _PMs.ref(pm, nw, :arcs_new)
#        constraint_activation_new(pm, nw, (l,i,j))
#    end
#end

#function constraint_active_line_active(pm::_PMs.AbstractPowerModel, nw::Int)
#    for (l,i,j) in _PMs.ref(pm, nw, :arcs)
#        (l,i,j) in _PMs.ref(pm, nw, :arcs_damaged) || (l,i,j) in _PMs.ref(pm, nw, :arcs_new) ? nothing : constraint_activation_active(pm, nw, (l,i,j))
#    end
#end

#function constraint_mc_ohms_yt_from(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
#    branch = _PMs.ref(pm, nw, :branch, i)
#    f_bus = branch["f_bus"]
#    t_bus = branch["t_bus"]
#    f_idx = (i, f_bus, t_bus)
#    t_idx = (i, t_bus, f_bus)

#    g, b = _PMs.calc_branch_y(branch)
#    tr, ti = _PMs.calc_branch_t(branch)
#    g_fr = branch["g_fr"]
#    b_fr = branch["b_fr"]
#    tm = branch["tap"]

#    constraint_mc_ohms_yt_from(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
#end


#"ohms constraint for branches on the to-side"
#function constraint_mc_ohms_yt_to(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
#    branch = _PMs.ref(pm, nw, :branch, i)
#    f_bus = branch["f_bus"]
#    t_bus = branch["t_bus"]
#    f_idx = (i, f_bus, t_bus)
#    t_idx = (i, t_bus, f_bus)

#    g, b = _PMs.calc_branch_y(branch)
#    tr, ti = _PMs.calc_branch_t(branch)
#    g_to = branch["g_to"]
#    b_to = branch["b_to"]
#    tm = branch["tap"]

#    constraint_mc_ohms_yt_to(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
#end

#function constraint_variable(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
#    constraint_variable_he(pm, nw)
#    constraint_variable_te(pm, nw)
#    constraint_variable_xe(pm, nw)
#end

#function constraint_variable_he(pm::_PMs.AbstractPowerModel, nw::Int=pm.cnw)
#    for (l,i,j) in _PMs.ref(pm, nw, :arcs_damaged)
#        constraint_variable_he_d(pm, nw, (l,i,j))
#    end
#    for (l,i,j) in _PMs.ref(pm, nw, :arcs_new)
#        constraint_variable_he_n(pm, nw, (l,i,j))
#    end
#end

#function constraint_variable_te(pm::_PMs.AbstractPowerModel, nw::Int=pm.cnw)
#    for (l,i,j) in _PMs.ref(pm, nw, :arcs_damaged)
#        constraint_variable_te_d(pm, nw, (l,i,j))
#    end
#    for (l,i,j) in _PMs.ref(pm, nw, :arcs_new)
#        constraint_variable_te_n(pm, nw, (l,i,j))
#    end
#end

#function constraint_variable_xe(pm::_PMs.AbstractPowerModel, nw::Int=pm.cnw)
#    for (l,i,j) in _PMs.ref(pm, nw, :arcs_damaged)
#        constraint_variable_xe_d(pm, nw, (l,i,j))
#    end
#    for (l,i,j) in _PMs.ref(pm, nw, :arcs_new)
#        constraint_variable_xe_n(pm, nw, (l,i,j))
#    end
#end

#function constraint_mc_voltage_angle_difference(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
#    branch = _PMs.ref(pm, nw, :branch, i)
#    f_bus = branch["f_bus"]
#    t_bus = branch["t_bus"]
#    f_idx = (i, f_bus, t_bus)
#    pair = (f_bus, t_bus)

#    constraint_mc_voltage_angle_difference(pm, nw, f_idx, branch["angmin"], branch["angmax"])
#end

#function constraint_cycle_function(pm::_PMs.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
#    branch = _PMs.ref(pm, nw, :branch, i)
#    f_bus = branch["f_bus"]
#    t_bus = branch["t_bus"]
#    arcs = (i, f_bus, t_bus)

#    constraint_cycle_function(pm, nw, arcs)
#end

#function constraint_cycle_elimination(pm, nw::Int, tour::Int)
#    tours = _PMs.ref(pm, nw, :arc_tour, tour)

#    constraint_cycle_elimination(pm, nw::Int, tours)
#end

# Switch constraints

"""
    constraint_mc_switch_state_open_close_inline_ne(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)


"""
function constraint_mc_switch_state_open_close_inline_ne(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)
    switch = _PMD.ref(pm, nw, :switch_inline_ne, i)

    f_bus = switch["f_bus"]
    t_bus = switch["t_bus"]

    f_connections = switch["f_connections"]
    t_connections = switch["t_connections"]

    constraint_mc_switch_voltage_open_close_inline_ne(pm, nw, i, f_bus, t_bus, f_connections, t_connections)
#    constraint_mc_switch_power_open_close_inline_ne(pm, nw, i, f_bus, t_bus, f_connections, t_connections)
end
