@doc raw"""
    constraint_mc_switch_thermal_limit(pm::AbstractUnbalancedActivePowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rating::Vector{<:Real})::Nothing

Active power only switch thermal limit constraint

math```
-S_{max} \leq p_{fr} \leq S_{max}
```
"""
function constraint_mc_switch_inline_ne_thermal_limit(pm::_PMD.AbstractUnbalancedActivePowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rating::Vector{<:Real})::Nothing
    psw = _PMD.var(pm, nw, :psw_inline_ne, f_idx)

    mu_sm_fr = JuMP.ConstraintRef[]
    for (idx, c) in enumerate(f_connections)
        if rating[idx] < Inf
            JuMP.lower_bound(psw[c]) < -rating[idx] && set_lower_bound(psw[c], -rating[idx])
            JuMP.upper_bound(psw[c]) >  rating[idx] && set_upper_bound(psw[c],  rating[idx])
        end
    end

    _PMD.con(pm, nw, :mu_sm_switch_inline_ne)[f_idx] = mu_sm_fr

    nothing
end
