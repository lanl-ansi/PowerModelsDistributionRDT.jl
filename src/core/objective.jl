
function objective_rdt(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default)
    JuMP.@objective(pm.model, Min,
        sum(ref(pm, nw, :branch, l)["harden_cost"] * var(pm, nw, :he, l) for l in ref(pm, nw, :branch_harden)) # harden cost
        + sum(branch["construction_cost"] * var(pm, nw, :xe, l) for (l, branch) in ref(pm, nw, :branch_ne)) # new line cost
        + sum(switch["switch_cost"] * var(pm, nw, :te, l) for (l, switch) in ref(pm, nw, :switch_inline_ne)) # new swith cost
        + sum(gen["microgrid_cost"] * var(pm, nw, :ue, l) for (l, gen) in ref(pm, nw, :gen_ne)) # new generator cost
    )
end
