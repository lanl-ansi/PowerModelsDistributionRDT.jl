
function objective_rdt(pm::_PMD.AbstractUnbalancedPowerModel; nw::Int=_PMD.nw_id_default)
    JuMP.@objective(pm.model, Min,
                    sum(_PMD.ref(pm, nw, :branch, l)["harden_cost"] * _PMD.var(pm, nw, :he, l) for l in _PMD.ref(pm,nw,:branch_harden)) # harden cost
                 +  sum(branch["construction_cost"] * _PMD.var(pm, nw, :xe, l) for (l, branch) in _PMD.ref(pm,nw,:branch_ne)) # new line cost
                 +  sum(switch["switch_cost"] * _PMD.var(pm, nw, :te, l) for (l, switch) in _PMD.ref(pm,nw,:switch_inline_ne)) # new swith cost
                 +  sum(gen["microgrid_cost"] * _PMD.var(pm, nw, :ue, l) for (l, gen) in _PMD.ref(pm,nw,:gen_ne)) # new generator cost
    )
end
