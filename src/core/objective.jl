
function objective(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw)
    harden_cost = []
    JuMP.@objective(pm.model, Min, 
                    sum((pm.ref[:nw][1][:branch][l]["harden_cost"] * _PMs.var(pm, nw, :he_d, (l,i,j))) for (l,i,j) in pm.ref[:arcs_damaged_all]) # harden cost 
                  + sum((pm.ref[:nw][1][:branch][l]["construction_cost"] * _PMs.var(pm, nw, :he_n, (l,i,j))) for (l,i,j) in pm.ref[:arcs_new_all]) # based off 6c
    )
end
