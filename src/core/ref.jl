
#function ref_add_critical_leve!(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, limit::Float64 = .90)
#    _PMs.ref(pm, nw)[:critical_level] = limit
#end

#function ref_add_demand_level!(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, limit::Float64 = .80)
#    _PMs.ref(pm, nw)[:demand_level] = limit
#end

#function ref_add_vm_imbalance!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
#    if _IM.ismultinetwork(data)
#        nws_data = data["nw"]
#    else
#        nws_data = Dict("0" => data)
#    end
#    for (n, nw_data) in nws_data
#        nw_id = parse(Int, n)
#        nw_ref = ref[:nw][nw_id]
#        nw_ref[:bus_bal] = []
#        for (i, bus) in nw_data["bus"]
#            haskey(bus, "vm_vuf_max") ? append!(nw_ref[:bus_bal], parse(Int,i)) : nothing
#        end
#    end
#end

#function ref_add_pq_imbalance!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
#    if _IM.ismultinetwork(data)
#        nws_data = data["nw"]
#    else
#        nws_data = Dict("0" => data)
#    end
#    for (n, nw_data) in nws_data
#        nw_id = parse(Int, n)
#        nw_ref = ref[:nw][nw_id]
#        nw_ref[:arcs_bal] = Tuple{Int64,Int64, Int64}[]
#        for (i, branch) in nw_data["branch"]
#            if haskey(branch, "pq_imbalance")
#                j = branch["f_bus"]
#                k = branch["t_bus"]
#                push!(nw_ref[:arcs_bal], (parse(Int,i),j,k))
#            end
#        end
#        for (i,trans) in nw_data["transformer"]
#            if haskey(trans, "pq_imbalance")
#                j = trans["f_bus"]
#                k = trans["t_bus"]
#                push!(nw_ref[:arcs_bal], (parse(Int,i),j,k))
#            end
#        end
#    end
#end

function ref_add_branch_ne!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    for (nw, nw_ref) in ref[:it][_PMD.pmd_it_sym][:nw]
        ### filter out inactive components ###
        nw_ref[:branch_ne] = Dict(x for x in get(nw_ref, :branch_ne, Dict()) if (x.second["br_status"] != 0 && x.second["f_bus"] in keys(nw_ref[:bus]) && x.second["t_bus"] in keys(nw_ref[:bus])))

        ### setup arcs from edges ###
        nw_ref[:arcs_branch_ne_from] = [(i,branch["f_bus"],branch["t_bus"]) for (i,branch) in nw_ref[:branch_ne]]
        nw_ref[:arcs_branch_ne_to]   = [(i,branch["t_bus"],branch["f_bus"]) for (i,branch) in nw_ref[:branch_ne]]
        nw_ref[:arcs_branch_ne]      = [nw_ref[:arcs_branch_ne_from]; nw_ref[:arcs_branch_ne_to]]

        ### bus connected component lookups ###
        bus_arcs = Dict((i, Tuple{Int,Int,Int}[]) for (i,bus) in nw_ref[:bus])
        for (l,i,j) in nw_ref[:arcs_branch_ne]
            push!(bus_arcs[i], (l,i,j))
        end
        nw_ref[:bus_arcs_branch_ne] = bus_arcs

        ### connections
        conns = Dict{Int,Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}}([(i, []) for (i, bus) in nw_ref[:bus]])
        for (i, obj) in nw_ref[:branch_ne]
            push!(conns[obj["f_bus"]], ((obj["index"], obj["f_bus"], obj["t_bus"]), obj["f_connections"]))
            if obj["f_bus"] != obj["t_bus"]
                push!(conns[obj["t_bus"]], ((obj["index"], obj["t_bus"], obj["f_bus"]), obj["t_connections"]))
            end
        end
        nw_ref[:bus_arcs_conns_branch_ne] = conns

        ### aggregate info for pairs of connected buses ###
        if !haskey(nw_ref, :buspairs_ne)
            nw_ref[:buspairs_ne] = _PMD.calc_buspair_parameters(nw_ref[:bus], nw_ref[:branch_ne])
        end
    end
end

function ref_add_switch_inline_ne!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    for (nw, nw_ref) in ref[:it][_PMD.pmd_it_sym][:nw]
        ### filter out inactive components ###
        nw_ref[:switch_inline_ne] = Dict(x for x in get(nw_ref, :switch_inline_ne, Dict()) if (x.second["status"] != 0 && x.second["f_bus"] in keys(nw_ref[:bus]) && x.second["t_bus"] in keys(nw_ref[:bus])))

        ### setup arcs from edges ###
        nw_ref[:arcs_switch_inline_ne_from] = [(i,branch["f_bus"],branch["t_bus"]) for (i,branch) in nw_ref[:switch_inline_ne]]
        nw_ref[:arcs_switch_inline_ne_to]   = [(i,branch["t_bus"],branch["f_bus"]) for (i,branch) in nw_ref[:switch_inline_ne]]
        nw_ref[:arcs_switch_inline_ne]      = [nw_ref[:arcs_switch_inline_ne_from]; nw_ref[:arcs_switch_inline_ne_to]]

        ### bus connected component lookups ###
        bus_arcs = Dict((i, Tuple{Int,Int,Int}[]) for (i,bus) in nw_ref[:bus])
        for (l,i,j) in nw_ref[:arcs_switch_inline_ne]
            push!(bus_arcs[i], (l,i,j))
        end
        nw_ref[:bus_arcs_switch_inline_ne] = bus_arcs

        ### connections
        conns = Dict{Int,Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}}([(i, []) for (i, bus) in nw_ref[:bus]])
        for (i, obj) in nw_ref[:switch_inline_ne]
            push!(conns[obj["f_bus"]], ((obj["index"], obj["f_bus"], obj["t_bus"]), obj["f_connections"]))
            if obj["f_bus"] != obj["t_bus"]
                push!(conns[obj["t_bus"]], ((obj["index"], obj["t_bus"], obj["f_bus"]), obj["t_connections"]))
            end
        end
        nw_ref[:bus_arcs_conns_switch_inline_ne] = conns
    end
end

function ref_add_gen_ne!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    for (nw, nw_ref) in ref[:it][_PMD.pmd_it_sym][:nw]
        nw_ref[:gen_ne] = Dict(x for x in get(nw_ref, :gen_ne, Dict()) if (x.second["gen_status"] != 0 && x.second["gen_bus"] in keys(nw_ref[:bus])))

        bus_objs = Dict((i, Int[]) for (i,bus) in nw_ref[:bus])
        for (i, obj) in nw_ref[:gen_ne]
            push!(bus_objs[obj["gen_bus"]], i)
        end
        nw_ref[Symbol("bus_gen_nes")] = bus_objs

        conns = Dict{Int,Vector{Tuple{Int,Vector{Int}}}}([(i, []) for (i, bus) in nw_ref[:bus]])
        for (i, obj) in nw_ref[:gen_ne]
            if obj["gen_status"] != 0
                push!(conns[obj["gen_bus"]], (i, obj["connections"]))
            end
        end
        nw_ref[Symbol("bus_conns_gen_ne")] = conns
    end
end


function ref_add_transformer_ne!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    for (nw, nw_ref) in ref[:it][_PMD.pmd_it_sym][:nw]
        ### filter out inactive components ###
        nw_ref[:transformer_ne] = Dict(x for x in get(nw_ref, :transformer_ne, Dict()) if (x.second["br_status"] != 0 && x.second["f_bus"] in keys(nw_ref[:bus]) && x.second["t_bus"] in keys(nw_ref[:bus])))

        ### setup arcs from edges ###
        nw_ref[:arcs_transformer_ne_from] = [(i,branch["f_bus"],branch["t_bus"]) for (i,branch) in nw_ref[:transformer_ne]]
        nw_ref[:arcs_transformer_ne_to]   = [(i,branch["t_bus"],branch["f_bus"]) for (i,branch) in nw_ref[:transformer_ne]]
        nw_ref[:arcs_transformer_ne]      = [nw_ref[:arcs_transformer_ne_from]; nw_ref[:arcs_transformer_ne_to]]

        ### bus connected component lookups ###
        bus_arcs = Dict((i, Tuple{Int,Int,Int}[]) for (i,bus) in nw_ref[:bus])
        for (l,i,j) in nw_ref[:arcs_transformer_ne]
            push!(bus_arcs[i], (l,i,j))
        end
        nw_ref[:bus_arcs_transformer_ne] = bus_arcs

        ### connections
        conns = Dict{Int,Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}}([(i, []) for (i, bus) in nw_ref[:bus]])
        for (i, obj) in nw_ref[:transformer_ne]
            push!(conns[obj["f_bus"]], ((obj["index"], obj["f_bus"], obj["t_bus"]), obj["f_connections"]))
            if obj["f_bus"] != obj["t_bus"]
                push!(conns[obj["t_bus"]], ((obj["index"], obj["t_bus"], obj["f_bus"]), obj["t_connections"]))
            end
        end
        nw_ref[:bus_arcs_conns_transformer_ne] = conns
    end
end


function ref_add_branch_harden!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    for (nw, nw_ref) in ref[:it][_PMD.pmd_it_sym][:nw]
        nw_ref[:branch_harden] = [i for (i,branch) in nw_ref[:branch] if branch["can_harden"] == true]
    end
end


function ref_add_undamaged_branch!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    for (nw, nw_ref) in ref[:it][_PMD.pmd_it_sym][:nw]
        nw_ref[:undamaged_branch] = [i for (i,branch) in nw_ref[:branch] if !(i in nw_ref[:damaged_branch])]
    end
end


# functions to enumerate sub tours:
#       - genSubTour: generate a list of cycles, each item as a set with the buses cycled through
#       - findCycle: with two paths forming a cycle, extracting the cycle
#       - node2arcs: convert the output of genSubTour to a list of cycles with arcs
#function findCycle(path1, path2)
    # once we get two paths that merge, get the cycle part
#    sInd = 1;
#    sameBool = true;
#    while sameBool
#        if path1[sInd] == path2[sInd]
#            if sInd < min(length(path1),length(path2))
#                sInd += 1;
#            else
#                sameBool = false;
#            end
#        else
#            sameBool = false;
#            sInd -= 1;
#        end
#    end
#    cycle = [path1[sInd]];
#    for i1 in (sInd+1):length(path1)
#        push!(cycle,path1[i1]);
#    end
#    for i2 in (sInd+1):length(path2)
#        push!(cycle,path2[i2]);
#    end
#    return unique(cycle);
#end

#function genSubTour(ref::Dict{Symbol,<:Any})
    # start searching from each node
    # get a list of all buses
    # all_bus = keys(pm, :bus)
    # println(all_bus)
#    cycle_list = Dict();
#    tourList = [];
#    for i in keys(ref[:bus])
#        cycle_list[i] = [];
#        marked_bus = [i];
#        history_bus = Dict();
#        from_bus = Dict();
#        for j in keys(ref[:bus])
#            history_bus[j] = [];
#            from_bus[j] = -1;
#        end
#        history_bus[i] = [i];
#        busList = [i];
        # BFS to get all cycles starts/ends at bus i
#        while busList != []
            # get the adjacent buses of current_bus
#            current_bus = busList[1];
            # println(current_bus, " ", from_bus[current_bus]);
            # get the adjacent buses of the current bus
#            adjacent_bus = unique([arcs[3] for arcs in ref[:bus_arcs][current_bus]]);
#            history_current = history_bus[current_bus];
#            for j in adjacent_bus
#                if (!(j in marked_bus))
#                    push!(marked_bus,j);
#                    push!(busList,j);
#                    from_bus[j] = current_bus;
#                    history_j = copy(history_current);
#                    history_bus[j] = push!(history_j,j);
#                elseif j != from_bus[current_bus]
                    # record the cycle
#                    cycle = Set(findCycle(history_bus[current_bus],history_bus[j]));
#                    if !(cycle in cycle_list[i])
#                        push!(cycle_list[i],cycle);
#                    end
#                end
#            end
#            deleteat!(busList,1);
#        end
#        append!(tourList,cycle_list[i]);
#    end
#    return unique(tourList);
#end

#function node2arcs(cycle, ref::Dict{Symbol,<:Any})
    # convert the output from genSubTour (node cycle) to arc cyle
#    nodeList = sort([item for item in cycle]);
    # get a dictionary with each pair of adjacent nodes
#    arcDict = Dict();
#    for i in 1:length(nodeList)
#        for j in i+1:length(nodeList)
#            arcDict[nodeList[i],nodeList[j]] = [a for a in ref[:bus_arcs][i] if j == a[3]];
#        end
#    end
#    arcIter = [];
#    for aKey in keys(arcDict)
#        if arcDict[aKey] != []
#            push!(arcIter,arcDict[aKey])
#        end
#    end
#    arc_cycle_list = [];
#    for i in Iterators.product(arcIter...)
#        push!(arc_cycle_list,i);
#    end

#    return arc_cycle_list;
#end


#function ref_add_subtour!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
#    if _IM.ismultinetwork(data)
#        nws_data = data["nw"]
#    else
#        nws_data = Dict("0" => data)
#    end
#    for (n, nw_data) in nws_data
#        nw_id = parse(Int, n)
#        nw_ref = ref[:nw][nw_id]
#        tourList = genSubTour(nw_ref)
#        nw_ref[:arc_tour] = []
#        for tour in tourList
#            hold = []
#            for arc_tour in node2arcs(tour, nw_ref)
#                append!(hold, arc_tour)
#            end
#            push!(nw_ref[:arc_tour], hold)
#        end
#    end
#end

function ref_add_rdt!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    ref_add_branch_ne!(ref, data)
    ref_add_transformer_ne!(ref, data)
    ref_add_branch_harden!(ref, data)
    ref_add_gen_ne!(ref, data)
    ref_add_switch_inline_ne!(ref,data)
    ref_add_undamaged_branch!(ref,data)
#    ref_add_vm_imbalance!(ref, data)
#    ref_add_pq_imbalance!(ref, data)
#    ref_add_subtour!(ref, data)
end
