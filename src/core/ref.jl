
function ref_add_critical_leve!(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, limit::Float64 = .90)
    _PMs.ref(pm, nw)[:critical_level] = limit
end

function ref_add_demand_level!(pm::_PMs.AbstractPowerModel; nw::Int=pm.cnw, limit::Float64 = .80)
    _PMs.ref(pm, nw)[:demand_level] = limit
end

function ref_add_vm_imbalance!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    if _INs.ismultinetwork(data)
        nws_data = data["nw"]
    else
        nws_data = Dict("0" => data)
    end
    for (n, nw_data) in nws_data
        nw_id = parse(Int, n)
        nw_ref = ref[:nw][nw_id]
        nw_ref[:bus_bal] = []
        for (i, bus) in nw_data["bus"]
            haskey(bus, "vm_vuf_max") ? append!(nw_ref[:bus_bal], parse(Int,i)) : nothing
        end
    end
end

function ref_add_pq_imbalance!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    if _INs.ismultinetwork(data)
        nws_data = data["nw"]
    else
        nws_data = Dict("0" => data)
    end
    for (n, nw_data) in nws_data
        nw_id = parse(Int, n)
        nw_ref = ref[:nw][nw_id]
        nw_ref[:arcs_bal] = Tuple{Int64,Int64, Int64}[]
        for (i, branch) in nw_data["branch"]
            if haskey(branch, "pq_imbalance")
                j = branch["f_bus"]
                k = branch["t_bus"]
                push!(nw_ref[:arcs_bal], (parse(Int,i),j,k))
            end
        end
        for (i,trans) in nw_data["transformer"]
            if haskey(trans, "pq_imbalance")
                j = trans["f_bus"]
                k = trans["t_bus"]
                push!(nw_ref[:arcs_bal], (parse(Int,i),j,k))
            end
        end
    end
end

function ref_add_damaged_lines!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    if _INs.ismultinetwork(data)
        nws_data = data["nw"]
    else
        nws_data = Dict("0" => data)
    end
    hold = Array{Tuple{Int64,Int64,Int64},1}()
    for (n, nw_data) in nws_data
        nw_id = parse(Int, n)
        nw_ref = ref[:nw][nw_id]
        nw_ref[:arcs_damaged] = Array{Tuple{Int64,Int64,Int64},1}()
        for i in nw_ref[:disabled_lines]
            branch = nw_ref[:branch][parse(Int, i)]
            branch["is_new"] ? nothing : push!(nw_ref[:arcs_damaged], (parse(Int, i),branch["t_bus"],branch["f_bus"]))
            if (parse(Int, i),branch["t_bus"],branch["f_bus"]) in hold
                nothing
            else
                branch["is_new"] ? nothing : push!(hold, (parse(Int, i),branch["t_bus"],branch["f_bus"]))
            end
        end
    end
    ref[:arcs_damaged_all] = hold
end

function ref_add_new_lines!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    if _INs.ismultinetwork(data)
        nws_data = data["nw"]
    else
        nws_data = Dict("0" => data)
    end
    hold = Array{Tuple{Int64,Int64,Int64},1}()
    for (n, nw_data) in nws_data
        nw_id = parse(Int, n)
        nw_ref = ref[:nw][nw_id]
        nw_ref[:arcs_new] = Array{Tuple{Int64,Int64,Int64},1}()
        for (i, branch) in nw_ref[:branch]
            if branch["is_new"]
                push!(nw_ref[:arcs_new], (i,branch["t_bus"],branch["f_bus"]))
                if (i,branch["t_bus"],branch["f_bus"]) in hold
                    nothing
                else
                   push!(hold, (i,branch["t_bus"],branch["f_bus"]))
                end
            end
        end
    end
    ref[:arcs_new_all] = hold
end

# functions to enumerate sub tours:
#       - genSubTour: generate a list of cycles, each item as a set with the buses cycled through
#       - findCycle: with two paths forming a cycle, extracting the cycle
#       - node2arcs: convert the output of genSubTour to a list of cycles with arcs
function findCycle(path1, path2)
    # once we get two paths that merge, get the cycle part
    sInd = 1;
    sameBool = true;
    while sameBool
        if path1[sInd] == path2[sInd]
            if sInd < min(length(path1),length(path2))
                sInd += 1;
            else
                sameBool = false;
            end
        else
            sameBool = false;
            sInd -= 1;
        end
    end
    cycle = [path1[sInd]];
    for i1 in (sInd+1):length(path1)
        push!(cycle,path1[i1]);
    end
    for i2 in (sInd+1):length(path2)
        push!(cycle,path2[i2]);
    end
    return unique(cycle);
end

function genSubTour(ref::Dict{Symbol,<:Any})
    # start searching from each node
    # get a list of all buses
    # all_bus = keys(pm, :bus)
    # println(all_bus)
    cycle_list = Dict();
    tourList = [];
    for i in keys(ref[:bus])
        cycle_list[i] = [];
        marked_bus = [i];
        history_bus = Dict();
        from_bus = Dict();
        for j in keys(ref[:bus])
            history_bus[j] = [];
            from_bus[j] = -1;
        end
        history_bus[i] = [i];
        busList = [i];
        # BFS to get all cycles starts/ends at bus i
        while busList != []
            # get the adjacent buses of current_bus
            current_bus = busList[1];
            # println(current_bus, " ", from_bus[current_bus]);
            # get the adjacent buses of the current bus
            adjacent_bus = unique([arcs[3] for arcs in ref[:bus_arcs][current_bus]]);
            history_current = history_bus[current_bus];
            for j in adjacent_bus
                if (!(j in marked_bus))
                    push!(marked_bus,j);
                    push!(busList,j);
                    from_bus[j] = current_bus;
                    history_j = copy(history_current);
                    history_bus[j] = push!(history_j,j);
                elseif j != from_bus[current_bus]
                    # record the cycle
                    cycle = Set(findCycle(history_bus[current_bus],history_bus[j]));
                    if !(cycle in cycle_list[i])
                        push!(cycle_list[i],cycle);
                    end
                end
            end
            deleteat!(busList,1);
        end
        append!(tourList,cycle_list[i]);
    end
    return unique(tourList);
end

function node2arcs(cycle, ref::Dict{Symbol,<:Any})
    # convert the output from genSubTour (node cycle) to arc cyle
    nodeList = sort([item for item in cycle]);
    # get a dictionary with each pair of adjacent nodes
    arcDict = Dict();
    for i in 1:length(nodeList)
        for j in i+1:length(nodeList)
            arcDict[nodeList[i],nodeList[j]] = [a for a in ref[:bus_arcs][i] if j == a[3]];
        end
    end
    arcIter = [];
    for aKey in keys(arcDict)
        if arcDict[aKey] != []
            push!(arcIter,arcDict[aKey])
        end
    end
    arc_cycle_list = []; 
    for i in Iterators.product(arcIter...)
        push!(arc_cycle_list,i);
    end

    return arc_cycle_list;
end


function ref_add_subtour!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    if _INs.ismultinetwork(data)
        nws_data = data["nw"]
    else
        nws_data = Dict("0" => data)
    end
    for (n, nw_data) in nws_data
        nw_id = parse(Int, n)
        nw_ref = ref[:nw][nw_id]
        tourList = genSubTour(nw_ref)
        nw_ref[:arc_tour] = []
        for tour in tourList
            hold = []
            for arc_tour in node2arcs(tour, nw_ref)
                append!(hold, arc_tour)
            end
            push!(nw_ref[:arc_tour], hold)
        end
    end
end

function ref_add_rdt!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    ref_add_new_lines!(ref, data)
    ref_add_damaged_lines!(ref, data)
    ref_add_vm_imbalance!(ref, data)
    ref_add_pq_imbalance!(ref, data)
    ref_add_subtour!(ref, data)
end
