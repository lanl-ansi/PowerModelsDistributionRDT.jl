
function add_load_weights!(data::Dict{String,Any})
    for (i, load) in data["load"]
        if !haskey(load, "weight")
            load["weight"] = 1
        end
    end
end

"Read in JSON file path and the nominal scenario; output a list of scenarios, each as a PowerModels data structure"
function gen_scenarios(filePath::String, data::Dict{String,Any})
    # only lines are affected by contingencies
    dataScen = JSON.Parser.parsefile(filePath)
    scenNo = length(dataScen["scenarios"])
    scenList = 1:scenNo
    scenData = []
    for ω in scenList
        dataLocal = deepcopy(data)
        scenItem = dataScen["scenarios"][ω]
        branchNew = Dict{String,Any}()
        for iKey in keys(dataLocal["branch"])
            if !(iKey in scenItem["disable_lines"])
                branchNew[iKey] = dataLocal["branch"][iKey]
            end
        end
        dataLocal["branch"] = branchNew
        push!(scenData, dataLocal)
    end
    return scenData
end


function correct_network_data!(data::Dict{String,Any})
    _PMD.check_connectivity(data)

    _PMD.correct_branch_directions!(data)
    _PMD.check_branch_loops(data)

    _PMD.correct_bus_types!(data)

    _PMD.propagate_network_topology!(data)

    _PMD.correct_mc_voltage_angle_differences!(data)
    _PMD.correct_mc_thermal_limits!(data)

    _PMD.correct_cost_functions!(data)
    _PMD.standardize_cost_terms!(data)
end


function replicate(sn_data::Dict{String,<:Any}, count::Int; global_keys::Set{String}=Set{String}())
    return _IM.replicate(sn_data, count, union(global_keys, _PMD._pmd_global_keys))
end


function calc_total_real_load(loads::Dict)
    pd = zeros(Float64, length(first(loads)[2]["pd"]))

    for (i, load) in loads
        pd = pd + load["pd"]
    end

    return pd
end


function calc_total_reactive_load(loads::Dict)
    qd = zeros(Float64, length(first(loads)[2]["qd"]))

    for (i, load) in loads
        qd = qd + load["qd"]
    end

    return qd
end


function calc_branch_current_max(branch::Dict{String,<:Any}, bus::Dict{String,<:Any}, total_real_load::Vector{Float64}, total_reactive_load::Vector{Float64})::Vector{Float64}
    connections = [findfirst(isequal(c), bus["terminals"]) for c in (branch["f_bus"] == bus["index"] ? branch["f_connections"] : branch["t_connections"])]
    max_current = sqrt.(total_real_load .^ 2 + total_reactive_load .^ 2) ./ bus["vmin"][connections]

    if haskey(branch, "c_rating_a")
        return min.(branch["c_rating_a"], max_current)
    else
        return max_current
    end
end


function calc_branch_power_max(branch::Dict{String,<:Any}, bus::Dict{String,<:Any}, total_real_load::Vector{Float64}, total_reactive_load::Vector{Float64})::Vector{Float64}
    connections = [findfirst(isequal(c), bus["terminals"]) for c in (branch["f_bus"] == bus["index"] ? branch["f_connections"] : branch["t_connections"])]
    max_power = sqrt.(total_real_load .^ 2 + total_reactive_load .^ 2)

    if haskey(branch, "rating_a")
        return min.(branch["rating_a"], max_power)
    else
        return max_power
    end
end


"TODO: This is a (likely) incorrect hack that should ultimate live in PowerModelsDistribution?"
function calc_theta_delta_bounds(data::Dict{String,<:Any})
    bus_count = length(data["bus"])
    branches = [branch for branch in values(data["branch"])]
    if haskey(data, "ne_branch")
        append!(branches, values(data["ne_branch"]))
    end

    angle_min = Real[]
    angle_max = Real[]

    conductors = 1
    if haskey(data, "conductors")
        conductors = data["conductors"]
    end
    conductor_ids = 1:conductors

    for c in conductor_ids
        angle_mins = [branch["angmin"][c] for branch in branches]
        angle_maxs = [branch["angmax"][c] for branch in branches]

        sort!(angle_mins)
        sort!(angle_maxs, rev=true)

        if length(angle_mins) > 1
            # note that, this can occur when dclines are present
            angle_count = min(bus_count - 1, length(branches))

            angle_min_val = sum(angle_mins[1:angle_count])
            angle_max_val = sum(angle_maxs[1:angle_count])
        else
            angle_min_val = angle_mins[1]
            angle_max_val = angle_maxs[1]
        end

        push!(angle_min, angle_min_val)
        push!(angle_max, angle_max_val)
    end

    if haskey(data, "conductors")
        return angle_min, angle_max
    else
        return angle_min[1], angle_max[1]
    end
end


"""
    calc_connected_components(data::Dict{String,<:Any}; edges::Union{Missing, Vector{<:String}}=missing, type::Union{Missing,String}=missing, check_enabled::Bool=true)::Set

computes the connected components of the network graph
returns a set of sets of bus ids, each set is a connected component

This function differs from the powermodelsdistribution implementation of this functionality in that it accounts for edges being damaged (with the option to harden them) - which functionally makes
    these edges act like switches in the context of computing load blocks

In computing connected components, the default is to ignore expansion edges (branch_ne, switch_inline_ne)

"""
function calc_connected_components(data::Dict{String,<:Any}; edges::Union{Missing,Vector{String}}=missing, type::Union{Missing,String}=missing, check_enabled::Bool=true)::Set{Set}
    pmd_data = _PMD.get_pmd_data(data)

    if ismultinetwork(pmd_data)
        error("multinetwork data is not yet supported, recommend to use on each subnetwork independently")
    end

    if get(pmd_data, "data_model", _PMD.MATHEMATICAL) == _PMD.ENGINEERING
        return _calc_connected_components_eng(pmd_data; edges=ismissing(edges) ? _PMD._eng_edge_elements : edges, type=type, check_enabled=check_enabled)
    elseif get(pmd_data, "data_model", _PMD.MATHEMATICAL) == _PMD.MATHEMATICAL
        return _calc_connected_components_math(pmd_data; edges=ismissing(edges) ? _PMD._math_edge_elements : edges, type=type, check_enabled=check_enabled)
    else
        error("data_model `$(get(pmd_data, "data_model", MATHEMATICAL))` is unrecongized")
    end
end


"""
computes the connected components of the network graph
returns a set of sets of bus ids, each set is a connected component
"""
function _calc_connected_components_eng(data; edges::Vector{<:String}=_eng_edge_elements, type::Union{Missing,String}=missing, check_enabled::Bool=true)::Set{Set{String}}
    @assert get(data, "data_model", _PMD.MATHEMATICAL) == _PMD.ENGINEERING

    active_bus = Dict{String,Dict{String,Any}}(x for x in data["bus"] if x.second["status"] == ENABLED || !check_enabled)
    active_bus_ids = Set{String}([i for (i, bus) in active_bus])

    neighbors = Dict{String,Vector{String}}(i => [] for i in active_bus_ids)
    for edge_type in edges
        for (id, edge_obj) in get(data, edge_type, Dict{Any,Dict{String,Any}}())
            if edge_obj["status"] == ENABLED || !check_enabled
                if edge_type == "transformer" && haskey(edge_obj, "bus")
                    for f_bus in edge_obj["bus"]
                        for t_bus in edge_obj["bus"]
                            if f_bus != t_bus
                                push!(neighbors[f_bus], t_bus)
                                push!(neighbors[t_bus], f_bus)
                            end
                        end
                    end
                elseif edge_type == "branch"
                    if type != "load_blocks" || get(edge_obj, "is_damaged", false) == false
                        push!(neighbors[edge_obj["f_bus"]], edge_obj["t_bus"])
                        push!(neighbors[edge_obj["t_bus"]], edge_obj["f_bus"])
                    end
                else
                    if (edge_type == "switch" || edge_type == "switch_inline_ne") && !ismissing(type)
                        if type == "load_blocks"
                            if edge_obj["dispatchable"] == NO && edge_obj["state"] == CLOSED
                                push!(neighbors[edge_obj["f_bus"]], edge_obj["t_bus"])
                                push!(neighbors[edge_obj["t_bus"]], edge_obj["f_bus"])
                            end
                        elseif type == "blocks"
                            if edge_obj["state"] == CLOSED
                                push!(neighbors[edge_obj["f_bus"]], edge_obj["t_bus"])
                                push!(neighbors[edge_obj["t_bus"]], edge_obj["f_bus"])
                            end
                        end
                    else
                        push!(neighbors[edge_obj["f_bus"]], edge_obj["t_bus"])
                        push!(neighbors[edge_obj["t_bus"]], edge_obj["f_bus"])
                    end
                end
            end
        end
    end

    component_lookup = Dict(i => Set{String}([i]) for i in active_bus_ids)
    touched = Set{String}()

    for i in active_bus_ids
        if !(i in touched)
            _PMD._cc_dfs(i, neighbors, component_lookup, touched)
        end
    end

    return Set{Set{String}}(values(component_lookup))
end


"""
computes the connected components of the network graph
returns a set of sets of bus ids, each set is a connected component
"""
function _calc_connected_components_math(data::Dict{String,<:Any}; edges::Vector{<:String}=_math_edge_elements, type::Union{Missing,String}=missing, check_enabled::Bool=true)::Set{Set{Int}}
    @assert get(data, "data_model", _PMD.MATHEMATICAL) == _PMD.MATHEMATICAL

    active_bus = Dict{String,Dict{String,Any}}(x for x in data["bus"] if x.second[_PMD.pmd_math_component_status["bus"]] != _PMD.pmd_math_component_status_inactive["bus"] || !check_enabled)
    active_bus_ids = Set{Int}([parse(Int, i) for (i, bus) in active_bus])

    neighbors = Dict{Int,Vector{Int}}(i => [] for i in active_bus_ids)
    for edge_type in edges
        for (id, edge_obj) in get(data, edge_type, Dict{Any,Dict{String,Any}}())
            if edge_obj[_PMD.pmd_math_component_status[edge_type]] != _PMD.pmd_math_component_status_inactive[edge_type] || !check_enabled
                if (edge_type == "switch" || edge_type == "switch_inline_ne") && !ismissing(type)
                    if type == "load_blocks"
                        if edge_obj["dispatchable"] != 1 && edge_obj["state"] == 1
                            push!(neighbors[edge_obj["f_bus"]], edge_obj["t_bus"])
                            push!(neighbors[edge_obj["t_bus"]], edge_obj["f_bus"])
                        end
                    elseif type == "blocks"
                        if edge_obj["state"] != 0
                            push!(neighbors[edge_obj["f_bus"]], edge_obj["t_bus"])
                            push!(neighbors[edge_obj["t_bus"]], edge_obj["f_bus"])
                        end
                    end
                elseif edge_type == "branch"
                    if type != "load_blocks" || get(edge_obj, "is_damaged", false) == false
                        push!(neighbors[edge_obj["f_bus"]], edge_obj["t_bus"])
                        push!(neighbors[edge_obj["t_bus"]], edge_obj["f_bus"])
                    end
                else
                    push!(neighbors[edge_obj["f_bus"]], edge_obj["t_bus"])
                    push!(neighbors[edge_obj["t_bus"]], edge_obj["f_bus"])
                end
            end
        end
    end

    component_lookup = Dict(i => Set{Int}([i]) for i in active_bus_ids)
    touched = Set{Int}()

    for i in active_bus_ids
        if !(i in touched)
            _PMD._cc_dfs(i, neighbors, component_lookup, touched)
        end
    end

    return Set{Set{Int}}(values(component_lookup))
end
