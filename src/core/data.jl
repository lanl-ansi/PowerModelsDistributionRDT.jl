
function add_load_weights!(data::Dict{String, Any})
    for (i, load) in data["load"]
        if !haskey(load, "weight")
            load["weight"] = 1
        end
    end
end

"Read in JSON file path and the nominal scenario; output a list of scenarios, each as a PowerModels data structure"
function gen_scenarios(filePath::String, data::Dict{String,Any})
    # only lines are affected by contingencies
    dataScen = JSON.Parser.parsefile(filePath);
    scenNo = length(dataScen["scenarios"]);
    scenList = 1:scenNo;
    scenData = [];
    for ω in scenList
        dataLocal = deepcopy(data);
        scenItem = dataScen["scenarios"][ω];
        branchNew = Dict{String,Any}();
        for iKey in keys(dataLocal["branch"])
            if !(iKey in scenItem["disable_lines"])
                branchNew[iKey] = dataLocal["branch"][iKey];
            end
        end
        dataLocal["branch"] = branchNew;
        push!(scenData,dataLocal);
    end
    return scenData;
end


function correct_network_data!(data::Dict{String,Any})
    _PMD.check_connectivity(data)
    _PM.correct_transformer_parameters!(data)
    _PM.correct_voltage_angle_differences!(data)
    _PM.correct_thermal_limits!(data)
    _PMD.correct_branch_directions!(data)
    _PMD.check_branch_loops(data)
    _PMD.correct_bus_types!(data)
    _PM.correct_dcline_limits!(data)
    _PMD.correct_cost_functions!(data)
    _PMD.standardize_cost_terms!(data)
end


function calc_total_real_load(loads::Dict)
    pd = zeros(Float64,length(first(loads)[2]["pd"]))

    for (i, load) in loads
        pd = pd + load["pd"]
    end

    return pd
end


function calc_total_reactive_load(loads::Dict)
    qd = zeros(Float64,length(first(loads)[2]["qd"]))

    for (i, load) in loads
        qd = qd + load["qd"]
    end

    return qd
end


function calc_branch_current_max(branch::Dict{String,<:Any}, bus::Dict{String,<:Any}, total_real_load::Vector{Float64}, total_reactive_load::Vector{Float64})::Vector{Float64}
    connections = [findfirst(isequal(c), bus["terminals"]) for c in (branch["f_bus"] == bus["index"] ? branch["f_connections"] : branch["t_connections"])]
    max_current = sqrt.(total_real_load.^2 + total_reactive_load.^2)./bus["vmin"][connections]

    if haskey(branch, "c_rating_a")
        return min.(branch["c_rating_a"],max_current)
    else
        return max_current
    end
end


function calc_branch_power_max(branch::Dict{String,<:Any}, bus::Dict{String,<:Any}, total_real_load::Vector{Float64}, total_reactive_load::Vector{Float64})::Vector{Float64}
    connections = [findfirst(isequal(c), bus["terminals"]) for c in (branch["f_bus"] == bus["index"] ? branch["f_connections"] : branch["t_connections"])]
    max_power = sqrt.(total_real_load.^2 + total_reactive_load.^2)

    if haskey(branch, "rating_a")
        return min.(branch["rating_a"],max_power)
    else
        return max_power
    end
end
