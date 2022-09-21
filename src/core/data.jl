
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
# ["model", "name", "status", "active_phases", "qd", "vnom_kv", "source_id", "load_bus", "index", "conn", "pd"]

function correct_network_data!(data::Dict{String,Any})
    _PMs.check_connectivity(data)
    _PMs.correct_transformer_parameters!(data)
    _PMs.correct_voltage_angle_differences!(data)
    _PMs.correct_thermal_limits!(data)
    _PMs.correct_branch_directions!(data)
    _PMs.check_branch_loops(data)
    _PMs.correct_bus_types!(data)
    _PMs.correct_dcline_limits!(data)
    _PMs.correct_cost_functions!(data)
    _PMs.standardize_cost_terms!(data)
end
