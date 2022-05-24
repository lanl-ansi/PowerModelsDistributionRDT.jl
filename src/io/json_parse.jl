import JSON
using LinearAlgebra

pms_keys = [
    "bus",
    "source_type",
    "name",
    "dcline",
    "source_version",
    "branch",
    "gen",
    "storage",
    "switch",
    "conductors",
    "shunt",
    "load",
    "per_unit",
    "baseMVA",
    "data_model",
    "transformer" ]

global_keys = Set{String}([
    "phase_variation",
    "total_load_met",
    "critical_load_met",
    "chance_constraint",
    "scenarios",
    "data_model",
    ])

function parse_json(io::IOStream; validate=false)
    json_data = JSON.parse(io; dicttype=Dict, inttype=Int64)
    return pm_data = json_to_powermodels(json_data)
end

function json_to_powermodels(data::Dict{String,Any})
    pm_data = Dict{String,Any}()
    pm_data["conductors"] = 3
    lookups = create_lookups(data)
    pm_data["phase_variation"] = convert(Float64,data["phase_variation"])
    pm_data["total_load_met"] = convert(Float64,data["total_load_met"])
    pm_data["critical_load_met"] = convert(Float64,data["critical_load_met"])
    pm_data["chance_constraint"] = convert(Float64,data["chance_constraint"])
    json2pm_branch!(data, pm_data, lookups)
    json2pm_gen!(data, pm_data, lookups)
    json2pm_bus!(data, pm_data, lookups)
    json2pm_load!(data, pm_data, lookups)
    pm_data["settings"] = get(pm_data,"settings", Dict{String,Any}("sbase_default" => get(data, "baseMVA", 1e6)))
    add_keys!(pm_data)
    correct_network_data!(pm_data)
    haskey(data,"scenarios") ? json2pm_scenarios!(data, pm_data, lookups) : println("No scenarios found in file")
    haskey(pm_data, "scenarios") ?  mn_data = _PMs.replicate(pm_data, length(keys(pm_data["scenarios"]))+1, global_keys=global_keys) : mn_data = pm_data

    mn_data["nw"]["0"] = mn_data["nw"][string(length(keys(pm_data["scenarios"]))+1)]
    delete!(mn_data["nw"], string(length(keys(pm_data["scenarios"]))+1))

    println("Need to actually use the damage information in the scenarios to define stuff!!!!!")

    # not sure why this isn't getting replicated
    for n in keys(mn_data["nw"])
        mn_data["nw"][n]["per_unit"] = pm_data["per_unit"]
    end

    return mn_data
end

function create_lookups(data::Dict{String,Any})
    lookups = Dict{Symbol,Any}()
    lookups[:bus] = Dict{String,Any}()
    lookups[:type] = Dict{String,Any}()
    for (i, bus) in enumerate(data["buses"])
        lookups[:bus][bus["id"]] = i
        occursin("source", bus["id"]) ? lookups[:type][bus["id"]] = 3 : lookups[:type][bus["id"]] = 1
    end
    lookups[:branch] = Dict{Int,Any}()
    for (i, branch) in enumerate(data["line_codes"])
        lookups[:branch][branch["line_code"]+1] = i
    end
    return lookups
end

function add_keys!(data::Dict{String,Any})
    for key in pms_keys
        !haskey(data, key) ? data[key] = Dict{String,Any}() : nothing
        "per_unit" == key ? data[key] = true : nothing
        #"data_model" == key ? data[key] = _PMD.ENGINEERING : nothing
        "data_model" == key ? data[key] = _PMD.MATHEMATICAL : nothing
        "transformer" == key ? data[key] = Dict{String, Any}() : nothing
    end
end

# The lines are defined with the following parameters
# id                        A string that uniquely identifies the line
# node1_id                  The id of one node (bus) the line is connected to
# node2_id                  The id of one node (bus) the line is connected to
# has_phase                 An array of whether or not a line carries a phase
# capacity                  per unit capacity (rating) of the line
# length                    per unit length of the line
# num_phases                The number of phases the line carries
# is_transformer            Whether or the line is a transformer
# line_code                 Id of the line code for the line
# construction_cost         The construction cost for building this line
# harden_cost               The hardening cost for the line
# switch_cost               The cost for building a switch at the line
# is_new                    Boolean for indicating if the line needs to be built
# can_harden                Boolean for indicating if the line can be hardened
# can_add_switch            Boolean for indicating if a switch can be added
# has_switch                Boolean for if a switch is already on the line
# '
# 'The line codes are defined with the following parameters
# line code                    A string that uniquely identifies the line code
# num_phases                   Number of phases for the line code
# rmatrix                      An array of the resistance terms for the line
#                              (per length, per unit).
# xmatrix                      An array of the reactance terms for the line
#                              (per length, per unit).
#
function json2pm_branch!(data::Dict{String,Any}, pm_data::Dict{String,Any}, lookups::Dict{Symbol,Any})
    branch_data            = Dict{String,Any}()
    branch_ne_data         = Dict{String,Any}()
    lookups[:branch_names] = Dict{String,Any}()
    for (i, branch) in enumerate(data["lines"])
        # information about line
        id = string(i)
        info = Dict{String,Any}()

        if branch["is_new"]
            branch_ne_data[id]  = info
        else
            branch_data[id] = info
        end

        info["name"] = branch["id"]
        info["index"] = i
        lookups[:branch_names][branch["id"]] = i
        info["source_id"] = ["branch", i]
        info["br_status"] = 1

        # transformer or not
        info["transformer"] = branch["is_transformer"]
        info["tap"] = Array{Float64,1}([1.0 for i in 1:pm_data["conductors"]])
        info["shift"] = Array{Float64,1}([0.0 for i in 1:pm_data["conductors"]])

        # connections
        info["f_bus"] = lookups[:bus][branch["node1_id"]]
        info["t_bus"] = lookups[:bus][branch["node2_id"]]

        # line parameters
        info["br_r"] = arrays_2_matrix(data["line_codes"][lookups[:branch][branch["line_code"]+1]], "rmatrix", convert(Float64,branch["length"]))
        info["br_x"] = arrays_2_matrix(data["line_codes"][lookups[:branch][branch["line_code"]+1]], "xmatrix", convert(Float64,branch["length"]))
        info["g_to"] = make_zeros()
        info["g_fr"] = make_zeros()
        info["b_to"] = make_zeros()
        info["b_fr"] = make_zeros()

        # active phases
        info["active_phases"] = [i for i in 1:pm_data["conductors"] if branch["has_phase"][i]]

        # switch
        info["switch"] = branch["has_switch"]

        # capacity
        info["rate_a"] = Array{Float64,1}([branch["capacity"] for i in 1:pm_data["conductors"]])
        info["rate_b"] = Array{Float64,1}([branch["capacity"] for i in 1:pm_data["conductors"]])
        info["rate_c"] = Array{Float64,1}([branch["capacity"] for i in 1:pm_data["conductors"]])

        # angle
        info["angmin"] = Array{Float64,1}([-.523599 for i in 1:pm_data["conductors"]])
        info["angmax"] = Array{Float64,1}([.523599 for i in 1:pm_data["conductors"]])

        info["f_connections"] = collect(1:pm_data["conductors"])
        info["t_connections"] = collect(1:pm_data["conductors"])

        # keys that maybe missing from some lines
        keys_ =["switch_cost", "harden_cost", "construction_cost", "num_poles", "can_harden", "can_add_switch", "is_new", "has_switch"]
        for k in keys_
            if haskey(branch, k)
                info[k] = branch[k]
            end
        end

        if !haskey(info, "can_harden") && info["is_new"] == false
            info["can_harden"] = haskey(info, "harden_cost") ? true : false
        end

    end
    pm_data["branch"]    = branch_data
    pm_data["branch_ne"] = branch_ne_data
end



# The generators (power consumption points) are defined with the following parameters
# id                        A string that uniquely identifies the generator
# node_id                   The id of the node (bus) the generator is connected to
# has_phase                 An array of whether or not a generator has a phase
# max_real_phase            An array of the maximum active generation in each phase
# max_reactive_phase        An array of the maximum reactive generation in each phase
# microgrid_cost            The cost for building the generator
# is_new                    Boolean to indicate whether or not this generator has to
#                           be built.
function json2pm_gen!(data::Dict{String,Any}, pm_data::Dict{String,Any}, lookups)
    gen_data    = Dict{String,Any}()
    gen_ne_data = Dict{String,Any}()

    for (i, gen) in enumerate(data["generators"])
        id = string(i)
        info = Dict{String,Any}()

        if gen["is_new"]
            gen_ne_data[id] = info
        else
            gen_data[id] = info
        end

        # information about gen
        info["name"] = gen["id"]
        info["index"] = i
        info["gen_bus"] = lookups[:bus][gen["node_id"]]
        info["source_id"] = ["gen", i]
        info["gen_status"] = 1
        info["model"] = 2

        # active phases
        info["active_phases"] = [i for i in 1:pm_data["conductors"] if gen["has_phase"][i]]

        # gen parameters
        info["pg"] = Array{Float64,1}([0 for i in 1:pm_data["conductors"]])
        info["qg"] = Array{Float64,1}([0 for i in 1:pm_data["conductors"]])
        info["pmax"] = Array{Float64,1}(gen["max_real_phase"])
        info["pmin"] = Array{Float64,1}([0.0 for i in 1:pm_data["conductors"]])
        info["qmax"] = Array{Float64,1}(gen["max_reactive_phase"])
        info["qmin"] = -Array{Float64,1}(gen["max_reactive_phase"])
        # voltage ref
        info["vg"] = Array{Float64,1}(convert(Array{Float64}, data["buses"][lookups[:bus][gen["node_id"]]]["ref_voltage"]))

        info["microgrid_cost"] = gen["microgrid_cost"]
        # applys generic cost
        info["ncost"] = 3
        info["cost"] = [0.0, 1.0, 0.0]

        info["connections"] = collect(1:pm_data["conductors"])
    end
    pm_data["gen"]    = gen_data
    pm_data["gen_ne"] = gen_ne_data
end



# The buses (nodes) are defined with the following parameters
# id                    A string that uniquely identifies the bus
# min_voltage           The minimum voltage level at the bus in p.u.
# max_voltage           The maximum voltage level at the bus in p.u.
# ref_voltage           The reference voltage level at the bus in p.u.
# has_phase             An array of whether or not a bus has a phase
#
function json2pm_bus!(data::Dict{String,Any}, pm_data::Dict{String,Any}, lookups)
    bus_data = Dict{String,Any}()
    for (i, bus) in enumerate(data["buses"])
        id = string(i)
        # information about bus
        bus_data[id] = Dict{String,Any}()
        bus_data[id]["name"] = bus["id"]
        bus_data[id]["index"] = i
        bus_data[id]["bus_i"] = lookups[:bus][bus["id"]]
        bus_data[id]["source_id"] = ["bus", i]
        bus_data[id]["bus_type"] = lookups[:type][bus["id"]]

        # active phases
        bus_data[id]["active_phases"] = [i for i in 1:pm_data["conductors"] if bus["has_phase"][i]]

        # bus parameters
        bus_data[id]["vmin"] = Array{Float64,1}([bus["min_voltage"] for i in 1:pm_data["conductors"]])
        bus_data[id]["vmax"] = Array{Float64,1}([bus["max_voltage"] for i in 1:pm_data["conductors"]])
        bus_data[id]["vm"] = Array{Float64,1}(convert(Array{Float64}, bus["ref_voltage"]))
        bus_data[id]["va"] = Array{Float64,1}([0.0 for i in 1:pm_data["conductors"]])

        # area and zone
        bus_data[id]["zone"] = 1
        bus_data[id]["area"] = 1


        #coordinates
        bus_data[id]["x"] = bus["x"]
        bus_data[id]["y"] = bus["y"]

        # Things needed for an engineerng model
        bus_data[id]["terminals"] = collect(1:pm_data["conductors"])
        bus_data[id]["grounded"]  = fill(false, pm_data["conductors"])
        bus_data[id]["status"]    = 1
        bus_data[id]["rg"]        = [0.0, 0.0, 0.0]
        bus_data[id]["xg"]        = [0.0, 0.0, 0.0]
    end
    pm_data["bus"] = bus_data
end



# The loads (power consumption points) are defined with the following parameters
# id                        A string that uniquely identifies the load
# node_id                   The id of the node (bus) the load is connected to
# has_phase                 An array of whether or not a load has a phase
# max_real_phase            An array of the desired active load in each phase
# max_reactive_phase        An array of the desired reactive load in each phase
# is_critical               Boolean for whether or not the load is critical
#
function json2pm_load!(data::Dict{String,Any}, pm_data::Dict{String,Any}, lookups)
    load_data = Dict{String,Any}()
    for (i, load) in enumerate(data["loads"])
        id = string(i)
        # information about load
        load_data[id] = Dict{String,Any}()
        load_data[id]["name"] = load["id"]
        load_data[id]["index"] = i
        load_data[id]["load_bus"] = lookups[:bus][load["node_id"]]
        load_data[id]["source_id"] = ["bus", i]
        load_data[id]["status"] = 1

        # active phases
        load_data[id]["active_phases"] = [i for i in 1:pm_data["conductors"] if load["has_phase"][i]]

        # bus parameters
        load_data[id]["pd"] = Array{Float64,1}(load["max_real_phase"])
        load_data[id]["qd"] = Array{Float64,1}(load["max_reactive_phase"])

        # critical
        load["is_critical"] ? load_data[id]["weight"] = 100 : load_data[id]["weight"] = 1

        load_data[id]["model"] = _PMD.POWER
        load_data[id]["connections"] = [i for i in 1:pm_data["conductors"] if load["has_phase"][i]]
        push!(load_data[id]["connections"], 4) # 4 = neutral
        load_data[id]["configuration"] = _PMD.WYE
        load_data[id]["bus"] = string(lookups[:bus][load["node_id"]])
        load_data[id]["pd_nom"] = load_data[id]["pd"]
        load_data[id]["qd_nom"] = load_data[id]["qd"]
        load_data[id]["vm_nom"] = pm_data["bus"][string(lookups[:bus][load["node_id"]])]["vm"]
    end
    pm_data["load"] = load_data
end



# The scenarios are defined with the following parameters
# id                             A string that uniquely identifies the scenario
# hardened_disabled_lines        An array of power lines that are damaged after being
#                                hardened
# disabled_communication_lines   An array of communication lines that are damaged
# disabled_lines                 An array of power lines that are damaged
#
function json2pm_scenarios!(data::Dict{String,Any}, pm_data::Dict{String,Any}, lookups)
    scenario_data = Dict{String,Any}()
    for (i, s) in enumerate(data["scenarios"])
        id = string(i)
        # information about senarios
        scenario_data[id] = Dict{String,Any}()
        scenario_data[id]["name"] = s["id"]
        haskey(s, "disable_lines") ? scenario_data[id]["disabled_lines"] = [string(lookups[:branch_names][i]) for i in s["disable_lines"]] : nothing
        haskey(s, "hardened_disabled_lines") ? scenario_data[id]["hardened_disabled_lines"] = [string(lookups[:branch_names][i]) for i in s["hardened_disabled_lines"]] : nothing
        haskey(s, "disabled_communication_lines") ? scenario_data[id]["disabled_communication_lines"] = [string(lookups[:branch_names][i]) for i in s["disabled_communication_lines"]] : nothing
    end
    pm_data["scenarios"] = scenario_data
end



function json2pm_shunt!(data::Dict{String,Any}, pm_data::Dict{String,Any}, lookups)
    shunt_data = Dict{String,Any}()
    pm_data["shunt"] = shunt_data
end

function json2pm_storage!(data::Dict{String,Any}, pm_data::Dict{String,Any}, lookups)
    storage_data = Dict{String,Any}()
    pm_data["storage"] = storage_data
end

function arrays_2_matrix(m::Dict{String,Any}, string::String, l::Float64)
    conductors = 3
    value = zeros(Float64, 0, conductors)
    for r in m[string]
        value = [value;transpose(convert(Array{Float64,1},r))]
    end
    value[1] < 1e-9 ? value[1] = 1e-9 : nothing
    value[5] < 1e-9 ? value[5] = 1e-9 : nothing
    value[9] < 1e-9 ? value[9] = 1e-9 : nothing
    return value
end

function make_zeros()
    conductors = 3
    value = zeros(Float64, 3, 3)
    return value
end
