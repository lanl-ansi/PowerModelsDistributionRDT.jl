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
    "baseMVA" ]        


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
    haskey(data,"scenarios") ? json2pm_scenarios!(data, pm_data, lookups) : println("No scenarios found in file")
    add_keys!(pm_data)
    return pm_data
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
    branch_data = Dict{String,Any}()
    lookups[:branch_names] = Dict{String,Any}()
    for (i, branch) in enumerate(data["lines"])
        # information about line 
        id = string(i)
        branch_data[id] = Dict{String,Any}()
        branch_data[id]["name"] = branch["id"]
        branch_data[id]["index"] = i
        lookups[:branch_names][branch["id"]] = i
        branch_data[id]["source_id"] = ["branch", i]
        !branch["is_new"] ? branch_data[id]["br_status"] = 1 : branch_data[id]["br_status"] = 0
      
        # transformer or not 
        branch_data[id]["transformer"] = branch["is_transformer"]
        branch_data[id]["tap"] = _PMs.MultiConductorVector([1.0 for i in 1:pm_data["conductors"]])
        branch_data[id]["shift"] = _PMs.MultiConductorVector([0.0 for i in 1:pm_data["conductors"]])

        # connections
        branch_data[id]["f_bus"] = lookups[:bus][branch["node1_id"]]
        branch_data[id]["t_bus"] = lookups[:bus][branch["node2_id"]]

        # line parameters 
        branch_data[id]["br_r"] = arrays_2_matrix(data["line_codes"][lookups[:branch][branch["line_code"]+1]], "rmatrix", convert(Float64,branch["length"])) 
        branch_data[id]["br_x"] = arrays_2_matrix(data["line_codes"][lookups[:branch][branch["line_code"]+1]], "xmatrix", convert(Float64,branch["length"]))
        branch_data[id]["g_to"] = make_zeros()
        branch_data[id]["g_fr"] = make_zeros()
        branch_data[id]["b_to"] = make_zeros()
        branch_data[id]["b_fr"] = make_zeros()

        # active phases
        branch_data[id]["active_phases"] = [i for i in 1:pm_data["conductors"] if branch["has_phase"][i]]

        # switch
        branch_data[id]["switch"] = branch["has_switch"]

        # capacity 
        branch_data[id]["rate_a"] = _PMs.MultiConductorVector([branch["capacity"] for i in 1:pm_data["conductors"]])
        branch_data[id]["rate_b"] = _PMs.MultiConductorVector([branch["capacity"] for i in 1:pm_data["conductors"]])
        branch_data[id]["rate_c"] = _PMs.MultiConductorVector([branch["capacity"] for i in 1:pm_data["conductors"]])

        # angle
        branch_data[id]["angmin"] = _PMs.MultiConductorVector([-.523599 for i in 1:pm_data["conductors"]])
        branch_data[id]["angmax"] = _PMs.MultiConductorVector([.523599 for i in 1:pm_data["conductors"]])

        # keys that maybe missing from some lines
        keys_ =["switch_cost", "harden_cost", "construction_cost", "num_poles"]
        for k in keys_
            if haskey(branch, k)
                branch_data[id][k] = branch[k]
            end
        end
    end
    pm_data["branch"] = branch_data
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
    gen_data = Dict{String,Any}()
    for (i, gen) in enumerate(data["generators"])
        id = string(i)
        # information about gen 
        gen_data[id] = Dict{String,Any}()
        gen_data[id]["name"] = gen["id"]
        gen_data[id]["index"] = i
        gen_data[id]["gen_bus"] = lookups[:bus][gen["node_id"]]
        gen_data[id]["source_id"] = ["gen", i]
        if !gen["is_new"]
            gen_data[id]["gen_status"] = 1
            if lookups[:type][gen["node_id"]] == 1
                lookups[:type][gen["node_id"]] = 2
            end
        else
            gen_data[id]["gen_status"] = 0
        end
        gen_data[id]["model"] = 2

        # active phases
        gen_data[id]["active_phases"] = [i for i in 1:pm_data["conductors"] if gen["has_phase"][i]]

        # gen parameters
        gen_data[id]["pg"] = _PMs.MultiConductorVector([0 for i in 1:pm_data["conductors"]])
        gen_data[id]["qg"] = _PMs.MultiConductorVector([0 for i in 1:pm_data["conductors"]])
        gen_data[id]["pmax"] = _PMs.MultiConductorVector(gen["max_real_phase"])
        gen_data[id]["pmin"] = _PMs.MultiConductorVector([0.0 for i in 1:pm_data["conductors"]])
        gen_data[id]["qmax"] = _PMs.MultiConductorVector(gen["max_reactive_phase"])
        gen_data[id]["qmin"] = -_PMs.MultiConductorVector(gen["max_reactive_phase"])
        # voltage ref
        gen_data[id]["vg"] = _PMs.MultiConductorVector(convert(Array{Float64}, data["buses"][lookups[:bus][gen["node_id"]]]["ref_voltage"]))

        #cost
        gen_data[id]["microgrid_cost"] = gen["microgrid_cost"]
    end
    pm_data["gen"] = gen_data
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
        bus_data[id]["vmin"] = _PMs.MultiConductorVector([bus["min_voltage"] for i in 1:pm_data["conductors"]])
        bus_data[id]["vmax"] = _PMs.MultiConductorVector([bus["max_voltage"] for i in 1:pm_data["conductors"]])
        bus_data[id]["vm"] = _PMs.MultiConductorVector(convert(Array{Float64}, bus["ref_voltage"]))
        bus_data[id]["va"] = _PMs.MultiConductorVector([0.0 for i in 1:pm_data["conductors"]])

        # area and zone
        bus_data[id]["zone"] = 1
        bus_data[id]["area"] = 1

        
        #coordinates
        bus_data[id]["x"] = bus["x"]
        bus_data[id]["y"] = bus["x"]
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
        load_data[id]["pd"] = _PMs.MultiConductorVector(load["max_real_phase"])
        load_data[id]["qd"] = _PMs.MultiConductorVector(load["max_reactive_phase"])

        # critical 
        load["is_critical"] ? load_data[id]["weight"] = 100 : load_data[id]["weight"] = 1
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
    value = _PMs.MultiConductorMatrix(value, conductors) * l
    # value[1] < 1e-9 ? value[1] = 1e-9 : nothing
    # value[5] < 1e-9 ? value[5] = 1e-9 : nothing
    # value[9] < 1e-9 ? value[9] = 1e-9 : nothing 
    return value
end

function make_zeros()
    conductors = 3
    value = zeros(Float64, 3, 3)
    return _PMs.MultiConductorMatrix(value, conductors)
end



