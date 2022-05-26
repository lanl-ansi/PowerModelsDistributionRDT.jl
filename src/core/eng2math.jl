
function calc_unique_virtual_bus_id_start(data::Dict{String,Any})
    virtual_bus_id = 0
    for i in keys(data["bus"])
        if parse(Int64,i) >= virtual_bus_id
            virtual_bus_id = parse(Int64,i)  + 1
        end
    end
    return virtual_bus_id
end

function calc_unique_switch_id_start(data::Dict{String,Any})
    switch_id = 0
    for i in keys(data["switch"])
        if parse(Int64,i) >= switch_id
            switch_id = parse(Int64,i)  + 1
        end
    end
    return switch_id
end

# Function looks for all lines where switches already exist and adds a switch
function transform_switch_inline!(data_math::Dict{String,Any}, data_eng::Dict{String,Any})
    switch_inline_data         = data_math["switch"]

    virtual_bus_id = calc_unique_virtual_bus_id_start(data_math)
    unique_switch_id = calc_unique_switch_id_start(data_math)

    for (i, branch) in data_math["branch"]
        if get(branch, "has_switch", false) == false
            # create the virtual bus
            bus_id = string(virtual_bus_id)

            # information about bus
            bus_info              = deepcopy(data_eng["bus"][string(branch["t_bus"])])
            bus_info["name"]      = virtual_bus_id
            bus_info["index"]     = virtual_bus_id
            bus_info["bus_i"]     = virtual_bus_id
            bus_info["source_id"] = ["bus", virtual_bus_id]

            # create the switch
            switch_info           = Dict{String,Any}()
            switch_id             = string(i)
            if haskey(data_math["switch"],i)
                switch_id = string(unique_switch_id)
                unique_switch_id = unique_switch_id + 1
            end

            # information about the switch
            switch_info["name"]          = switch_id
            switch_info["index"]         = parse(Int64,i)
            switch_info["status"]        = _PMD.ENABLED
            switch_info["f_bus"]         = virtual_bus_id
            switch_info["t_bus"]         = branch["t_bus"]
            switch_info["f_connections"] = deepcopy(branch["f_connections"])
            switch_info["t_connections"] = deepcopy(branch["t_connections"])
            switch_info["rate_a"]        = branch["rate_a"]
            switch_info["rate_b"]        = branch["rate_b"]
            switch_info["rate_c"]        = branch["rate_c"]
            switch_info["dispatchable"]  = _PMD.YES
            switch_info["state"]         = _PMD.CLOSED
            switch_info["switch_branch"] = i

            # final updates to everything
            branch["t_bus"] = virtual_bus_id
            switch_inline_data[switch_id] = switch_info
            data_math["bus"][bus_id]         = bus_info
            virtual_bus_id                   = virtual_bus_id + 1
        end
    end
end



# Function looks for all lines where switches cab be added (expanded) and creates ne_switches in these
# locations, which are always closed if unbuilt
function transform_switch_inline_ne!(data_math::Dict{String,Any}, data_eng::Dict{String,Any})
    switch_inline_ne_data         = Dict{String,Any}()

    virtual_bus_id = calc_unique_virtual_bus_id_start(data_math)

    for (i, branch) in data_math["branch"]
        if get(branch, "can_add_switch", false) == true || (get(branch, "switch_cost", 0.0) > 0.0 && get(branch, "has_switch", false) == false)
            # create the virtual bus
            bus_id = string(virtual_bus_id)

            # information about bus
            bus_info              = deepcopy(data_eng["bus"][string(branch["t_bus"])])
            bus_info["name"]      = virtual_bus_id
            bus_info["index"]     = virtual_bus_id
            bus_info["bus_i"]     = virtual_bus_id
            bus_info["source_id"] = ["bus", virtual_bus_id]

            # create the switch
            switch_info           = Dict{String,Any}()
            switch_id             = string(i)

            # information about the switch
            switch_info["name"]          = switch_id
            switch_info["index"]         = parse(Int64,i)
            switch_info["status"]        = _PMD.ENABLED
            switch_info["f_bus"]         = virtual_bus_id
            switch_info["t_bus"]         = branch["t_bus"]
            switch_info["f_connections"] = deepcopy(branch["f_connections"])
            switch_info["t_connections"] = deepcopy(branch["t_connections"])
            switch_info["rate_a"]        = branch["rate_a"]
            switch_info["rate_b"]        = branch["rate_b"]
            switch_info["rate_c"]        = branch["rate_c"]
            switch_info["dispatchable"]  = _PMD.YES
            switch_info["state"]         = _PMD.CLOSED
            switch_info["switch_cost"]   = branch["switch_cost"]
            switch_info["switch_branch"] = i

            # final updates to everything
            branch["t_bus"] = virtual_bus_id
            switch_inline_ne_data[switch_id] = switch_info
            data_math["bus"][bus_id]         = bus_info
            virtual_bus_id                   = virtual_bus_id + 1
        end
    end

    data_math["switch_inline_ne"] = switch_inline_ne_data
end