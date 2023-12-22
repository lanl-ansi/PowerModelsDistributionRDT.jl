
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


function ref_add_damaged_tag!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    for (nw, nw_ref) in ref[:it][_PMD.pmd_it_sym][:nw]
        for i in nw_ref[:damaged_branch]
            nw_ref[:branch][i]["is_damaged"] = true
        end
    end
end

function ref_add_global_constants!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    for (nw, nw_ref) in ref[:it][_PMD.pmd_it_sym][:nw]
        nw_ref[:total_real_load] = calc_total_real_load(nw_ref[:load])
        nw_ref[:total_reactive_load] = calc_total_reactive_load(nw_ref[:load])
        nw_ref[:off_angmin], nw_ref[:off_angmax] = calc_theta_delta_bounds(data["nw"][string(nw)])
    end
end

"""
    _ref_add_load_blocks!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})

Ref extension to add load blocks to ref for a single element of a multi-network
    The key difference with the implementation of PowerModelsONM is that the observation that some edges can be damaged potentially increases the number of load blocks
"""
function _ref_add_load_blocks!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    ref[:blocks] = Dict{Int,Set}(i => block.second for (i,block) in enumerate(sort([sum(map(x->SHA.sha1(string(ref[:bus][x]["name"])), collect(b)))=>b for b in calc_connected_components(data; type="load_blocks", check_enabled=true)]; by=x->x.first)))
    ref[:bus_block_map] = Dict{Int,Int}(bus => b for (b,block) in ref[:blocks] for bus in block)
#    ref[:block_branches] = Dict{Int,Set}(b => Set{Int}() for (b,_) in ref[:blocks])
#    ref[:block_loads] = Dict{Int,Set}(i => Set{Int}() for (i,_) in ref[:blocks])
#    ref[:block_weights] = Dict{Int,Real}(i => 1.0 for (i,_) in ref[:blocks])
#    ref[:block_shunts] = Dict{Int,Set{Int}}(i => Set{Int}() for (i,_) in ref[:blocks])
#    ref[:block_gens] = Dict{Int,Set{Int}}(i => Set{Int}() for (i,_) in ref[:blocks])
#    ref[:block_storages] = Dict{Int,Set{Int}}(i => Set{Int}() for (i,_) in ref[:blocks])
#    ref[:microgrid_blocks] = Dict{Int,String}()
#    ref[:substation_blocks] = Vector{Int}()
#    ref[:bus_inverters] = Dict{Int,Set{Tuple{Symbol,Int}}}(i => Set{Tuple{Symbol,Int}}() for (i,_) in ref[:bus])
#    ref[:block_inverters] = Dict{Int,Set{Tuple{Symbol,Int}}}(b => Set{Tuple{Symbol,Int}}() for (b,_) in ref[:blocks])
#    ref[:dispatchable_loads] = Dict{Int, Dict}(i => load for (i,load) in ref[:load] if Int(load["dispatchable"]) == Int(_PMD.YES))
#    ref[:nondispatchable_loads] = Dict{Int, Dict}(i => load for (i,load) in ref[:load] if Int(load["dispatchable"]) == Int(_PMD.NO))
#    ref[:block_dispatchable_loads] = Dict{Int,Set}(i => Set{Int}() for (i,_) in ref[:blocks])

#    for (b,bus) in ref[:bus]
#        if !isempty(get(bus, "microgrid_id", ""))
#            ref[:block_weights][ref[:bus_block_map][b]] = 10.0
#            ref[:microgrid_blocks][ref[:bus_block_map][b]] = bus["microgrid_id"]
#        end
#    end

#    for (br,branch) in ref[:branch]
#        push!(ref[:block_branches][ref[:bus_block_map][branch["f_bus"]]], br)
#    end
#    ref[:block_line_losses] = Dict{Int,Float64}(i => sum(Float64[LinearAlgebra.norm(ref[:branch][br]["br_r"].+1im*ref[:branch][br]["br_x"]) for br in branches if ref[:branch][br][_PMD.pmd_math_component_status["branch"]] != _PMD.pmd_math_component_status_inactive["branch"]]) for (i,branches) in ref[:block_branches])

#    for (l,load) in ref[:load]
#        push!(ref[:block_loads][ref[:bus_block_map][load["load_bus"]]], l)
#        ref[:block_weights][ref[:bus_block_map][load["load_bus"]]] += 1e-2 * get(load, "priority", 1)
#        Int(load["dispatchable"]) == Int(_PMD.YES) && push!(ref[:block_dispatchable_loads][ref[:bus_block_map][load["load_bus"]]], l)
#    end
#    ref[:load_block_map] = Dict{Int,Int}(load => b for (b,block_loads) in ref[:block_loads] for load in block_loads)

#    for (s,shunt) in ref[:shunt]
#        push!(ref[:block_shunts][ref[:bus_block_map][shunt["shunt_bus"]]], s)
#    end
#    ref[:shunt_block_map] = Dict{Int,Int}(shunt => b for (b,block_shunts) in ref[:block_shunts] for shunt in block_shunts)

#    for (g,gen) in ref[:gen]
#        push!(ref[:block_gens][ref[:bus_block_map][gen["gen_bus"]]], g)
#        startswith(gen["source_id"], "voltage_source") && push!(ref[:substation_blocks], ref[:bus_block_map][gen["gen_bus"]])
#        push!(ref[:bus_inverters][gen["gen_bus"]], (:gen, g))
#        push!(ref[:block_inverters][ref[:bus_block_map][gen["gen_bus"]]], (:gen, g))
#    end
#    ref[:gen_block_map] = Dict{Int,Int}(gen => b for (b,block_gens) in ref[:block_gens] for gen in block_gens)

#    for (s,strg) in ref[:storage]
#        push!(ref[:block_storages][ref[:bus_block_map][strg["storage_bus"]]], s)
#        push!(ref[:bus_inverters][strg["storage_bus"]], (:storage, s))
#        push!(ref[:block_inverters][ref[:bus_block_map][strg["storage_bus"]]], (:storage, s))
#    end
#    ref[:storage_block_map] = Dict{Int,Int}(strg => b for (b,block_storages) in ref[:block_storages] for strg in block_storages)

#    for (i,_) in ref[:blocks]
#        if isempty(ref[:block_loads][i]) && isempty(ref[:block_shunts][i]) && isempty(ref[:block_gens][i]) && isempty(ref[:block_storages][i])
#            ref[:block_weights][i] = 0.0
#        end
#    end

#    ref[:block_graph] = Graphs.SimpleGraph(length(ref[:blocks]))
#    ref[:block_graph_edge_map] = Dict{Graphs.Edge,Int}()
#    ref[:block_switches] = Dict{Int,Set{Int}}(b => Set{Int}() for (b,_) in ref[:blocks])

#    for (s,switch) in ref[:switch]
#        f_block = ref[:bus_block_map][switch["f_bus"]]
#        t_block = ref[:bus_block_map][switch["t_bus"]]
#        Graphs.add_edge!(ref[:block_graph], f_block, t_block)
#        ref[:block_graph_edge_map][Graphs.Edge(f_block, t_block)] = s
#        ref[:block_graph_edge_map][Graphs.Edge(t_block, f_block)] = s

#        if Int(switch["dispatchable"]) == Int(_PMD.YES) && Int(switch["status"]) == Int(_PMD.ENABLED)
#            push!(ref[:block_switches][f_block], s)
#            push!(ref[:block_switches][t_block], s)
#        end
#    end

    # Build block pairs for radiality constraints
    ref[:block_pairs] = filter(((x,y),)->x!=y, Set{Tuple{Int,Int}}(
            Set([(ref[:bus_block_map][sw["f_bus"]],ref[:bus_block_map][sw["t_bus"]]) for (_,sw) in ref[:switch]]),
    ))

#    ref[:neighbors] = Dict{Int,Vector{Int}}(i => Graphs.neighbors(ref[:block_graph], i) for i in Graphs.vertices(ref[:block_graph]))

#    ref[:switch_scores] = Dict{Int,Float64}(s => 0.0 for (s,_) in ref[:switch])
#    total_line_losses = sum(values(ref[:block_line_losses]))
#    for type in ["storage", "gen"]
#        for (id,obj) in ref[Symbol(type)]
#            if obj[_PMD.pmd_math_component_status[type]] != _PMD.pmd_math_component_status_inactive[type]
#                start_block = ref[:bus_block_map][obj["$(type)_bus"]]
#                paths = Graphs.enumerate_paths(Graphs.dijkstra_shortest_paths(ref[:block_graph], start_block))

#                for path in paths
#                    cumulative_weight = 0.0
#                    for (i,b) in enumerate(reverse(path[2:end]))
#                        block_line_losses = 1e-2 * ref[:block_line_losses][b]
#                        cumulative_weight += 1e-2 * ref[:block_weights][b]

#                        adjusted_cumulative_weight = cumulative_weight - (total_line_losses == 0.0 ? 0.0 : block_line_losses / total_line_losses)
#                        ref[:switch_scores][ref[:block_graph_edge_map][Graphs.Edge(path[end-i],b)]] += adjusted_cumulative_weight < 0 ? 0.0 : adjusted_cumulative_weight
#                    end
#                end
#            end
#        end
#    end
end

"""
    ref_add_load_blocks!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})

Ref extension to add load blocks to ref for all time steps
"""
function ref_add_load_blocks!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    _PMD.apply_pmd!(_ref_add_load_blocks!, ref, data; apply_to_subnetworks=true)
end

function ref_add_rdt!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    ref_add_branch_ne!(ref, data)
    ref_add_transformer_ne!(ref, data)
    ref_add_branch_harden!(ref, data)
    ref_add_gen_ne!(ref, data)
    ref_add_switch_inline_ne!(ref,data)
    ref_add_undamaged_branch!(ref,data)
    ref_add_damaged_tag!(ref,data)
    ref_add_global_constants!(ref,data)
    ref_add_load_blocks!(ref,data)
    _PMONM.ref_add_options!(ref,data)
end
