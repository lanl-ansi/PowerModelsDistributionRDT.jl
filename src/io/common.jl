
function parse_file(file::String; kwargs...)
    filetype = split(lowercase(file), '.')[end]
    io = open(file, "r")
    if filetype == "m"
        pmd_data = _PMD.parse_matlab(io)
    elseif filetype == "dss"
        pmd_data = _PMD.parse_file(file; kwargs...)
    elseif filetype == "json"
        pmd_data = parse_json(io)
    end
    # add_load_weights!(pmd_data)
    return pmd_data
end
