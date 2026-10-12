function read_species_basics(yaml)
    n_species = length(yaml["phases"][1]["species"])
    n_reactions = length(get(yaml,"reactions",Any[]))
    species_names = yaml["phases"][1]["species"]
    elements = yaml["phases"][1]["elements"]
    n_elements = length(elements)

    ele_matrix = zeros(n_elements, n_species)

    _species_names =
        [yaml["species"][i]["name"] for i = 1:length(yaml["species"])]

    for (i, species) in enumerate(species_names)
        spec = yaml["species"][findfirst(x -> x == species, _species_names)]
        for j = 1:n_elements
            if haskey(spec["composition"], elements[j])
                ele_matrix[j, i] = spec["composition"][elements[j]]
            end
        end
    end

    return n_species, n_reactions, species_names, elements, n_elements, ele_matrix
end

function _sidecar_string(npz, key)
    return String(vec(UInt8.(npz[key])))
end

function _matches_mechanism_hash(expected, mechanism)
    raw = read(mechanism)
    expected == bytes2hex(SHA.sha256(raw)) && return true
    # Git may change only the YAML line endings when checking out a sidecar
    # generated on another OS. Keep all other bytes significant.
    lf = replace(String(raw), "\r\n" => "\n")
    expected == bytes2hex(SHA.sha256(lf)) && return true
    return expected == bytes2hex(SHA.sha256(replace(lf, "\n" => "\r\n")))
end

function _phase_imports(phase)
    files = String[]
    for kind in ("species", "reactions")
        selectors = get(phase,kind,Any[])
        selectors isa AbstractVector || continue
        for item in selectors
            item isa AbstractDict || continue
            for key in keys(item)
                key isa AbstractString || throw(ArgumentError("phase section names must be strings"))
                occursin('/',key) && push!(files,first(rsplit(String(key),'/';limit=2)))
            end
        end
    end
    sort!(unique!(files))
end

function _sidecar_selection(npz)
    get(npz,"sidecar_format_utf8",nothing) === nothing && return nothing
    _sidecar_string(npz,"sidecar_format_utf8") == "arrhenius-sidecar-v4" || return nothing
    haskey(npz,"phase_selection_utf8") || throw(ArgumentError("v4 sidecar needs phase selection metadata"))
    selection = try
        YAML.load(_sidecar_string(npz,"phase_selection_utf8"))
    catch
        throw(ArgumentError("invalid phase selection metadata"))
    end
    selection isa AbstractDict || throw(ArgumentError("phase selection metadata must be a mapping"))
    names = get(selection,"species_names",nothing)
    names isa AbstractVector && !isempty(names) && all(x -> x isa AbstractString,names) &&
        length(unique(names)) == length(names) || throw(ArgumentError("invalid sidecar species names"))
    nr = get(selection,"n_reactions",nothing)
    nr isa Integer && !(nr isa Bool) && nr >= 0 || throw(ArgumentError("invalid selected reaction count"))
    deps = get(selection,"dependencies",nothing)
    deps isa AbstractDict && all(k isa AbstractString && v isa AbstractString &&
        occursin(r"^[0-9a-f]{64}$",v) for (k,v) in deps) ||
        throw(ArgumentError("invalid sidecar dependency hashes"))
    selection
end

function _validate_sidecar_metadata(npz, mechanism; data_paths=String[])
    if haskey(npz, "sidecar_format_utf8")
        format = _sidecar_string(npz, "sidecar_format_utf8")
        format in ("arrhenius-sidecar-v2", "arrhenius-sidecar-v3", "arrhenius-sidecar-v4") ||
            throw(ArgumentError(
                "unsupported Arrhenius sidecar format: $format",
            ))
    end
    if haskey(npz, "source_sha256_utf8")
        expected = _sidecar_string(npz, "source_sha256_utf8")
        _matches_mechanism_hash(expected, mechanism) || throw(ArgumentError(
            "the Arrhenius sidecar does not match $mechanism; regenerate " *
            "$(mechanism).npz from the current YAML mechanism",
        ))
    end
    selection = _sidecar_selection(npz)
    if selection !== nothing
        haskey(npz,"source_sha256_utf8") || throw(ArgumentError("v4 sidecar needs the source hash"))
        root = YAML.load_file(mechanism)
        deps = selection["dependencies"]
        Set(keys(deps)) == Set(_phase_imports(root["phases"][1])) ||
            throw(ArgumentError("sidecar dependency list does not match the selected phase; regenerate it"))
        for (file,digest) in deps
            imported = _plasma_import(file,mechanism,data_paths)
            _matches_mechanism_hash(digest,imported) || throw(ArgumentError(
                "sidecar dependency $file has changed; regenerate $(mechanism).npz"))
        end
    end
    return nothing
end

function _solution_species_yaml(yaml, mechanism, data_paths)
    sources = String[]
    names,defs = _plasma_species(yaml,yaml["phases"][1],mechanism,data_paths; source_paths=sources)
    isempty(names) && throw(ArgumentError("at least one species is required"))
    # Imported species retain their source unit context, including SI defaults;
    # the receiving phase's units must not reinterpret their numeric thermo data.
    roots = Dict{String,Any}()
    for i in eachindex(defs)
        abspath(sources[i]) == abspath(mechanism) && continue
        root = get!(roots,sources[i]) do
            YAML.load_file(sources[i])
        end
        units = merge(Dict("temperature"=>"K","pressure"=>"Pa","energy"=>"J","quantity"=>"kmol"),
            get(root,"units",Dict()),get(defs[i],"units",Dict()))
        defs[i] = merge(defs[i],Dict("units"=>units))
    end
    normalized = copy(yaml)
    normalized["phases"] = copy(yaml["phases"])
    normalized["phases"][1] = merge(yaml["phases"][1],Dict("species"=>names))
    normalized["species"] = defs
    normalized
end

function _validate_kinetics_dimensions(npz, ns, nr)
    shapes = (("molecular_weights",(ns,)),("efficiencies_coeffs",(ns,nr)),
        ("product_stoich_coeffs",(ns,nr)),("reactant_stoich_coeffs",(ns,nr)),
        ("reactant_orders",(ns,nr)),("is_reversible",(nr,)),("Arrhenius_coeffs",(nr,3)))
    for (key,shape) in shapes
        haskey(npz,key) && size(npz[key]) == shape ||
            throw(ArgumentError("sidecar $key must have shape $shape; regenerate the sidecar"))
    end
    nothing
end

"""
    CreateSolution(mech; data_paths=String[])
    
Load the first gas phase from a YAML mechanism and its preprocessed NPZ sidecar.
Imported species/reaction sections require a v4 sidecar. Imports are resolved
relative to the mechanism first, then in `data_paths`. Thermochemistry and
kinetics are evaluated in Julia after loading.
"""
function CreateSolution(mech; data_paths=String[])
    yaml = YAML.load_file(mech)
    npz = npzread("$mech.npz")
    _validate_sidecar_metadata(npz, mech; data_paths)
    selection = _sidecar_selection(npz)
    phase = yaml["phases"][1]
    if selection === nothing && (!(get(phase,"species",nothing) isa AbstractVector) ||
            !all(x -> x isa AbstractString,phase["species"]) ||
            get(phase,"reactions","all") != "all")
        throw(ArgumentError("selected/imported sections require a v4 sidecar; regenerate $(mech).npz"))
    end
    yaml = _solution_species_yaml(yaml,mech,data_paths)

    #### Basics
    n_species, n_reactions, species_names,
    elements, n_elements, ele_matrix = read_species_basics(yaml)
    if selection !== nothing
        species_names == selection["species_names"] || throw(ArgumentError("sidecar species order does not match the phase"))
        n_reactions = Int(selection["n_reactions"])
    end

    #### Thermo
    if yaml["phases"][1]["thermo"] == "ideal-gas" # switch to work with Cantera standard
        thermo = IdealGasThermo(yaml)
    else
        constructorThermo = Symbol(yaml["phases"][1]["thermo"], :Thermo)
        thermo = @eval($constructorThermo)(yaml)
    end


    #### Kinetic data

    # Inert phases have no reaction arrays. Empty NPZ entries are omitted by the
    # exporter for compatibility with NPZ.jl; restore their unambiguous shapes.
    if n_reactions == 0
        npz = Dict{String,Any}(npz)
        for key in ("efficiencies_coeffs","product_stoich_coeffs","reactant_stoich_coeffs","reactant_orders")
            haskey(npz,key) || (npz[key] = zeros(n_species,0))
        end
        haskey(npz,"is_reversible") || (npz["is_reversible"] = Bool[])
        haskey(npz,"Arrhenius_coeffs") || (npz["Arrhenius_coeffs"] = zeros(0,3))
    end
    _validate_kinetics_dimensions(npz,n_species,n_reactions)
    MW = vec(Float64.(npz["molecular_weights"]))
    efficiencies_coeffs_full = Matrix{Float64}(npz["efficiencies_coeffs"])
    product_stoich_coeffs = sparse(Float64.(npz["product_stoich_coeffs"]))
    reactant_stoich_coeffs = sparse(Float64.(npz["reactant_stoich_coeffs"]))
    reactant_orders = sparse(Float64.(npz["reactant_orders"]))
    is_reversible = Vector{Bool}(vec(Bool.(npz["is_reversible"])))
    Arrhenius_coeffs = Matrix{Float64}(npz["Arrhenius_coeffs"])
    if haskey(npz, "Arrhenius_A0")
        Arrhenius_A0 = npz["Arrhenius_A0"]
        Arrhenius_b0 = npz["Arrhenius_b0"]
        Arrhenius_Ea0 = npz["Arrhenius_Ea0"]
    else
        Arrhenius_A0 = []
        Arrhenius_b0 = []
        Arrhenius_Ea0 = []
    end

    if haskey(npz, "Troe_A")
        Troe_A = npz["Troe_A"]
        Troe_T1 = npz["Troe_T1"]
        Troe_T2 = npz["Troe_T2"]
        Troe_T3 = npz["Troe_T3"]
    else
        Troe_A = []
        Troe_T1 = []
        Troe_T2 = []
        Troe_T3 = []
    end
    Arrhenius_0 = Matrix{Float64}(
        hcat(Arrhenius_A0, Arrhenius_b0, Arrhenius_Ea0),
    )
    Troe_ = Matrix{Float64}(hcat(Troe_A, Troe_T1, Troe_T2, Troe_T3))

    has_plog = selection === nothing && any(
        get(reaction, "type", "") == "pressure-dependent-Arrhenius"
        for reaction in get(yaml,"reactions",Any[])
    )
    if haskey(npz, "Plog_reaction_indices")
        plog_reaction_indices = vec(Int64.(npz["Plog_reaction_indices"]))
        plog_collider_indices = haskey(npz, "Plog_collider_indices") ?
            vec(Int64.(npz["Plog_collider_indices"])) :
            zeros(Int64, length(plog_reaction_indices))
        length(plog_collider_indices) == length(plog_reaction_indices) ||
            throw(ArgumentError(
                "Plog_collider_indices must contain one entry per PLOG reaction",
            ))
        plog = PlogData(
            plog_reaction_indices,
            plog_collider_indices,
            vec(Int64.(npz["Plog_group_offsets"])),
            vec(Float64.(npz["Plog_pressures"])),
            vec(Int64.(npz["Plog_rate_offsets"])),
            Matrix{Float64}(npz["Plog_Arrhenius"]),
        )
    elseif has_plog
        throw(ArgumentError(
            "mechanism contains pressure-dependent-Arrhenius reactions, but " *
            "its sidecar has no PLOG data; regenerate it with the Cantera 3.2 " *
            "Arrhenius sidecar exporter",
        ))
    else
        plog = PlogData(
            Int64[], Int64[], Int64[1], Float64[], Int64[1], zeros(0, 3),
        )
    end

    if selection !== nothing || any(haskey(npz, key) for key in ("index_three_body", "index_falloff", "index_falloff_Troe"))
        # The Cantera API also identifies implicit or explicit colliders whose
        # YAML reaction lacks a `type: three-body` field.
        index_three_body = vec(Int64.(get(npz, "index_three_body", Int64[])))
        index_falloff = vec(Int64.(get(npz, "index_falloff", Int64[])))
        index_falloff_Troe = vec(Int64.(get(npz, "index_falloff_Troe", Int64[])))
        length(index_falloff) == length(index_falloff_Troe) ||
            throw(ArgumentError("one Troe index is required per falloff reaction"))
    else
        index_three_body = Int64[]
        index_falloff = Int64[]
        index_falloff_Troe = Int64[]
        j = 1
        for i = 1:n_reactions
            reaction = yaml["reactions"][i]
            if haskey(reaction, "type")
                if reaction["type"] == "three-body"
                    push!(index_three_body, i)
                end
                if reaction["type"] == "falloff"
                    push!(index_falloff, i)
                    if haskey(reaction, "Troe")
                        push!(index_falloff_Troe, j)
                        j = j + 1
                    else
                        push!(index_falloff_Troe, -1)
                    end
                end
            end
        end
    end

    all(i -> 1 <= i <= n_reactions,[index_three_body;index_falloff]) &&
        all(i -> i == -1 || 1 <= i <= size(Troe_,1),index_falloff_Troe) &&
        size(Arrhenius_0,1) == length(index_falloff) ||
        throw(ArgumentError("invalid sidecar third-body/falloff indices or dimensions"))

    i_reactant = Vector{Vector{Int64}}()
    i_product = Vector{Vector{Int64}}()
    for i = 1:n_reactions
        push!(i_reactant, findall(!iszero, reactant_orders[:, i]))
        push!(i_product, findall(product_stoich_coeffs[:, i] .> 0.01))
    end

    vk = sparse(product_stoich_coeffs - reactant_stoich_coeffs)
    vk_sum = vec(Float64.(sum(vk, dims=1)))

    for i in 1:n_reactions
        if !((i in index_three_body) | (i in index_falloff))
            efficiencies_coeffs_full[:, i] .= 0.0
        end
    end
    efficiencies_coeffs = sparse(efficiencies_coeffs_full)

    bm_indices = vec(Int64.(get(npz,"BlowersMasel_reaction_indices",Int64[])))
    bm_coefficients = Matrix{Float64}(get(npz,"BlowersMasel_coefficients",zeros(0,4)))
    size(bm_coefficients) == (length(bm_indices),4) &&
        all(i -> 1 <= i <= n_reactions,bm_indices) || throw(ArgumentError("invalid Blowers–Masel sidecar dimensions"))
    for row in eachrow(bm_coefficients)
        BlowersMaselRate(row...)
    end
    if selection === nothing && isempty(bm_indices) && any(get(r,"type","") == "Blowers-Masel" for r in get(yaml,"reactions",Any[]))
        throw(ArgumentError("regenerate the sidecar to include Blowers–Masel parameters"))
    end
    reaction = Reaction(
        product_stoich_coeffs,
        reactant_stoich_coeffs,
        reactant_orders,
        is_reversible,
        Arrhenius_coeffs,
        Arrhenius_0,
        Troe_,
        index_three_body,
        index_falloff,
        index_falloff_Troe,
        efficiencies_coeffs,
        i_reactant,
        i_product,
        n_reactions,
        vk,
        vk_sum,
        plog,
        BlowersMaselData(bm_indices,bm_coefficients),
    )


    #### Transport data

    ionized = get(yaml["phases"][1],"transport","") == "ionized-gas"
    if ionized && haskey(npz,"species_viscosities_poly")
        haskey(npz,"transport_model_utf8") &&
            _sidecar_string(npz,"transport_model_utf8") == "ionized-gas" ||
            throw(ArgumentError("ionized transport requires model-specific fits; regenerate the sidecar"))
    end
    if haskey(npz, "species_viscosities_poly")
        species_viscosities_poly =
            Matrix{Float64}(npz["species_viscosities_poly"])
        thermal_conductivity_poly =
            Matrix{Float64}(npz["thermal_conductivity_poly"])
        binary_diff_coeffs_poly =
            Matrix{Float64}(npz["binary_diff_coeffs_poly"])
        poly_order = size(species_viscosities_poly)[1]
    else
        species_viscosities_poly = zeros(2, 2)
        thermal_conductivity_poly = zeros(2, 2)
        binary_diff_coeffs_poly = zeros(2, 2)
        poly_order = 6
    end

    trans = Transport(poly_order,
                      species_viscosities_poly, 
                      thermal_conductivity_poly, 
                      binary_diff_coeffs_poly,
                      ionized ? :ionized_gas : :mixture_averaged)

    gas = Solution(
        n_species,
        n_reactions,
        MW,
        species_names,
        elements,
        ele_matrix,
        thermo,
        trans,
        reaction,
    )
    return gas
end

export CreateSolution
