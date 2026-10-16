# Native reader for the full Boltzmann-two-term plasma mechanism path. The
# smaller isotropic reader remains in PlasmaIO.jl; this file reuses its
# validated scalar, species, equation, unit, and composition helpers.

function _plasma_rate_sections(root, phase, path, data_paths)
    selected = get(phase, "reactions", nothing)
    selected isa AbstractVector || _plasma_error("phase reactions must select named sections")
    out = Tuple{Any,Tuple{Float64,Float64,Float64,String},String}[]
    for (selection_index, item) in enumerate(selected)
        item isa AbstractDict && length(item) == 1 ||
            _plasma_error("phase reactions item $selection_index must select one section")
        raw_section, wanted = first(item)
        section = String(raw_section)
        source, key, context = root, section, "mechanism section $section"
        if occursin('/', section)
            endswith(section, "/reactions") ||
                _plasma_error("unsupported reaction section $section; expected <file>/reactions")
            file = section[1:end-length("/reactions")]
            isempty(file) && _plasma_error("invalid imported reaction section $section")
            imported = _plasma_import(file, path, data_paths)
            source = YAML.load_file(imported)
            source isa AbstractDict || _plasma_error("$imported is not a YAML mapping")
            key, context = "reactions", "$imported reactions"
        end
        entries = get(source, key, nothing)
        entries isa AbstractVector || _plasma_error("missing reaction section $context")
        chosen = if wanted == "all"
            entries
        else
            wanted isa AbstractVector ||
                _plasma_error("$context selection must be 'all' or a list of reaction IDs")
            picked = Any[]
            for id in wanted
                hits = findall(e -> e isa AbstractDict && get(e, "id", nothing) == id, entries)
                isempty(hits) && _plasma_error("reaction ID $id not found in $context")
                length(hits) == 1 || _plasma_error("reaction ID $id is ambiguous in $context")
                push!(picked, entries[only(hits)])
            end
            picked
        end
        units = _plasma_units(source)
        for (local_index, entry) in enumerate(chosen)
            entry isa AbstractDict || _plasma_error("$context entry $local_index must be a mapping")
            push!(out, (entry, units, "$context entry $local_index"))
        end
    end
    return out
end

function _plasma_rate_sides(equation, index, electron, context)
    eq = String(equation)
    if occursin("<=>", eq)
        sides, reversible = split(eq, "<=>"), true
    elseif occursin("=>", eq)
        sides, reversible = split(eq, "=>"), false
    else
        _plasma_error("$context needs a => or <=> reaction arrow")
    end
    length(sides) == 2 || _plasma_error("$context needs exactly one reaction arrow")
    # Expand parenthesized colliders so the shared helper can remove them from
    # stoichiometry while retaining generic or specific collider semantics.
    expanded = [replace(String(side), r"\(\+\s*([^)]+?)\s*\)" => s"+ \1") for side in sides]
    return _plasma_side(expanded[1], index, electron, context),
           _plasma_side(expanded[2], index, electron, context), reversible
end

function _plasma_rate_triplet(rate, order, units, context)
    rate isa AbstractDict || _plasma_error("$context must be a rate mapping")
    haskey(rate, "A") || _plasma_error("$context is missing A")
    order >= 1 || _plasma_error("$context has invalid total reaction order $order")
    lf, qf, tf, activation_unit = units
    A = _plasma_float(rate["A"], "$context A")
    A >= 0 || _plasma_error("negative A in $context")
    A *= (lf^3 / qf)^(order - 1) / tf
    b = _plasma_float(get(rate, "b", 0), "$context b")
    Ea_si = _plasma_activation(get(rate, "Ea", 0), activation_unit, "$context Ea")
    all(isfinite, (A, b, Ea_si)) || _plasma_error("nonfinite parameters in $context")
    return (A, b, Ea_si)
end

function _plasma_chebyshev(r, order, units, context)
    raw_temperature = get(r, "temperature-range", nothing)
    raw_pressure = get(r, "pressure-range", nothing)
    raw_data = get(r, "data", nothing)
    raw_temperature isa AbstractVector && length(raw_temperature) == 2 ||
        _plasma_error("$context needs a two-value temperature-range")
    raw_pressure isa AbstractVector && length(raw_pressure) == 2 ||
        _plasma_error("$context needs a two-value pressure-range")
    raw_data isa AbstractVector && !isempty(raw_data) ||
        _plasma_error("$context needs a nonempty Chebyshev data matrix")
    all(row -> row isa AbstractVector && !isempty(row), raw_data) ||
        _plasma_error("$context Chebyshev data must contain nonempty rows")
    n_pressure = length(first(raw_data))
    all(row -> length(row) == n_pressure, raw_data) ||
        _plasma_error("$context Chebyshev rows must have equal length")
    coefficients = Matrix{Float64}(undef, length(raw_data), n_pressure)
    for i in axes(coefficients, 1), j in axes(coefficients, 2)
        coefficients[i, j] = _plasma_float(raw_data[i][j], "$context coefficient [$i,$j]")
    end
    Tmin = _plasma_temperature(raw_temperature[1], "$context minimum temperature")
    Tmax = _plasma_temperature(raw_temperature[2], "$context maximum temperature")
    Tmin < Tmax || _plasma_error("$context temperature-range must be increasing")
    Pmin = _plasma_pressure(raw_pressure[1], "$context minimum pressure")
    Pmax = _plasma_pressure(raw_pressure[2], "$context maximum pressure")
    if Pmin == Pmax
        n_pressure == 1 ||
            _plasma_error("$context may use equal pressure bounds only with one coefficient column")
    else
        Pmin < Pmax || _plasma_error("$context pressure-range must be increasing")
    end
    order >= 1 || _plasma_error("$context has invalid total reaction order $order")
    lf, qf, tf, _ = units
    rate_factor = (lf^3 / qf)^(order - 1) / tf
    isfinite(rate_factor) && rate_factor > 0 ||
        _plasma_error("$context has an invalid rate-unit conversion")
    coefficients[1, 1] += log10(rate_factor)
    return coefficients, (Tmin, Tmax), (Pmin, Pmax)
end

function _plasma_troe_parameters(raw, context)
    raw isa AbstractDict || _plasma_error("$context Troe parameters must be a mapping")
    haskey(raw, "A") && haskey(raw, "T1") && haskey(raw, "T3") ||
        _plasma_error("$context Troe parameters require A, T1, and T3")
    a = _plasma_float(raw["A"], "$context Troe A")
    T1 = _plasma_temperature(raw["T1"], "$context Troe T1")
    T2 = haskey(raw, "T2") ? _plasma_temperature(raw["T2"], "$context Troe T2") : Inf
    T3 = _plasma_temperature(raw["T3"], "$context Troe T3")
    0 <= a <= 1 || _plasma_error("$context Troe A must be between zero and one")
    return (a, T1, T2, T3)
end

function _plasma_initial_boltzmann_energy(phase, eedf, state)
    temperatures = Tuple{Any,String}[]
    for (mapping, label) in ((state, "phase state"), (phase, "phase"))
        for key in ("electron-temperature", "Te")
            haskey(mapping, key) && push!(temperatures, (mapping[key], "$label $key"))
        end
    end
    has_mean = haskey(eedf, "mean-electron-energy")
    (has_mean ? 1 : 0) + length(temperatures) <= 1 ||
        _plasma_error("specify only one initial electron energy or temperature")
    has_mean && return _plasma_electron_energy(eedf["mean-electron-energy"], "mean-electron-energy")
    Te = isempty(temperatures) ? 0.001 :
        _plasma_temperature(first(only(temperatures)), last(only(temperatures)))
    return 1.5 * Te * _EEDF_BOLTZMANN / _EEDF_ELECTRON_CHARGE
end

function _plasma_matrix(rows::Vector{NTuple{N,Float64}}) where N
    out = zeros(Float64, length(rows), N)
    for i in eachindex(rows), j in 1:N
        out[i, j] = rows[i][j]
    end
    return out
end

function _plasma_boltzmann_mechanism(root, ph, path, data_paths, atomic_weights)
    data_paths isa AbstractVector || _plasma_error("data_paths must be a list")
    all(p -> p isa AbstractString, data_paths) || _plasma_error("data_paths entries must be strings")
    get(ph, "thermo", nothing) == "plasma" || _plasma_error("selected phase is not plasma")
    eedf = get(ph, "electron-energy-distribution", nothing)
    eedf isa AbstractDict || _plasma_error("phase has no EEDF")
    get(eedf, "type", nothing) == "Boltzmann-two-term" ||
        _plasma_error("selected phase does not use Boltzmann-two-term EEDF")
    raw_elements = get(ph, "elements", nothing)
    raw_elements isa AbstractVector && !isempty(raw_elements) || _plasma_error("phase needs elements")
    elements = String.(raw_elements)
    species_root = haskey(root, "species") ? root : merge(root, Dict("species" => Any[]))
    names, definitions = _plasma_species(species_root, ph, path, data_paths)
    elemental_matrix, molecular_weights, electron_index =
        _plasma_species_data(names, definitions, elements, atomic_weights)
    species_index = Dict(name => i for (i, name) in enumerate(names))
    electron = names[electron_index]
    selected = _plasma_rate_sections(root, ph, path, data_paths)
    n_species, n_reactions = length(names), length(selected)

    reactants = zeros(Int, n_species, n_reactions)
    products = zeros(Int, n_species, n_reactions)
    orders = zeros(Int, n_species, n_reactions)
    reversible = zeros(Bool, n_reactions)
    simple_thirdbody = falses(n_reactions)
    efficiencies = zeros(Float64, n_species, n_reactions)
    rate_types = zeros(UInt8, n_reactions)
    rate_parameters = zeros(Float64, 6, n_reactions)
    collision_energy = [Float64[] for _ in 1:n_reactions]
    cross_sections = [Float64[] for _ in 1:n_reactions]
    arrhenius = zeros(Float64, n_reactions, 3)
    low_rows, troe_rows = NTuple{3,Float64}[], NTuple{4,Float64}[]
    index_three_body, index_falloff, index_falloff_troe = Int64[], Int64[], Int64[]
    chebyshev_indices = Int[]
    chebyshev_coefficients = Matrix{Float64}[]
    chebyshev_temperature_ranges = Tuple{Float64,Float64}[]
    chebyshev_pressure_ranges = Tuple{Float64,Float64}[]
    supported = ("arrhenius", "elementary", "three-body", "falloff",
                 "two-temperature-plasma", "electron-collision-plasma", "chebyshev")

    for (j, (r, units, source_context)) in enumerate(selected)
        equation = get(r, "equation", nothing)
        equation isa AbstractString || _plasma_error("$source_context has no equation")
        context = "reaction $j $(repr(equation))"
        raw_kind = get(r, "type", "elementary")
        raw_kind isa AbstractString || _plasma_error("$context type must be a string")
        kind = lowercase(String(raw_kind))
        kind in supported || _plasma_error("unsupported type $kind in $context")
        lhs, rhs, arrow_reversible = _plasma_rate_sides(equation, species_index, electron, context)
        if haskey(r, "reversible")
            declared = r["reversible"]
            declared isa Bool || _plasma_error("$context reversible must be boolean")
            declared == arrow_reversible || _plasma_error("$context reversible conflicts with its arrow")
        end
        special = kind in ("two-temperature-plasma", "electron-collision-plasma", "chebyshev")
        arrow_reversible && special && _plasma_error("reversible $kind reaction is unsupported in $context")
        # Repeated reactants/products in specialized plasma equations are
        # stoichiometric participants. Only an explicit M may make those rows
        # third-body aware; general specific-collider inference remains useful
        # for the declared gas three-body and falloff families.
        thirdbody_kind = if kind in ("two-temperature-plasma", "chebyshev") &&
                get(lhs, "M", 0) == 0 && get(rhs, "M", 0) == 0
            "elementary"
        elseif kind == "arrhenius"
            "elementary"
        else
            kind
        end
        has_thirdbody, efficiency =
            _plasma_thirdbody!(lhs, rhs, r, thirdbody_kind, species_index, electron, context)
        kind in ("falloff", "three-body") && !has_thirdbody &&
            _plasma_error("$context $kind reaction has no resolvable collider")
        simple = has_thirdbody && kind != "falloff"
        simple_thirdbody[j] = simple
        has_thirdbody && (efficiencies[:, j] .= efficiency)
        reversible[j] = arrow_reversible
        orders[:, j] .= _plasma_orders(lhs, r, species_index, electron, context)
        for (name, coefficient) in lhs
            reactants[species_index[name], j] = coefficient
        end
        for (name, coefficient) in rhs
            products[species_index[name], j] = coefficient
        end
        forward_order = sum(@view orders[:, j])

        if kind in ("arrhenius", "elementary", "three-body")
            total_order = forward_order + (simple ? 1 : 0)
            parameters = _plasma_parameters(r, "arrhenius", total_order, units, context)
            rate_types[j] = 0x01
            rate_parameters[:, j] .= parameters
            arrhenius[j, :] .= (parameters[1], parameters[2], parameters[3] / 4184.0)
        elseif kind == "falloff"
            broadening = [key for key in ("SRI", "Tsang") if haskey(r, key)]
            isempty(broadening) || _plasma_error(
                "$context uses unsupported falloff parameterization $(join(broadening, ", "))")
            haskey(r, "rate-constant") &&
                _plasma_error("$context falloff reaction must use low/high-P rate constants")
            high = _plasma_rate_triplet(get(r, "high-P-rate-constant", nothing),
                                        forward_order, units, "$context high-P rate")
            low = _plasma_rate_triplet(get(r, "low-P-rate-constant", nothing),
                                       forward_order + 1, units, "$context low-P rate")
            rate_types[j] = 0x01
            rate_parameters[1:3, j] .= high
            arrhenius[j, :] .= (high[1], high[2], high[3] / 4184.0)
            push!(index_falloff, j)
            push!(low_rows, (low[1], low[2], low[3] / 4184.0))
            if haskey(r, "Troe")
                push!(troe_rows, _plasma_troe_parameters(r["Troe"], context))
                push!(index_falloff_troe, length(troe_rows))
            else
                push!(index_falloff_troe, -1)
            end
        elseif kind == "two-temperature-plasma"
            total_order = forward_order + (simple ? 1 : 0)
            rate_types[j] = 0x02
            rate_parameters[:, j] .= _plasma_parameters(r, kind, total_order, units, context)
        elseif kind == "electron-collision-plasma"
            has_thirdbody && _plasma_error("$context collision cannot use a third body")
            reactants[electron_index, j] > 0 || _plasma_error("$context needs an incident electron")
            rate_types[j] = 0x03
            collision_energy[j], cross_sections[j] = _eedf_table(r, context)
        else
            has_thirdbody && _plasma_error("$context Chebyshev reaction cannot use a third body")
            rate_types[j] = 0x04
            coefficients, temperature_range, pressure_range =
                _plasma_chebyshev(r, forward_order, units, context)
            push!(chebyshev_indices, j)
            push!(chebyshev_coefficients, coefficients)
            push!(chebyshev_temperature_ranges, temperature_range)
            push!(chebyshev_pressure_ranges, pressure_range)
        end
        kind != "electron-collision-plasma" &&
            (haskey(r, "energy-levels") || haskey(r, "cross-sections")) &&
            _plasma_error("noncollision $context cannot contain collision tables")
        simple && push!(index_three_body, j)
    end

    stoichiometry = Float64.(products - reactants)
    state = get(ph, "state", Dict{Any,Any}())
    state isa AbstractDict || _plasma_error("phase state must be a mapping")
    initial_x = _plasma_initial_x(state, names, species_index, electron, molecular_weights)
    eedf_model = read_eedf_model(path; phase=ph["name"], data_paths=data_paths)
    _plasma_validate(elemental_matrix, molecular_weights, stoichiometry,
                     rate_parameters, eedf_model.energy_edges, 2.0, initial_x)

    product_sparse, reactant_sparse = sparse(Float64.(products)), sparse(Float64.(reactants))
    order_sparse, vk = sparse(Float64.(orders)), sparse(stoichiometry)
    i_reactant = [Int64.(findall(!iszero, @view orders[:, j])) for j in 1:n_reactions]
    i_product = [Int64.(findall(>(0), @view products[:, j])) for j in 1:n_reactions]
    reaction = Reaction(product_sparse, reactant_sparse, order_sparse, reversible,
        arrhenius, _plasma_matrix(low_rows), _plasma_matrix(troe_rows),
        index_three_body, index_falloff, index_falloff_troe, sparse(efficiencies),
        i_reactant, i_product, n_reactions, vk, vec(Float64.(sum(vk, dims=1))),
        PlogData(Int64[], Int64[], Int64[1], Float64[], Int64[1], zeros(0, 3)),
        BlowersMaselData(Int64[], zeros(0, 4)))

    for (name, definition) in zip(names, definitions)
        species_thermo = get(definition, "thermo", nothing)
        species_thermo isa AbstractDict || _plasma_error("species $name has no thermo mapping")
        get(species_thermo, "model", "NASA7") == "NASA7" ||
            _plasma_error("Boltzmann plasma import currently requires NASA7 thermo for species $name")
    end
    thermo_document = Dict{Any,Any}(
        "phases" => Any[Dict{Any,Any}("species" => names)],
        "species" => definitions,
        "units" => get(root, "units", Dict{Any,Any}()),
    )
    thermo = IdealGasThermo(thermo_document)
    thermal = _PlasmaThermalData(thermo, reaction, chebyshev_indices,
        chebyshev_coefficients, chebyshev_temperature_ranges,
        chebyshev_pressure_ranges, eedf_model)
    initial_temperature = _plasma_temperature(get(state, "T", 300), "initial temperature")
    initial_pressure = _plasma_pressure(get(state, "P", 101325), "initial pressure")
    initial_energy = _plasma_initial_boltzmann_energy(ph, eedf, state)

    return PlasmaMechanism(String(get(ph, "name", "plasma")), n_species, n_reactions,
        names, elements, molecular_weights, elemental_matrix, electron_index,
        reactants, products, orders, stoichiometry, simple_thirdbody, efficiencies,
        rate_types, rate_parameters, collision_energy, cross_sections,
        copy(eedf_model.energy_edges), 2.0, initial_temperature, initial_pressure,
        initial_energy, initial_x, thermal)
end
