# Loader for electron-collision cross-section data in Cantera-style plasma YAML.
# Tables are used verbatim: energy levels in eV, cross sections in m^2 as written
# in the source file (no top-level unit conversion is applied here).

const _EEDF_KINDS = Dict{String,ElectronCollisionKind}(
    "effective" => EffectiveCollision,
    "elastic" => ElasticCollision,
    "excitation" => ExcitationCollision,
    "ionization" => IonizationCollision,
    "attachment" => AttachmentCollision,
)

_eedf_error(msg) = throw(ArgumentError(msg))

function _eedf_float(x, what)
    x isa Real || _eedf_error("$what must be a real number, got $(repr(x))")
    v = Float64(x)
    isfinite(v) || _eedf_error("$what must be finite")
    return v
end

function _eedf_table(entry, ctx)
    e = get(entry, "energy-levels", nothing)
    s = get(entry, "cross-sections", nothing)
    e isa AbstractVector || _eedf_error("$ctx: missing or invalid energy-levels table")
    s isa AbstractVector || _eedf_error("$ctx: missing or invalid cross-sections table")
    length(e) == length(s) ||
        _eedf_error("$ctx: energy-levels and cross-sections length mismatch " *
                    "($(length(e)) vs $(length(s)))")
    length(e) >= 2 || _eedf_error("$ctx: at least two cross-section table points required")
    energy = [_eedf_float(v, "$ctx energy level") for v in e]
    cross = [_eedf_float(v, "$ctx cross section") for v in s]
    any(<(0.0), energy) && _eedf_error("$ctx: energies must be nonnegative")
    any(<(0.0), cross) && _eedf_error("$ctx: cross sections must be nonnegative")
    all(>(0.0), diff(energy)) ||
        _eedf_error("$ctx: energy levels must be strictly increasing")
    return energy, cross
end

function _eedf_threshold(entry, ctx)
    thr = _eedf_float(get(entry, "threshold", 0.0), "$ctx threshold")
    thr < 0.0 && _eedf_error("$ctx: threshold must be nonnegative")
    return thr
end

# Parse one side of a reaction equation into species => stoichiometric count.
# Splits on the spaced " + " separator so ionic names such as "O2+" survive.
function _eedf_side(side, ctx)
    counts = Dict{String,Int}()
    for term in split(side, " + "; keepempty=false)
        term = strip(term)
        isempty(term) && continue
        m = match(r"^(\d+)\s+(.+)$", term)
        name = m === nothing ? term : strip(m.captures[2])
        n = m === nothing ? 1 : parse(Int, m.captures[1])
        isempty(name) && _eedf_error("$ctx: empty species name in '$side'")
        counts[name] = get(counts, name, 0) + n
    end
    isempty(counts) && _eedf_error("$ctx: no species found in '$side'")
    return counts
end

function _eedf_reaction_kind(react, prod, electron, compositions, ctx)
    react == prod && return EffectiveCollision

    # Files predating phase species declarations conventionally used Electron.
    # Keep their stoichiometric inference as a compatibility path; phase-aware
    # files use the selected species compositions, matching Cantera's charge test.
    if compositions === nothing
        pe = get(prod, electron, 0)
        re = get(react, electron, 0)
        pe > re && return IonizationCollision
        pe < re && return AttachmentCollision
        return ExcitationCollision
    end

    positive = false
    negative = false
    for name in keys(prod)
        name == electron && continue
        haskey(compositions, name) ||
            _eedf_error("$ctx: product species $(repr(name)) is not selected by the phase")
        charge = -get(compositions[name], "E", 0.0)
        positive |= charge > 0.0
        negative |= charge < 0.0
    end
    positive && negative &&
        _eedf_error("$ctx: both positive and negative product ions make collision kind ambiguous")
    positive && return IonizationCollision
    negative && return AttachmentCollision
    return ExcitationCollision
end

function _eedf_species_context(root, phase, path, data_paths)
    haskey(phase, "species") || return "Electron", nothing
    species_root = root
    if !haskey(root, "species")
        species_root = copy(root)
        species_root["species"] = Any[]
    end
    names, defs = _plasma_species(species_root, phase, path, data_paths)
    compositions = Dict{String,Dict{String,Float64}}()
    electrons = String[]
    for (name, def) in zip(names, defs)
        comp = get(def, "composition", nothing)
        comp isa AbstractDict ||
            _eedf_error("selected species '$name' has no composition mapping")
        parsed = Dict{String,Float64}()
        nonzero = 0
        for (raw_element, raw_count) in comp
            element = String(raw_element)
            count = _eedf_float(raw_count, "species '$name' composition $element")
            parsed[element] = count
            count != 0.0 && (nonzero += 1)
        end
        compositions[name] = parsed
        nonzero == 1 && get(parsed, "E", 0.0) == 1.0 && push!(electrons, name)
    end
    length(electrons) == 1 ||
        _eedf_error("phase must select exactly one pure {E: 1} electron species")
    return only(electrons), compositions
end

function _eedf_section(root, path, data_paths, section)
    source = root
    key = section
    context = "mechanism section '$section'"
    if occursin('/', section)
        endswith(section, "/reactions") ||
            _eedf_error("unsupported imported reaction section '$section'; expected <file>/reactions")
        file = section[1:end-length("/reactions")]
        isempty(file) && _eedf_error("invalid imported reaction section '$section'")
        imported = _plasma_import(file, path, data_paths)
        source = YAML.load_file(imported)
        source isa AbstractDict || _eedf_error("$imported is not a YAML mapping")
        key = "reactions"
        context = "$imported reactions"
    end
    entries = get(source, key, nothing)
    entries isa AbstractVector || _eedf_error("missing or invalid $context")
    return entries, context
end

function _eedf_select(entries, wanted, context)
    wanted == "all" && return collect(entries)
    wanted isa AbstractVector ||
        _eedf_error("$context selection must be 'all' or a list of reaction IDs")
    selected = Any[]
    for raw_id in wanted
        raw_id isa AbstractString || _eedf_error("$context reaction IDs must be strings")
        hits = findall(e -> e isa AbstractDict && get(e, "id", nothing) == raw_id,
                       entries)
        isempty(hits) && _eedf_error("reaction ID $(repr(raw_id)) not found in $context")
        length(hits) == 1 ||
            _eedf_error("reaction ID $(repr(raw_id)) is ambiguous in $context")
        push!(selected, entries[only(hits)])
    end
    return selected
end

function _eedf_selected_sections(root, phase, path, data_paths)
    selected = get(phase, "reactions", nothing)
    if selected === nothing
        # Compatibility with the standalone Phelps-style layout, which predates
        # phase reaction selectors.
        return [(get(root, "electron-collisions", Any[]), :root,
                 "electron-collisions"),
                (get(root, "reactions", Any[]), :reaction, "reactions")]
    end
    selected isa AbstractVector || _eedf_error("phase reactions must be a list")
    # PlasmaPhase installs root collision tables before selecting kinetics.
    root_entries = get(root, "electron-collisions", Any[])
    root_entries isa AbstractVector || _eedf_error("electron-collisions must be a list")
    sections = Tuple{Vector,Symbol,String}[(collect(root_entries), :root, "electron-collisions")]
    for (i, item) in enumerate(selected)
        item isa AbstractDict && length(item) == 1 ||
            _eedf_error("phase reactions item $i must select one section")
        raw_section, wanted = first(item)
        raw_section isa AbstractString ||
            _eedf_error("phase reactions item $i has a non-string section name")
        section = String(raw_section)
        entries, context = _eedf_section(root, path, data_paths, section)
        chosen = _eedf_select(entries, wanted, context)
        style = section == "electron-collisions" ? :root : :reaction
        push!(sections, (chosen, style, context))
    end
    return sections
end

function _eedf_grid(phase, ctx)
    eedf = get(phase, "electron-energy-distribution", nothing)
    eedf isa AbstractDict ||
        _eedf_error("$ctx: phase has no electron-energy-distribution")
    type = get(eedf, "type", nothing)
    type == "Boltzmann-two-term" ||
        _eedf_error("$ctx: unsupported electron-energy-distribution type " *
                    "$(repr(type)); only Boltzmann-two-term is supported")
    levels = get(eedf, "energy-levels", nothing)
    levels isa AbstractVector ||
        _eedf_error("$ctx: missing or invalid EEDF energy-levels")
    edges = [_eedf_float(v, "$ctx EEDF energy level") for v in levels]
    length(edges) >= 3 ||
        _eedf_error("$ctx: EEDF grid needs at least 3 energy levels")
    any(<(0.0), edges) &&
        _eedf_error("$ctx: EEDF energy levels must be nonnegative")
    all(>(0.0), diff(edges)) ||
        _eedf_error("$ctx: EEDF energy levels must be strictly increasing")
    return edges
end

function _eedf_phase(root, phase)
    phases = get(root, "phases", nothing)
    phases isa AbstractVector && !isempty(phases) ||
        _eedf_error("file has no phases section")
    if phase === nothing
        idx = findfirst(p -> p isa AbstractDict &&
                             haskey(p, "electron-energy-distribution"), phases)
        idx === nothing &&
            _eedf_error("no phase declares an electron-energy-distribution; " *
                        "pass phase=<name> to select one explicitly")
        return phases[idx]
    end
    idx = findfirst(p -> p isa AbstractDict && get(p, "name", nothing) == phase,
                    phases)
    idx === nothing && _eedf_error("phase $(repr(phase)) not found")
    return phases[idx]
end

"""
    read_eedf_model(path; phase=nothing, data_paths=String[]) -> EEDFModel

Read electron-collision cross sections from a Cantera-style plasma YAML file.
Only `Boltzmann-two-term` EEDF definitions are supported. Collision reaction
sections are resolved from the selected phase, including `<file>/reactions`
imports searched through `data_paths`. Root-level `electron-collisions` entries
require an explicit `kind`. Reaction kinds use an explicit `kind` when present;
otherwise they are inferred from selected species composition and product charge.
For reaction entries, zero thresholds are inferred from the first positive
tabulated energy only for excitation, ionization, and attachment collisions.
The selector-free `electron-collisions` layout keeps its legacy `Electron`
identity and explicit root thresholds without resolving thermodynamic imports.
"""
function read_eedf_model(path; phase=nothing, data_paths=String[])
    data_paths isa AbstractVector || _eedf_error("data_paths must be a list")
    all(p -> p isa AbstractString, data_paths) ||
        _eedf_error("data_paths entries must be strings")
    root = YAML.load_file(path)
    root isa AbstractDict || _eedf_error("$(path): not a YAML mapping")
    selected_phase = _eedf_phase(root, phase)
    edges = _eedf_grid(selected_phase, "read_eedf_model")
    legacy_root = !haskey(selected_phase, "reactions") &&
                  haskey(root, "electron-collisions")
    electron, compositions = legacy_root ? ("Electron", nothing) :
        _eedf_species_context(root, selected_phase, path, data_paths)

    collisions = ElectronCollision[]
    targets = String[]
    seen_background = Set{String}()

    function add!(target, kind, thr, energy, cross, origin)
        if kind === EffectiveCollision || kind === ElasticCollision
            key = target
            key in seen_background &&
                _eedf_error("duplicate $(kind) record for target '$target'")
            push!(seen_background, key)
        end
        target in targets || push!(targets, target)
        push!(collisions,
              ElectronCollision(target, kind, thr, energy, cross, origin))
    end

    for (entries, style, section_context) in
            _eedf_selected_sections(root, selected_phase, path, data_paths)
        entries isa AbstractVector || _eedf_error("$section_context must be a list")
        for entry in entries
            entry isa AbstractDict || _eedf_error("$section_context entry must be a mapping")
            if style === :root
                ctx = "$section_context entry $(repr(get(entry, "target", nothing)))"
                target = get(entry, "target", nothing)
                target isa AbstractString && !isempty(strip(target)) && target != electron ||
                    _eedf_error("$section_context entry: missing or invalid target")
                compositions === nothing || haskey(compositions, String(target)) ||
                    _eedf_error("$ctx: target is not selected by the phase")
                raw_kind = get(entry, "kind", nothing)
                haskey(_EEDF_KINDS, raw_kind) ||
                    _eedf_error("$ctx: explicit kind required, one of " *
                                "$(join(sort(collect(keys(_EEDF_KINDS))), ", "))")
                kind = _EEDF_KINDS[raw_kind]
                energy, cross = _eedf_table(entry, ctx)
                threshold = _eedf_threshold(entry, ctx)
            else
                get(entry, "type", nothing) == "electron-collision-plasma" || continue
                eq = get(entry, "equation", nothing)
                eq isa AbstractString ||
                    _eedf_error("electron-collision-plasma reaction: missing equation")
                ctx = "reaction '$eq'"
                occursin("<=>", eq) &&
                    _eedf_error("$ctx: reversible electron collisions are unsupported")
                sides = split(eq, "=>")
                length(sides) == 2 || _eedf_error("$ctx: malformed equation")
                react = _eedf_side(sides[1], ctx)
                prod = _eedf_side(sides[2], ctx)
                if compositions !== nothing
                    for name in union(keys(react), keys(prod))
                        haskey(compositions, name) ||
                            _eedf_error("$ctx: species $(repr(name)) is not selected by the phase")
                    end
                end
                get(react, electron, 0) == 1 ||
                    _eedf_error("$ctx: one incident $electron required")
                target_names = filter(!=(electron), collect(keys(react)))
                length(target_names) == 1 ||
                    _eedf_error("$ctx: expected exactly one collision target besides " *
                                "$electron, got $(join(target_names, ", "))")
                target = only(target_names)
                react[target] == 1 || _eedf_error("$ctx: one target particle required")
                raw_kind = get(entry, "kind", nothing)
                if raw_kind === nothing
                    kind = _eedf_reaction_kind(react, prod, electron, compositions, ctx)
                else
                    haskey(_EEDF_KINDS, raw_kind) ||
                        _eedf_error("$ctx: unsupported collision kind $(repr(raw_kind))")
                    kind = _EEDF_KINDS[raw_kind]
                end
                energy, cross = _eedf_table(entry, ctx)
                threshold = _eedf_threshold(entry, ctx)
            end

            if style === :reaction && threshold == 0.0 &&
                    kind in (ExcitationCollision, IonizationCollision,
                             AttachmentCollision)
                idx = findfirst(>(0.0), energy)
                idx === nothing &&
                    _eedf_error("$ctx: cannot infer threshold; no positive energy level")
                threshold = energy[idx]
            end
            add!(String(target), kind, threshold, energy, cross, style)
        end
    end

    isempty(targets) &&
        _eedf_error("$(path): no electron collisions found; no participating targets")
    return EEDFModel(edges, collisions, targets)
end
