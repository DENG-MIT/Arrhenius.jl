"Species mole fractions, with species along rows and grid points along columns."
function mole_fractions(f::AbstractPremixedFlame)
    X = mass_fractions(f) ./ f.gas.MW
    X ./= sum(X;dims=1)
    return X
end

"Mass density [kg/m³] at each grid point."
density(f::AbstractPremixedFlame) = f.pressure .* vec(transpose(f.gas.MW)*mole_fractions(f)) ./ (R .* temperature(f))

"Volumetric chemical heat release rate [W/m³], positive for heat release."
function heat_release_rate(f::AbstractPremixedFlame)
    w = FlameWorkspace(f)
    _flame_properties!(w,f,f.state;update_transport=false)
    return -vec(sum(w.h .* w.source;dims=1))
end

# Includes species order, chemistry, thermochemistry and transport fits. Display
# without truncation provides a deterministic fingerprint within a Julia version.
_flame_mechanism_signature(gas) = bytes2hex(SHA.sha256(repr(
    (gas.species_names,gas.MW,gas.thermo,gas.reaction,gas.trans);context=:limit=>false)))
_flame_utf8(value) = collect(codeunits(string(value)))

"""
    save_flame(path, flame; overwrite=false, basis=:mass)

Save a `.npz` restart snapshot or a `.csv` table. Snapshots retain the grid,
boundary conditions, temperature profile, transport options, spatial discretization,
profile-grid policy, and native state.
CSV columns contain position [m], velocity [m/s], temperature [K], density
[kg/m³], heat release [W/m³] and species fractions (`basis=:mass` or `:mole`).
"""
function save_flame(path::AbstractString,f::AbstractPremixedFlame;overwrite=false,basis=:mass)
    ispath(path) && !overwrite && throw(ArgumentError("output exists; use overwrite=true"))
    basis in (:mass,:mole) || throw(ArgumentError("basis must be :mass or :mole"))
    extension = lowercase(splitext(path)[2])
    if extension == ".csv"
        fractions = basis == :mass ? mass_fractions(f) : mole_fractions(f)
        data = hcat(f.grid,velocity(f),temperature(f),density(f),heat_release_rate(f),transpose(fractions))
        prefix = basis == :mass ? "Y_" : "X_"
        open(path,"w") do io
            println(io,join(vcat(["z_m","velocity_m_s","temperature_K","density_kg_m3","heat_release_W_m3"],
                prefix .* f.gas.species_names),","))
            for row in eachrow(data)
                println(io,join(row,","))
            end
        end
    elseif extension == ".npz"
        arrays = Dict{String,Any}(
            "format_utf8"=>_flame_utf8("arrhenius-premixed-flame-v1"),
            "kind_utf8"=>_flame_utf8(f isa FreeFlame ? "free" : "burner"),
            "mechanism_signature_utf8"=>_flame_utf8(_flame_mechanism_signature(f.gas)),
            "julia_version_utf8"=>_flame_utf8(VERSION),
            "species_names_utf8"=>_flame_utf8(join(f.gas.species_names,"\n")),
            "grid"=>f.grid,"state"=>f.state,"inlet_Y"=>f.inlet_Y,
            "conditions"=>[f.pressure,f.inlet_temperature,f.anchor,f.fixed_temperature],
            "transport_utf8"=>_flame_utf8(f.transport_model),
            "discretization_utf8"=>_flame_utf8(f.discretization),
            "gradient_basis_utf8"=>_flame_utf8(f.flux_gradient_basis),
            "soret"=>[Int(f.soret_enabled)])
        if f isa BurnerFlame
            arrays["mass_flux"] = [f.mass_flux]
            arrays["profile_grid_policy_utf8"] = _flame_utf8(f.profile_grid_policy)
            if !isempty(f.profile_positions)
                arrays["profile_positions"] = f.profile_positions
                arrays["profile_temperatures"] = f.profile_temperatures
            end
        end
        npzwrite(path,arrays)
    else
        throw(ArgumentError("flame output extension must be .csv or .npz"))
    end
    return path
end

"""
    restore_flame!(flame, path)

Restore a native `.npz` snapshot into a flame using the same mechanism and kind.
For multicomponent or Soret snapshots, supply matching `MultiTransportData` on
the destination flame. A restored state is marked converged only if its current
native steady residual satisfies the solver tolerance. Snapshots fingerprint
the numerical mechanism representation and are intended for the same Julia version.
Older archives without a profile-grid policy restore `:full_knots`; missing
spatial discretization retains historical finite-difference behavior. Invalid
archives leave the destination unchanged.
"""
function restore_flame!(f::AbstractPremixedFlame,path::AbstractString)
    # Stage on a shallow shell sharing the immutable mechanism/transport inputs.
    # All mutable snapshot arrays are replaced before any profile/state mutation.
    staged = typeof(f)(ntuple(i->getfield(f,i),fieldcount(typeof(f)))...)
    _restore_flame_snapshot!(staged,path)
    for name in fieldnames(typeof(f))
        setfield!(f,name,getfield(staged,name))
    end
    return f
end

function _restore_flame_snapshot!(f::AbstractPremixedFlame,path::AbstractString)
    arrays = npzread(path)
    text(key) = String(vec(UInt8.(arrays[key])))
    text("format_utf8") == "arrhenius-premixed-flame-v1" || throw(ArgumentError("unknown flame snapshot format"))
    text("kind_utf8") == (f isa FreeFlame ? "free" : "burner") || throw(ArgumentError("flame kind mismatch"))
    text("mechanism_signature_utf8") == _flame_mechanism_signature(f.gas) ||
        throw(ArgumentError("snapshot mechanism fingerprint does not match"))
    z,u,Y = vec(arrays["grid"]),arrays["state"],vec(arrays["inlet_Y"])
    n,N = f.gas.n_species,length(z)
    size(u) == (n+2,N) && length(Y)==n || throw(DimensionMismatch("snapshot state dimensions do not match"))
    N>=5 && all(isfinite,z) && all(>(0),diff(z)) && all(isfinite,u) &&
        all(t->.2<=t<=6,@view(u[1,:])) && minimum(@view(u[2:n+1,:]))>=-1e-7 &&
        all(>(0),@view(u[end,:])) || throw(ArgumentError("invalid snapshot state or grid"))
    all(isfinite,Y) && all(>=(0),Y) && isapprox(sum(Y),1;atol=1e-10) || throw(ArgumentError("invalid snapshot inlet composition"))
    conditions = vec(arrays["conditions"])
    length(conditions)==4 && all(isfinite,conditions) && conditions[1]>0 &&
        200<=conditions[2]<=6000 && isinteger(conditions[3]) && 2<=conditions[3]<N ||
        throw(ArgumentError("invalid snapshot boundary conditions"))
    discretization = haskey(arrays,"discretization_utf8") ? Symbol(text("discretization_utf8")) : :finite_difference
    discretization in (:finite_difference,:conservative) || throw(ArgumentError("invalid snapshot discretization"))
    f.grid,f.state,f.inlet_Y = z,u,Y
    f.discretization = discretization
    f.pressure,f.inlet_temperature = conditions[1:2]
    f.anchor,f.fixed_temperature = Int(conditions[3]),conditions[4]
    f.dependent_species = argmax(Y)
    if f isa BurnerFlame
        f.mass_flux = only(arrays["mass_flux"])
        isfinite(f.mass_flux) && f.mass_flux>0 || throw(ArgumentError("invalid snapshot mass flux"))
        # Never clear arrays shared with the destination shell.
        f.profile_positions,f.profile_temperatures,f.imposed_temperature = Float64[],Float64[],Float64[]
        f.profile_grid_policy = :full_knots
        haskey(arrays,"profile_positions") == haskey(arrays,"profile_temperatures") ||
            throw(ArgumentError("snapshot requires both profile arrays"))
    end
    set_transport!(f,Symbol(text("transport_utf8"));soret=Bool(only(arrays["soret"])),
        flux_gradient_basis=Symbol(text("gradient_basis_utf8")))
    if f isa BurnerFlame
        policy = haskey(arrays,"profile_grid_policy_utf8") ? Symbol(text("profile_grid_policy_utf8")) : :full_knots
        policy in (:full_knots,:adaptive) || throw(ArgumentError("invalid snapshot profile grid policy"))
        if haskey(arrays,"profile_positions")
            set_temperature_profile!(f,vec(arrays["profile_positions"]),
                vec(arrays["profile_temperatures"]);relative=false,grid_policy=policy)
        elseif policy != :full_knots
            throw(ArgumentError("adaptive snapshot requires a temperature profile"))
        end
    end
    f.converged = norm(flame_residual!(similar(f.state),f),Inf) < 1e-8
    return f
end

# Recognized column spellings for flame CSV tables: the native `save_flame`
# output and Cantera's `SolutionArray.write_csv`/`Sim1D.save` CSV format.
const _flame_csv_columns = (
    grid = ("z_m", "grid"),
    velocity = ("velocity_m_s", "velocity"),
    temperature = ("temperature_K", "T"),
    density = ("density_kg_m3", "D"),
)

"""
    read_flame_csv(path) -> (; columns, data)

Read a flame CSV into tabulated data. Both layouts are recognized: the native
`save_flame` table (`z_m,velocity_m_s,temperature_K,density_kg_m3,...`) and the
equivalent Cantera table (`grid,velocity,T,D,X_<species>...`). `columns` holds
the header entries and `data` is a `Matrix{Float64}` with one row per saved
grid point, in file order. Rows or columns may be pruned before the table is
passed to `set_initial_guess!`.
"""
function read_flame_csv(path::AbstractString)
    isfile(path) || throw(ArgumentError("no such flame CSV: $path"))
    lines = collect(eachline(path))
    length(lines) >= 2 || throw(ArgumentError("flame CSV needs a header and at least one data row"))
    columns = String[strip(field) for field in split(lines[1], ",")]
    all(!isempty, columns) || throw(ArgumentError("flame CSV contains an empty header entry"))
    data = zeros(length(lines) - 1, length(columns))
    for (i, line) in enumerate(@view lines[2:end])
        fields = split(line, ",")
        length(fields) == length(columns) ||
            throw(ArgumentError("flame CSV row $i has $(length(fields)) fields for $(length(columns)) columns"))
        for j in eachindex(fields)
            value = tryparse(Float64, strip(fields[j]))
            value === nothing &&
                throw(ArgumentError("flame CSV row $i column $(columns[j]) is not a number: $(fields[j])"))
            data[i, j] = value
        end
    end
    return (; columns, data)
end

"""
    set_initial_guess!(flame, data)

Load a restart initial guess from tabulated flame data, following Cantera's
`FlameBase.set_initial_guess(data=...)`. `data` is a CSV path (see
`read_flame_csv`) or a `(; columns, data)` table, for example a pruned subset
of a previous `save_flame` export. The flame adopts the data grid and takes its
temperature, velocity and species fractions from the table; the temperature
profile is shifted so its inlet matches the flame inlet, species absent from
the table start at zero (renormalizing the supplied fractions per grid point,
as Cantera's `SolutionArray` restore does; a table without any species columns
keeps the inlet composition), and boundary conditions are not overwritten. For
a `FreeFlame` the saved velocity becomes the mass-flux profile, and the
fixed-temperature point is reselected on the loaded profile; for a `BurnerFlame`
the prescribed mass flux is retained. The flame is marked unconverged; call
`solve!` to continue from the guess.
"""
function set_initial_guess!(f::AbstractPremixedFlame, data::AbstractString)
    return set_initial_guess!(f, read_flame_csv(data))
end

function set_initial_guess!(f::AbstractPremixedFlame, table::NamedTuple)
    f isa Union{FreeFlame,BurnerFlame} ||
        throw(ArgumentError("set_initial_guess! supports FreeFlame and BurnerFlame"))
    hasproperty(table, :columns) && hasproperty(table, :data) ||
        throw(ArgumentError("restart data must provide columns and data; see read_flame_csv"))
    columns = String[string(column) for column in table.columns]
    raw = table.data
    raw isa AbstractMatrix ||
        throw(ArgumentError("restart data must be a matrix with grid points along rows"))
    size(raw, 2) == length(columns) ||
        throw(DimensionMismatch("restart data has $(size(raw,2)) columns for $(length(columns)) headers"))
    gas, n = f.gas, f.gas.n_species
    select(names) = findfirst(column -> column in names, columns)
    iz, iT = select(_flame_csv_columns.grid), select(_flame_csv_columns.temperature)
    iu = select(_flame_csv_columns.velocity)
    iz === nothing && throw(ArgumentError("restart data lacks a grid column (z_m or grid)"))
    iT === nothing && throw(ArgumentError("restart data lacks a temperature column (temperature_K or T)"))
    f isa FreeFlame && iu === nothing &&
        throw(ArgumentError("restart data lacks a velocity column (velocity_m_s or velocity)"))
    z = Float64.(vec(raw[:, iz]))
    N = length(z)
    N >= 5 && all(isfinite, z) && all(>(0), diff(z)) ||
        throw(ArgumentError("restart grid needs at least five strictly increasing coordinates"))
    T = Float64.(vec(raw[:, iT]))
    u = iu === nothing ? zeros(N) : Float64.(vec(raw[:, iu]))
    all(isfinite, T) && all(isfinite, u) ||
        throw(ArgumentError("restart data contains non-finite values"))
    species_index = Dict(name => k for (k, name) in enumerate(gas.species_names))
    mole_cols = Tuple{Int,Int}[]
    mass_cols = Tuple{Int,Int}[]
    for colidx in eachindex(columns)
        column = columns[colidx]
        if startswith(column, "X_") || startswith(column, "Y_")
            haskey(species_index, column[3:end]) ||
                throw(ArgumentError("restart data contains species $(column[3:end]) not present in the mechanism"))
            push!(startswith(column, "X_") ? mole_cols : mass_cols,
                (colidx, species_index[column[3:end]]))
        end
    end
    # Cantera's restore prefers mole fractions when both spellings appear.
    cols = !isempty(mole_cols) ? mole_cols : mass_cols
    fractions = zeros(n, N)
    for (colidx, k) in cols
        fractions[k, :] = Float64.(vec(raw[:, colidx]))
    end
    # Cantera's SolutionArray restore renormalizes the supplied fractions per
    # grid point; species absent from the table stay zero. Without any species
    # columns the inlet composition represents the gas, as a fresh Cantera
    # flame keeps its phase composition.
    for j in 1:N
        total = sum(@view fractions[:, j])
        total > 0 && (fractions[:, j] ./= total)
    end
    Y = fractions
    if !isempty(mole_cols)
        mean_MW = vec(transpose(gas.MW) * fractions)
        Y = fractions .* gas.MW ./ mean_MW'
    end
    if isempty(cols)
        Y = repeat(f.inlet_Y, 1, N)
    end
    # The inlet boundary shifts the whole profile; the flame state derives the
    # density from the shifted temperature, data composition and pressure, as
    # Cantera's flow domain does on restore.
    T = T .+ (f.inlet_temperature - T[1])
    all(t -> isfinite(t) && 200 <= t <= 6000, T) ||
        throw(ArgumentError("restart temperature profile must remain within 200–6000 K after the inlet shift"))
    mdot = if f isa FreeFlame
        inverseMW = vec(transpose(Y) * (1 ./ gas.MW))
        mean_MW = 1 ./ inverseMW
        u .* f.pressure .* mean_MW ./ (R .* T)
    else
        fill(f.mass_flux, N)
    end
    state = zeros(n + 2, N)
    state[1, :] = T ./ 1000
    state[2:n+1, :] = Y
    state[end, :] = mdot
    f.grid, f.state = z, state
    if f isa FreeFlame
        # Mirror Cantera set_initial_guess: prefer a fixed-temperature point
        # close to Tmid (0.75 T0 + 0.25 Teq) without pinning the boundary
        # knot at the inlet as the anchor.
        T0, Teq, Tmid = T[1], T[end], .75 * T[1] + .25 * T[end]
        below = findall(<(Tmid), T)
        i = isempty(below) ? 1 : last(below)
        if Tmid - T[i] < 0.2 * (Tmid - T0)
            anchor, fT = i, T[i]
        elseif T[i+1] - Tmid < 0.2 * (Teq - Tmid)
            anchor, fT = i + 1, T[i+1]
        else
            anchor, fT = i, Tmid
        end
        f.anchor, f.fixed_temperature = clamp(anchor, 2, N - 1), fT
    else
        f.anchor = clamp(f.anchor, 2, N - 1)
    end
    f.converged = false
    return f
end

export density, heat_release_rate, save_flame, restore_flame!, read_flame_csv, set_initial_guess!
