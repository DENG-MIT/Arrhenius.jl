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

export density, heat_release_rate, save_flame, restore_flame!
