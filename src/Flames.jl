abstract type AbstractPremixedFlame end

# Temporary fixed-temperature stage used by the Cantera-compatible recovery
# ladder. The lock is removed before the energy-enabled solve is attempted.
const FLAME_TEMPERATURE_LOCKS = WeakKeyDict{Any,Vector{Float64}}()

"A freely propagating, adiabatic, planar ideal-gas premixed flame."
mutable struct FreeFlame{G<:Solution} <: AbstractPremixedFlame
    gas::G
    grid::Vector{Float64}
    pressure::Float64
    inlet_temperature::Float64
    inlet_Y::Vector{Float64}
    state::Matrix{Float64}
    anchor::Int
    fixed_temperature::Float64
    dependent_species::Int
    converged::Bool
    transport_model::Symbol
    soret_enabled::Bool
    multicomponent_data::Union{Nothing,MultiTransportData}
    flux_gradient_basis::Symbol
    discretization::Symbol
    workspace::Any
end

"A planar premixed flame with prescribed inlet mass flux."
mutable struct BurnerFlame{G<:Solution} <: AbstractPremixedFlame
    gas::G
    grid::Vector{Float64}
    pressure::Float64
    inlet_temperature::Float64
    inlet_Y::Vector{Float64}
    state::Matrix{Float64}
    anchor::Int
    fixed_temperature::Float64
    dependent_species::Int
    converged::Bool
    mass_flux::Float64
    profile_positions::Vector{Float64}
    profile_temperatures::Vector{Float64}
    imposed_temperature::Vector{Float64}
    profile_grid_policy::Symbol
    transport_model::Symbol
    soret_enabled::Bool
    multicomponent_data::Union{Nothing,MultiTransportData}
    flux_gradient_basis::Symbol
    discretization::Symbol
    workspace::Any
end

"""
    FreeFlame(gas; T=300, P=one_atm, X, width=0.03, grid=nothing)

Construct a native Julia premixed flame with mixture-averaged diffusion using
mole-fraction gradients. The inlet mass flux (and flame speed) is an eigenvalue.
Initialization uses native constant-enthalpy chemical equilibrium. Call `solve!`
to solve the species and energy equations with adaptive grid refinement.
`discretization=:finite_difference` is the default and matches Cantera's
original discretization. `:conservative` uses shared finite-volume species and
total enthalpy fluxes.
"""
function FreeFlame(gas::Solution; T=300.0, P=one_atm, X, width=0.03, grid=nothing,
        transport_model=gas.trans.model==:ionized_gas ? :ionized_gas : :mixture_averaged,
        multicomponent_data=nothing, soret=false,flux_gradient_basis=:mole,
        discretization=:finite_difference)
    discretization = Symbol(discretization)
    discretization in (:finite_difference,:conservative) ||
        throw(ArgumentError("discretization must be :finite_difference or :conservative"))
    isfinite(T) && 200 <= T <= 6000 && isfinite(P) && P > 0 ||
        throw(ArgumentError("finite inlet temperature in 200–6000 K and positive pressure required"))
    if gas.trans.model == :ionized_gas
        Symbol(transport_model)==:ionized_gas && !soret && multicomponent_data===nothing && flux_gradient_basis==:mole ||
            throw(ArgumentError("ionized flames require ionized transport with mole gradients and no Soret or external transport data"))
        return IonizedFlame(gas;kind=:free,T,P,X,width,grid,discretization)
    end
    gas.trans.poly_order == 5 || throw(ArgumentError("regenerate the sidecar to include native transport fits"))
    x = mole_fractions(gas,X)
    eq = equilibrate(gas; T, P, X=x, mode=:HP)
    eq.T > T+10 || throw(ArgumentError("the inlet does not support an exothermic adiabatic flame"))
    z = isnothing(grid) ? width .* [0,.2,.3,.35,.4,.5,.6,.8,1] : Float64.(grid)
    length(z) >= 5 && all(isfinite,z) && all(>(0),diff(z)) ||
        throw(ArgumentError("the grid needs at least five strictly increasing finite coordinates"))
    y = x .* gas.MW ./ dot(x,gas.MW)
    n = gas.n_species
    state = zeros(n+2,length(z))
    rho = P*dot(x,gas.MW)/(R*T)
    for j in eachindex(z)
        fraction = clamp(((z[j]-z[1])/(z[end]-z[1])-.3)/.2, 0, 1)
        state[1,j] = (T + fraction*(eq.T-T))/1000
        state[2:n+1,j] = y + fraction*(eq.Y-y)
        state[end,j] = rho
    end
    anchor = clamp(argmin(abs.(state[1,:] .- (.75*T+.25*eq.T)/1000)),2,length(z)-1)
    flame = FreeFlame(gas,z,Float64(P),Float64(T),y,state,anchor,
        1000*state[1,anchor],argmax(y),false,:mixture_averaged,false,nothing,:mole,discretization,nothing)
    return set_transport!(flame,transport_model; data=multicomponent_data,soret,flux_gradient_basis)
end

"""
    BurnerFlame(gas; mdot, T=300, P=one_atm, X, width=0.03, grid=nothing)

Construct a burner-stabilized premixed flame. `mdot` is the inlet mass flux in
kg/(m² s). The energy equation is enabled unless `set_temperature_profile!`
specifies a measured temperature profile.
The default `discretization=:finite_difference` matches Cantera's original
discretization. Use `:conservative` for finite-volume species and total enthalpy
balances.
"""
function BurnerFlame(gas::Solution; mdot,T=300.0,P=one_atm,X,width=.03,grid=nothing,
        transport_model=gas.trans.model==:ionized_gas ? :ionized_gas : :mixture_averaged,
        multicomponent_data=nothing,soret=false,flux_gradient_basis=:mole,
        discretization=:finite_difference)
    isfinite(mdot) && mdot > 0 || throw(ArgumentError("mass flux must be finite and positive"))
    if gas.trans.model == :ionized_gas
        Symbol(transport_model)==:ionized_gas && !soret && multicomponent_data===nothing && flux_gradient_basis==:mole ||
            throw(ArgumentError("ionized flames require ionized transport with mole gradients and no Soret or external transport data"))
        return IonizedFlame(gas;kind=:burner,T,P,X,width,grid,mdot,discretization)
    end
    initial_grid = isnothing(grid) ? width .* [0,.1,.2,.3,.5,.7,1] : grid
    base = FreeFlame(gas; T,P,X,width,grid=initial_grid,transport_model,multicomponent_data,soret,flux_gradient_basis,discretization)
    Teq = base.state[1,end]
    Yeq = copy(base.state[2:end-1,end])
    transition_length = (base.grid[end]-base.grid[1])/5
    for j in eachindex(base.grid)
        fraction = clamp((base.grid[j]-base.grid[1])/transition_length,0,1)
        base.state[1,j] = T/1000+fraction*(Teq-T/1000)
        base.state[2:end-1,j] = base.inlet_Y+fraction*(Yeq-base.inlet_Y)
        base.state[end,j] = mdot
    end
    return BurnerFlame(base.gas,base.grid,base.pressure,base.inlet_temperature,
        base.inlet_Y,base.state,base.anchor,base.fixed_temperature,
        base.dependent_species,false,Float64(mdot),Float64[],Float64[],Float64[],:full_knots,
        base.transport_model,base.soret_enabled,base.multicomponent_data,base.flux_gradient_basis,base.discretization,nothing)
end

"""
    set_transport!(flame, model; data=flame.multicomponent_data, soret=false)

Select `:mixture_averaged` or `:multicomponent` diffusion. Multicomponent
diffusion and either model's Soret diffusion require `MultiTransportData`
exported for the same mechanism. `soret=true` includes thermal diffusion.
Use `flux_gradient_basis=:mass` or `:mole` for mixture-averaged diffusion.
An existing solution is retained as the initial guess for the next `solve!`.
Enabling Soret on an adaptive prescribed-temperature grid first inserts all
profile knots and interpolates the state, changing its policy to `:full_knots`.
Invalid transport arguments leave the flame unchanged.
"""
function set_transport!(f::AbstractPremixedFlame,model;data=f.multicomponent_data,soret=false,
        flux_gradient_basis=f.flux_gradient_basis)
    model = Symbol(replace(String(model),"-"=>"_"))
    model in (:mixture_averaged,:multicomponent) || throw(ArgumentError("unknown flame transport model"))
    flux_gradient_basis = Symbol(flux_gradient_basis)
    flux_gradient_basis in (:mole,:mass) || throw(ArgumentError("flux gradient basis must be :mole or :mass"))
    soret isa Bool || throw(ArgumentError("soret must be true or false"))
    data === nothing || data isa MultiTransportData ||
        throw(ArgumentError("transport data must be MultiTransportData or nothing"))
    if model == :multicomponent || soret
        data isa MultiTransportData || throw(ArgumentError("multicomponent transport requires MultiTransportData"))
        data.species_names == f.gas.species_names && data.molecular_weights ≈ f.gas.MW ||
            throw(ArgumentError("multicomponent data must match the flame mechanism"))
    end
    if soret && f isa BurnerFlame && f.profile_grid_policy == :adaptive
        # Validate transport first; profile/state transfer then preserves every knot.
        set_temperature_profile!(f,f.profile_positions,f.profile_temperatures;
            relative=false,grid_policy=:full_knots)
    end
    f.transport_model,f.soret_enabled,f.multicomponent_data = model,soret,data
    f.flux_gradient_basis = flux_gradient_basis
    f.converged = false
    return f
end

function _interpolate_profile(positions,values,z)
    j = clamp(searchsortedlast(positions,z),1,length(positions)-1)
    fraction = (z-positions[j])/(positions[j+1]-positions[j])
    return (1-fraction)*values[j]+fraction*values[j+1]
end

"""
    set_temperature_profile!(flame::BurnerFlame, positions, temperatures;
                             relative=true, grid_policy=:full_knots)

Prescribe temperature [K] using piecewise-linear interpolation. Relative positions
span 0–1; absolute positions are in meters. The profile must cover the domain and
match the burner temperature at its inlet. All supplied pairs are copied and
retained when the solution grid changes.

For conservative flames, `grid_policy=:full_knots` inserts every interior profile
knot. `:adaptive` retains the current grid and adds cells during `solve!` wherever
the exact profile-to-cell-chord error exceeds 1% of the largest supplied temperature,
in addition to species and spacing refinement. Adaptive profiles require conservative
discretization and no Soret diffusion. Enabling Soret later restores full knots.
Finite-difference flames retain their historical profile sampling behavior.
"""
function set_temperature_profile!(f::BurnerFlame,positions,temperatures;
        relative=true,grid_policy=:full_knots)
    grid_policy in (:full_knots,:adaptive) ||
        throw(ArgumentError("profile grid policy must be :full_knots or :adaptive"))
    relative isa Bool || throw(ArgumentError("relative must be true or false"))
    grid_policy == :adaptive && (!_conservative_flame(f) || f.soret_enabled) &&
        throw(ArgumentError("adaptive profiles require conservative discretization without Soret"))
    z,T = Float64.(positions),Float64.(temperatures)
    z isa AbstractVector && T isa AbstractVector ||
        throw(ArgumentError("profile positions and temperatures must be vectors"))
    length(z) == length(T) && length(z) >= 2 || throw(DimensionMismatch("at least two position/temperature pairs required"))
    all(isfinite,z) && all(>(0),diff(z)) && all(t->isfinite(t)&&200<=t<=6000,T) ||
        throw(ArgumentError("profile requires increasing positions and temperatures in 200–6000 K"))
    if relative
        z = f.grid[1] .+ z .* (f.grid[end]-f.grid[1])
    end
    z[1] <= f.grid[1] && z[end] >= f.grid[end] || throw(ArgumentError("profile must cover the entire domain"))
    isapprox(_interpolate_profile(z,T,f.grid[1]),f.inlet_temperature;atol=1e-6,rtol=0) ||
        throw(ArgumentError("profile inlet must match burner temperature"))
    # Prepare the entire transfer before mutating the flame.
    oldz,oldu = f.grid,f.state
    newz = _conservative_flame(f) && grid_policy == :full_knots ?
        sort!(unique(vcat(oldz,filter(x->oldz[1]<x<oldz[end],z)))) : oldz
    newu = length(newz) == length(oldz) ? copy(oldu) :
        [_interpolate_profile(oldz,@view(oldu[k,:]),x) for k in axes(oldu,1), x in newz]
    anchor = findfirst(==(oldz[f.anchor]),newz)
    imposed = [_interpolate_profile(z,T,x) for x in newz]
    newu[1,:] .= imposed ./ 1000
    f.profile_positions,f.profile_temperatures = z,T
    f.grid,f.state,f.anchor = newz,newu,anchor
    f.imposed_temperature = imposed
    f.profile_grid_policy = grid_policy
    f.converged = false
    return f
end

# The difference of two piecewise-linear functions attains its extrema at knots.
function _profile_chord_defects(grid,positions,temperatures)
    errors=zeros(length(grid)-1)
    for j in eachindex(errors)
        left,right=grid[j],grid[j+1]
        Tl=_interpolate_profile(positions,temperatures,left)
        Tr=_interpolate_profile(positions,temperatures,right)
        lo=searchsortedfirst(positions,left)
        hi=searchsortedlast(positions,right)
        for k in lo:hi
            q=(positions[k]-left)/(right-left)
            errors[j]=max(errors[j],abs(temperatures[k]-((1-q)*Tl+q*Tr)))
        end
    end
    return errors
end

function _mark_profile_chord_defects!(insert,f,spacing)
    f isa BurnerFlame && f.profile_grid_policy == :adaptive || return nothing
    _conservative_flame(f) && !f.soret_enabled ||
        throw(ArgumentError("adaptive profiles require conservative discretization without Soret"))
    errors=_profile_chord_defects(f.grid,f.profile_positions,f.profile_temperatures)
    marks=errors .> .01*maximum(f.profile_temperatures)
    any(marks .& (spacing .< 2e-10)) &&
        error("profile error exceeds its temperature budget at minimum grid spacing")
    insert .|= marks
    return nothing
end

"Temperature [K] at each grid point."
temperature(f::AbstractPremixedFlame) = 1000 .* vec(f.state[1,:])
"Species mass fractions, with species along rows and grid points along columns."
mass_fractions(f::AbstractPremixedFlame) = copy(f.state[2:end-1,:])
"Axial velocity [m/s] at each grid point."
function velocity(f::AbstractPremixedFlame)
    return [f.state[end,j]*R*(1000*f.state[1,j])*
        sum(f.state[k+1,j]/f.gas.MW[k] for k in 1:f.gas.n_species)/f.pressure
        for j in eachindex(f.grid)]
end
flame_speed(f::AbstractPremixedFlame) = velocity(f)[1]

# Property adapters used by opposed-flow flames retain their own discretization.
_conservative_flame(f) = false
_conservative_flame(f::Union{FreeFlame,BurnerFlame}) = f.discretization == :conservative

include("FlameDerivatives.jl")

struct ConservativeFlameWorkspace
    density_diffusion::Matrix{Float64}
    species_flux::Matrix{Float64}
    enthalpy_flux::Vector{Float64}
    enthalpy::Vector{Float64}
    cp::Vector{Float64}
    previous_enthalpy::Vector{Float64}
    previous_species_enthalpy::Vector{Float64}
    derivatives::_ConservativeFlameDerivatives
end
ConservativeFlameWorkspace(n,N,gas) = ConservativeFlameWorkspace(zeros(n,N-1),
    zeros(n,N),zeros(N),zeros(N),zeros(N),zeros(N),zeros(n),_ConservativeFlameDerivatives(gas,N))

mutable struct PremixedStats
    newton_calls::Int
    newton_iterations::Int
    jacobians::Int
    residual_evals::Int
    property_nodes::Int
    time_steps::Int
    time_step_failures::Int
    steady_failures::Int
    regrids::Int
    grids::Vector{Int}
end
PremixedStats()=PremixedStats(0,0,0,0,0,0,0,0,0,Int[])

mutable struct _FlamePropertyNodes
    indices::Vector{Int}
    count::Int
end
Base.IteratorSize(::Type{_FlamePropertyNodes}) = Base.HasLength()
Base.eltype(::Type{_FlamePropertyNodes}) = Int
Base.length(nodes::_FlamePropertyNodes) = nodes.count
@inline Base.iterate(nodes::_FlamePropertyNodes,state=1) =
    state > nodes.count ? nothing : (nodes.indices[state],state+1)

struct FlameWorkspace{K}
    kinetics::K
    transport::TransportWorkspace
    stats::Base.RefValue{PremixedStats}
    X::Matrix{Float64}
    rho::Vector{Float64}
    cp::Matrix{Float64}
    h::Matrix{Float64}
    source::Matrix{Float64}
    source_scratch::Vector{Float64}
    flux::Matrix{Float64}
    conductivity::Vector{Float64}
    diffusion_prefactor::Matrix{Float64}
    C::Vector{Float64}
    entropy::Vector{Float64}
    xmid::Vector{Float64}
    ymid::Vector{Float64}
    band::Matrix{Float64}
    perturbed::Matrix{Float64}
    residual_perturbed::Matrix{Float64}
    steps::Vector{Float64}
    property_nodes::_FlamePropertyNodes
    correction::Vector{Float64}
    base_X::Matrix{Float64}
    base_rho::Vector{Float64}
    base_cp::Matrix{Float64}
    base_h::Matrix{Float64}
    base_source::Matrix{Float64}
    multi_transport::Union{Nothing,MultiTransportWorkspace}
    multi_prefactor::Array{Float64,3}
    thermal_diffusion::Matrix{Float64}
    mixture_thermal::Union{Nothing,MixtureThermalDiffusionWorkspace}
    rate_caches::Vector{_KineticsTemperatureCache}
    conservative::Union{Nothing,ConservativeFlameWorkspace}
    # Fallback for grids where signed species stall at their lower bound:
    # clipped node species with a finite-difference Jacobian.
    clip_species::Base.RefValue{Bool}
    # Cantera-policy persistent Jacobian state (SystemJacobian::m_age and
    # SteadyStateSystem::m_jac_ok), used only by _flame_persistent_jacobian
    # flames: the unfactored steady band, its LU pivots, the age of the last
    # evaluation, whether the stored matrix is usable, the 1/dt of the current
    # factorization, and the evaluation counter behind the growth decision.
    steady_jacobian::Matrix{Float64}
    jacobian_pivots::Base.RefValue{Vector{LinearAlgebra.BlasInt}}
    jacobian_age::Base.RefValue{Int}
    jacobian_valid::Base.RefValue{Bool}
    jacobian_rdt::Base.RefValue{Float64}
    jacobian_evals::Base.RefValue{Int}
    configuration::Tuple
end
function FlameWorkspace(f::AbstractPremixedFlame)
    n, N = f.gas.n_species, length(f.grid)
    return FlameWorkspace(KineticsWorkspace(f.gas.reaction),TransportWorkspace(f.gas),Ref(_premixed_stats(f)),
        zeros(n,N),zeros(N),zeros(n,N),zeros(n,N),zeros(n,N),zeros(n),zeros(n,N-1),zeros(N-1),
        zeros(n,N-1),zeros(n),zeros(n),zeros(n),zeros(n),
        zeros(6*(n+2)-2,(n+2)*N),zeros(n+2,N),zeros(n+2,N),zeros(N),_FlamePropertyNodes(Vector{Int}(undef,N),0),zeros((n+2)*N),
        zeros(n,N),zeros(N),zeros(n,N),zeros(n,N),zeros(n,N),
        f.transport_model == :multicomponent ? MultiTransportWorkspace(f.multicomponent_data) : nothing,
        zeros(n,n,f.transport_model == :multicomponent ? N-1 : 0),zeros(n,N-1),
        f.transport_model == :mixture_averaged && f.soret_enabled ?
            MixtureThermalDiffusionWorkspace(f.multicomponent_data) : nothing,
        [_KineticsTemperatureCache(f.gas.reaction) for _ in 1:N],
        _conservative_flame(f) ? ConservativeFlameWorkspace(n,N,f.gas) : nothing,Ref(false),
        zeros(6*(n+2)-2,(n+2)*N),Ref(LinearAlgebra.BlasInt[]),Ref(100000),Ref(false),Ref(NaN),Ref(0),_flame_workspace_configuration(f))
end

# A workspace contains mutable buffers and belongs to exactly one flame.
# Keep its full structural configuration beside it; compare references directly.
_flame_workspace_configuration(f) = (typeof(f), f.gas, f.gas.n_species, f.transport_model,
    f.soret_enabled, f.multicomponent_data, f.flux_gradient_basis,
    _conservative_flame(f), f.gas.species_names, f.gas.MW)

function _flame_workspace_matches(w::FlameWorkspace, f)
    c = w.configuration
    return c[1] === typeof(f) && c[2] === f.gas && c[3] == f.gas.n_species &&
        c[4] === f.transport_model && c[5] == f.soret_enabled &&
        c[6] === f.multicomponent_data && c[7] === f.flux_gradient_basis &&
        c[8] == _conservative_flame(f) && c[9] === f.gas.species_names && c[10] === f.gas.MW
end

function _take_flame_workspace(f)
    hasproperty(f,:transport_model) || return FlameWorkspace(f)
    pools = f.workspace
    if !(pools isa Vector{FlameWorkspace})
        pools = FlameWorkspace[]
        f.workspace = pools
    end
    n, N = f.gas.n_species,length(f.grid)
    idx = findfirst(w -> _flame_workspace_matches(w,f) && size(w.X)==(n,N),pools)
    idx === nothing || return pools[idx]
    # Keep a small flame-local set of grid sizes for retries and ladder stages.
    filter!(w -> _flame_workspace_matches(w,f),pools)
    w = FlameWorkspace(f)
    push!(pools,w)
    length(pools) > 4 && deleteat!(pools,1)
    return w
end

_return_flame_workspace!(f,w) = nothing
_flame_clipped(w) = w isa FlameWorkspace && w.clip_species[]
# Node species value and derivative weight under the active species form.
@inline _flame_node_species(w,y) = _flame_clipped(w) ? max(y,0.0) : y
@inline _flame_node_weight(w,y) = _flame_clipped(w) ? _flame_species_weight(y) : 1.0

# Keep density evaluation in one place: residuals, public flame density and
# refinement interpolation must use the same equation of state and species
# treatment for a given state.
@inline function _flame_density(f,T,Y;clip_species=false)
    sumY,inverseMW=0.0,0.0
    @inbounds for k in eachindex(f.gas.MW)
        y=clip_species ? max(Y[k],0.0) : Y[k]
        sumY+=y
        inverseMW+=y/f.gas.MW[k]
    end
    density=clip_species ? f.pressure*(sumY/inverseMW)/(R*T) :
        f.pressure/(inverseMW*R*T)
    return (density,inverse_mw=inverseMW)
end

@inline function _flame_midpoint_composition!(ymid,xmid,Yleft,Yright,MW,positive)
    sumY,inverseMW=0.0,0.0
    @inbounds for k in eachindex(MW)
        value=.5*(Yleft[k]+Yright[k])
        y=positive ? max(value,0.0) : value
        ymid[k]=y
        sumY+=y
        inverseMW+=y/MW[k]
    end
    meanMW=positive ? sumY/inverseMW : 1/inverseMW
    if positive
        @inbounds for k in eachindex(MW)
            ymid[k]/=sumY
        end
    end
    @inbounds for k in eachindex(MW)
        xmid[k]=ymid[k]*meanMW/MW[k]
    end
    return meanMW
end

@inline _flame_property_stats(w::FlameWorkspace,f::Union{FreeFlame,BurnerFlame}) = w.stats[]
@inline _flame_property_stats(w,f) = nothing

function _flame_properties!(w,f,u; update_transport=true,nodes=eachindex(f.grid),freeze_source=false,
        viscosity=nothing)
    gas, n, N = f.gas, f.gas.n_species, length(f.grid)
    conservative = _conservative_flame(f)
    stats = _flame_property_stats(w,f)
    @inbounds for j in nodes
        if stats !== nothing
            stats.property_nodes += 1
            PREMIXED_STATS.property_nodes += 1
        end
        T = 1000*u[1,j]
        # Node thermochemistry uses signed mass fractions, as Cantera's flow
        # residual does; only transport inputs at faces are clipped below.
        clip = w.clip_species[]
        density_state=_flame_density(f,T,@view(u[2:n+1,j]);clip_species=clip)
        w.rho[j]=density_state.density
        inverseMW=density_state.inverse_mw
        pRT = f.pressure/(R*T)
        for k in 1:n
            w.X[k,j] = (clip ? max(u[k+1,j],0.0) : u[k+1,j])/gas.MW[k]/inverseMW
            w.C[k] = pRT*w.X[k,j]
        end
        x = @view w.X[:,j]
        cal_h_cp_s0!(@view(w.h[:,j]), @view(w.cp[:,j]), w.entropy, gas, T, f.pressure, x)
        @view(w.h[:,j]) .*= R*T
        @view(w.cp[:,j]) .*= R
        w.entropy .*= R
        # A frozen source keeps the base chemistry for hybrid Jacobian columns.
        if !freeze_source
            wdot!(w.source_scratch,gas.reaction,T,w.C,w.entropy,@view(w.h[:,j]),w.kinetics;
                temperature_cache=w.rate_caches[j],pressure=f.pressure)
            for k in 1:n
                w.source[k,j] = w.source_scratch[k]
            end
        end
    end
    @inbounds for j in 1:N-1
        Tmid = 500*(u[1,j]+u[1,j+1])
        meanMW = _flame_midpoint_composition!(w.ymid,w.xmid,
            @view(u[2:n+1,j]),@view(u[2:n+1,j+1]),gas.MW,
            conservative || f.transport_model == :multicomponent || f.soret_enabled)
        if update_transport
            if f.transport_model == :multicomponent
                w.conductivity[j] = multicomponent_transport!(w.multi_transport,
                    f.multicomponent_data,gas,f.pressure,Tmid,w.xmid)
                for k in 1:n, l in 1:n
                    w.multi_prefactor[k,l,j] = f.pressure/(R*Tmid*meanMW)*
                        gas.MW[k]*gas.MW[l]*w.multi_transport.diffusion[k,l]
                end
                # A row-wise constant vanishes against normalized mole-fraction
                # gradients. Center on the fixed dependent species to avoid
                # cancellation in both the primal and analytic contractions.
                for k in 1:n
                    common = w.multi_prefactor[k,f.dependent_species,j]
                    for l in 1:n
                        w.multi_prefactor[k,l,j] -= common
                    end
                end
                w.thermal_diffusion[:,j] .= w.multi_transport.thermal_diffusion
                if conservative
                    # These positive scalar mixture diffusivities determine only
                    # the advective interpolation. Physical diffusion remains
                    # the full multicomponent matrix below.
                    _mixture_transport!(w.transport,gas,f.pressure,Tmid,w.xmid,:mole)
                    for k in 1:n
                        w.conservative.density_diffusion[k,j] =
                            f.pressure*meanMW/(R*Tmid)*w.transport.diffusion[k]
                    end
                end
            else
                mu, w.conductivity[j] = _mixture_transport!(w.transport,gas,f.pressure,Tmid,w.xmid,
                    f.flux_gradient_basis)
                viscosity === nothing || (viscosity[j] = mu)
                pRTm = f.pressure/(R*Tmid)
                pmRTm = f.pressure*meanMW/(R*Tmid)
                for k in 1:n
                    weight = f.flux_gradient_basis == :mass ? meanMW : gas.MW[k]
                    w.diffusion_prefactor[k,j] = pRTm*weight*w.transport.diffusion[k]
                    if conservative
                        w.conservative.density_diffusion[k,j] =
                            pmRTm*w.transport.diffusion[k]
                    end
                end
                if f.soret_enabled
                    w.thermal_diffusion[:,j] .= mixture_thermal_diffusion!(w.mixture_thermal,
                        f.multicomponent_data,f.pressure,Tmid,w.xmid)
                end
            end
        end
        dz = f.grid[j+1]-f.grid[j]
        if f.transport_model == :multicomponent
            for k in 1:n
                flux = 0.0
                for l in 1:n
                    flux += w.multi_prefactor[k,l,j]*(w.X[l,j+1]-w.X[l,j])/dz
                end
                w.flux[k,j] = flux
            end
        else
            fluxsum = 0.0
            for k in 1:n
                gradient = f.flux_gradient_basis == :mass ?
                    (u[k+1,j+1]-u[k+1,j])/dz : (w.X[k,j+1]-w.X[k,j])/dz
                w.flux[k,j] = -w.diffusion_prefactor[k,j]*gradient
                fluxsum += w.flux[k,j]
            end
            for k in 1:n
                w.flux[k,j] -= (conservative ? w.ymid[k] : u[k+1,j])*fluxsum
            end
        end
        if f.soret_enabled
            for k in 1:n
                w.flux[k,j] -= w.thermal_diffusion[k,j]*1000*(u[1,j+1]-u[1,j])/(dz*Tmid)
            end
        end
    end
end

const _flame_timescale = 1e-4


# Cantera SteadyStateSystem/MultiNewton solver policy for finite-difference
# FreeFlame/BurnerFlame solves: one Jacobian is stored per grid and reused
# across Newton calls (re-evaluated only past age 20 or after a damping
# failure with an old Jacobian, MultiNewton.cpp:223-283), refactored with the
# current 1/dt instead (SteadyStateSystem::initTimeInteg), and a successful
# time step grows dt only when it evaluated no Jacobian
# (SteadyStateSystem.cpp:309-311).
_flame_persistent_jacobian(f) = false
_flame_jacobian_count(w,f) = w isa FlameWorkspace ? w.jacobian_evals[] : _premixed_stats(f).jacobians
_flame_persistent_jacobian(f::Union{FreeFlame,BurnerFlame}) =
    f.discretization == :finite_difference

# Cantera SystemJacobian::updateTransient + SteadyStateSystem::factorizeJacobian:
# the stored steady Jacobian gains the diagonal transient term of the current
# time step and is refactorized without a Jacobian re-evaluation. The mask is
# the rows where flame_residual! adds the transient term: temperature and
# species at interior points, except the temperature row of a
# prescribed-temperature burner.
_flame_recoverable_factorization_error(e) = e isa SingularException ||
    (e isa LinearAlgebra.LAPACKException && e.info > 0)

function _flame_factorize!(f::Union{FreeFlame,BurnerFlame},w::FlameWorkspace,rdt)
    B,N = size(f.state)
    # A failed factorization overwrites the band before LAPACK reports its
    # singular pivot. Invalidate every cache field before touching that band;
    # callers may catch the exception and retry the solve later.
    w.jacobian_valid[] = false
    empty!(w.jacobian_pivots[])
    w.jacobian_rdt[] = NaN
    band = w.band
    copyto!(band,w.steady_jacobian)
    if rdt != 0
        diagonal = 4B-1
        # Both Cantera (EigenSparseJacobian::updateTransient) and the native
        # residual subtract the transient term: J_transient = J_steady - M/dt.
        # Cantera applies the sum-Y constraint only at boundary nodes; every
        # species row is transient at interior nodes in the FD formulation.
        # Prescribed-T burner overrides the temperature row at interior points.
        coefficient = _flame_timescale*rdt
        fixed_temperature = f isa BurnerFlame && !isempty(f.imposed_temperature)
        dep_row = f.dependent_species+1
        @inbounds for j in 2:N-1, k in 1:B-1
            fixed_temperature && k == 1 && continue
            k == dep_row && _conservative_flame(f) && continue
            band[diagonal,(j-1)*B+k] -= coefficient
        end
    end
    _,pivots = LinearAlgebra.LAPACK.gbtrf!(2B-1,2B-1,B*N,band)
    w.jacobian_rdt[] = rdt
    w.jacobian_pivots[] = pivots
    w.jacobian_valid[] = true
    return pivots
end

# Exponential-fit interpolation approaches centered interpolation on resolved
# cells and upwind interpolation in convection-dominated cells. One weight is
# shared by every species and total enthalpy, preserving elemental balances.
_flame_face_centering(Pe) = Pe < 1e-4 ? 1-Pe/6+Pe^3/360 : 1+2/Pe-1/tanh(Pe/2)

_flame_node_mass_flux(f,u,w,j) = u[end,j]

function _conservative_flame_fluxes!(c,f,u,w)
    n,N = f.gas.n_species,length(f.grid)
    MW = f.gas.MW
    @inbounds for j in 1:N
        h,cp = 0.0,0.0
        for k in 1:n
            h += u[k+1,j]*w.h[k,j]/MW[k]
            cp += u[k+1,j]*w.cp[k,j]/MW[k]
        end
        c.enthalpy[j],c.cp[j] = h,cp
    end
    @inbounds for j in 1:N-1
        cpface = 0.0
        density_diffusion = Inf
        for k in 1:n
            cpface += .25*(u[k+1,j]+u[k+1,j+1])*(w.cp[k,j]+w.cp[k,j+1])/MW[k]
            density_diffusion = min(density_diffusion,c.density_diffusion[k,j])
        end
        density_diffusion = min(density_diffusion,w.conductivity[j]/cpface)
        dz = f.grid[j+1]-f.grid[j]
        mdot = .5*(_flame_node_mass_flux(f,u,w,j)+_flame_node_mass_flux(f,u,w,j+1))
        Pe = abs(mdot)*dz/density_diffusion
        right_weight = .5*_flame_face_centering(Pe)
        mdot < 0 && (right_weight = 1-right_weight)
        left_weight = 1-right_weight
        hflux = mdot*(left_weight*c.enthalpy[j]+right_weight*c.enthalpy[j+1])-
            1000*w.conductivity[j]*(u[1,j+1]-u[1,j])/dz
        for k in 1:n
            c.species_flux[k,j] = mdot*(left_weight*u[k+1,j]+right_weight*u[k+1,j+1])+w.flux[k,j]
            hflux += .5*(w.h[k,j]+w.h[k,j+1])*w.flux[k,j]/MW[k]
        end
        c.enthalpy_flux[j] = hflux
    end
    # Natural outflow: zero diffusive/conductive boundary flux. Reactions and
    # accumulation in the final half-cell remain in the finite-volume balance.
    @inbounds for k in 1:n
        c.species_flux[k,N] = _flame_node_mass_flux(f,u,w,N)*u[k+1,N]
    end
    c.enthalpy_flux[N] = _flame_node_mass_flux(f,u,w,N)*c.enthalpy[N]
    return c
end

function _flame_previous_enthalpy!(c,f,previous)
    @inbounds for j in eachindex(f.grid)
        T = 1000*previous[1,j]
        # Ideal-gas species enthalpies are composition independent.
        cal_h_RT!(c.previous_species_enthalpy,f.gas,T,f.pressure,f.inlet_Y)
        h = 0.0
        for k in 1:f.gas.n_species
            h += previous[k+1,j]*c.previous_species_enthalpy[k]*R*T/f.gas.MW[k]
        end
        c.previous_enthalpy[j] = h
    end
    return c.previous_enthalpy
end

function _conservative_flame_residual!(residual,f,u,w;previous=nothing,dt=Inf,
        previous_enthalpy=nothing)
    n,N = f.gas.n_species,length(f.grid)
    c = _conservative_flame_fluxes!(w.conservative,f,u,w)
    if previous !== nothing && previous_enthalpy === nothing
        previous_enthalpy = _flame_previous_enthalpy!(c,f,previous)
    end
    @inbounds for j in 1:N
        if f isa BurnerFlame
            residual[end,j] = j == 1 ? u[end,j]-f.mass_flux : u[end,j]-u[end,j-1]
        elseif j == f.anchor
            residual[end,j] = u[1,j]-f.fixed_temperature/1000
        elseif j < f.anchor
            residual[end,j] = u[end,j+1]-u[end,j]
        else
            residual[end,j] = u[end,j]-u[end,j-1]
        end
        if j == 1
            residual[1,j] = u[1,j]-f.inlet_temperature/1000
            for k in 1:n
                residual[k+1,j] = u[end,j]*f.inlet_Y[k]-c.species_flux[k,j]
            end
        else
            cell = .5*(f.grid[min(j+1,N)]-f.grid[j-1])
            factor = _flame_timescale/w.rho[j]
            for k in 1:n
                residual[k+1,j] = factor*(f.gas.MW[k]*w.source[k,j]-
                    (c.species_flux[k,j]-c.species_flux[k,j-1])/cell)
                if previous !== nothing
                    residual[k+1,j] -= _flame_timescale/dt*(u[k+1,j]-previous[k+1,j])
                end
            end
            residual[1,j] = -factor/(1000*c.cp[j])*
                (c.enthalpy_flux[j]-c.enthalpy_flux[j-1])/cell
            if previous !== nothing
                # Formation enthalpy changes with composition. Storing only T
                # here would suppress chemical heating during pseudo-time steps.
                residual[1,j] -= _flame_timescale/(1000*c.cp[j]*dt)*
                    (c.enthalpy[j]-previous_enthalpy[j])
            end
        end
        if f isa BurnerFlame && !isempty(f.imposed_temperature)
            residual[1,j] = u[1,j]-f.imposed_temperature[j]/1000
        end
        if f isa Union{FreeFlame,BurnerFlame} &&
                !(f isa BurnerFlame && !isempty(f.imposed_temperature))
            temperature_lock = get(FLAME_TEMPERATURE_LOCKS, f, nothing)
            temperature_lock === nothing || (residual[1,j] = u[1,j]-temperature_lock[j])
        end
        total = 0.0
        for k in 2:n+1
            total += u[k,j]
        end
        residual[f.dependent_species+1,j] = total - 1
    end
    return residual
end

"Evaluate the discretized steady species, energy, and mass-flow residual."
function flame_residual!(residual, f::AbstractPremixedFlame, u=f.state, w=FlameWorkspace(f);
        previous=nothing, dt=Inf, update_transport=true,nodes=eachindex(f.grid),
        previous_enthalpy=nothing)
    # Counted at entry, before any work: every call to the public residual
    # evaluator increments residual_evals, mirroring counterflow_residual!.
    w.stats[].residual_evals+=1; PREMIXED_STATS.residual_evals+=1
    n,N = f.gas.n_species,length(f.grid)
    _flame_properties!(w,f,u; update_transport,nodes)
    if _conservative_flame(f)
        return _conservative_flame_residual!(residual,f,u,w;previous,dt,previous_enthalpy)
    end
    z, MW = f.grid, f.gas.MW
    @inbounds for j in 1:N
        if f isa BurnerFlame
            residual[end,j] = j == 1 ? u[end,j]-f.mass_flux : u[end,j]-u[end,j-1]
        elseif j == f.anchor
            residual[end,j] = u[1,j] - f.fixed_temperature/1000
        elseif j < f.anchor
            residual[end,j] = u[end,j+1]-u[end,j]
        else
            residual[end,j] = u[end,j]-u[end,j-1]
        end
        if j == 1
            residual[1,j] = u[1,j]-f.inlet_temperature/1000
            for k in 1:n
                residual[k+1,j] = u[end,j]*(f.inlet_Y[k]-u[k+1,j])-w.flux[k,1]
            end
        elseif j == N
            residual[1,j] = u[1,j]-u[1,j-1]
            for k in 1:n
                residual[k+1,j] = u[k+1,j]-u[k+1,j-1]
            end
        else
            left = z[j]-z[j-1]
            right = z[j+1]-z[j]
            cell = .5*(left+right)
            mdot = u[end,j]
            # Positive axial flow: first-order upwind convection, centered diffusion.
            dTdz = 1000*(u[1,j]-u[1,j-1])/left
            cpmean, chemical, enthalpyflux = 0.0,0.0,0.0
            for k in 1:n
                cpmean += u[k+1,j]*w.cp[k,j]/MW[k]
                chemical += w.h[k,j]*w.source[k,j]
                enthalpyflux += .5*(w.flux[k,j-1]+w.flux[k,j])*
                    (w.h[k,j]-w.h[k,j-1])/(left*MW[k])
                residual[k+1,j] = _flame_timescale/w.rho[j]*(MW[k]*w.source[k,j] -
                    (w.flux[k,j]-w.flux[k,j-1])/cell - mdot*(u[k+1,j]-u[k+1,j-1])/left)
            end
            conduction = 1000*(w.conductivity[j]*(u[1,j+1]-u[1,j])/right -
                w.conductivity[j-1]*(u[1,j]-u[1,j-1])/left)/cell
            residual[1,j] = _flame_timescale/(1000*w.rho[j]*cpmean)*
                (conduction-chemical-enthalpyflux-mdot*cpmean*dTdz)
            if previous !== nothing
                for k in 1:n+1
                    residual[k,j] -= _flame_timescale/dt*(u[k,j]-previous[k,j])
                end
            end
        end
        if f isa BurnerFlame && !isempty(f.imposed_temperature)
            residual[1,j] = u[1,j]-f.imposed_temperature[j]/1000
        end
        if f isa Union{FreeFlame,BurnerFlame} &&
                !(f isa BurnerFlame && !isempty(f.imposed_temperature))
            temperature_lock = get(FLAME_TEMPERATURE_LOCKS, f, nothing)
            temperature_lock === nothing || (residual[1,j] = u[1,j]-temperature_lock[j])
        end
        total = 0.0
        for k in 2:n+1
            total += u[k,j]
        end
        # Cantera replaces the dependent-species equation with sum(Y)-1 at
        # boundaries only. Interior FD rows retain every species balance.
        if j == 1 || j == N
            residual[f.dependent_species+1,j] = total - 1
        end
    end
    return residual
end

# Model-specific state contracts share the same banded Newton engine.
_flame_previous_enthalpy(f,w,previous) = _conservative_flame(f) ?
    _flame_previous_enthalpy!(w.conservative,f,previous) : nothing
_flame_perturbation(f,u,k,j)=1e-7*max(abs(u[k,j]),k==1 ? .1 : 1e-5)
@inline function _flame_cantera_species_perturbation(u,k,j)
    x=u[k,j]
    h=1e-5*abs(x)+1e-10
    return x<0 ? -h : h
end
# Fraction of the forward derivative that a positive finite-difference step of
# `_flame_perturbation` samples on a clipped face-midpoint species value: one
# at nonnegative values, the chord fraction within one step below zero.
@inline function _flame_species_weight(y)
    y>=0 && return 1.0
    h=1e-7*max(abs(y),1e-5)
    return max(y+h,0.0)/h
end
_flame_species_forward(y)=_flame_species_weight(y)>0
_flame_needs_correction_check(f)=false
_flame_step_sequencing(f)=_conservative_flame(f)
_flame_checkpoint_steady(f)=_conservative_flame(f)
_flame_reuse_trial(f)=_conservative_flame(f)
_flame_refine_threshold(f,k)=k==1 ? sqrt(eps(Float64))/1000 : sqrt(eps(Float64))
function _flame_bounds(f,k,B)
    low=k==1 ? .2 : k==B ? 1e-6 : -1e-7
    high=k==1 ? 6.0 : k==B ? 100.0 : 1.00001
    return low,high
end
# IdealGasThermo stores each species maximum temperature in the third Trange column.
# Cantera bounds T by twice the minimum species maximum.
_cantera_flame_temperature_bounds(f) = (0.2,2*minimum(@view(f.gas.thermo.Trange[:,end]))/1000)
_cantera_flame_species_bounds() = (-1e-7,1e5)
_cantera_flame_velocity_bounds() = (-1e20,1e20)
@inline _cantera_flame_flow_state_bounds(f,u,j) =
    _cantera_flame_velocity_bounds() .* _flame_state_density(f,u,j)

# Three grid colors exploit the nearest-neighbor block stencil. Each residual
# evaluation perturbs one component at every third point without overlap.
function _flame_jacobian(f,u,w,r; previous=nothing,dt=Inf,previous_enthalpy=nothing,analytic=true,
        update_transport=false,local_property_updates=true,cantera_omission=false)
    if analytic && !update_transport && _conservative_flame(f) && _conservative_flame_jacobian!(f,u,w,r;previous,dt,previous_enthalpy)
        return w.band
    end
    B,N = size(u)
    band = w.band
    fill!(band,0)
    kl = ku = 2*B-1
    perturbed = w.perturbed
    perturbed .= u
    rp = w.residual_perturbed
    steps = w.steps
    copyto!(w.base_X,w.X)
    copyto!(w.base_rho,w.rho)
    copyto!(w.base_cp,w.cp)
    copyto!(w.base_h,w.h)
    copyto!(w.base_source,w.source)
    if _conservative_flame(f) && previous !== nothing && previous_enthalpy === nothing
        previous_enthalpy = _flame_previous_enthalpy!(w.conservative,f,previous)
    end
    # Burner inlet/outlet residuals include face-flux terms from the adjacent
    # interior node. Keep only those endpoint-row finite differences; the
    # analytic helper fills the interior stencil. FreeFlame follows Cantera
    # auto mode and omits the corresponding endpoint entries.
    for k in 1:B, color in 1:3
        perturbed_any=false
        node_count=0
        for j in color:3:N
            auto_species_column = analytic && !update_transport && f isa Union{FreeFlame,BurnerFlame} &&
                f.transport_model == :mixture_averaged && !f.soret_enabled &&
                f.flux_gradient_basis == :mole && 2 <= k <= B-1 && 2 <= j < N
            if auto_species_column && !(f isa BurnerFlame && (j == 2 || j == N-1))
                continue
            end
            cantera_fd = analytic && !update_transport && f isa Union{FreeFlame,BurnerFlame} &&
                f.transport_model == :mixture_averaged && !f.soret_enabled &&
                f.flux_gradient_basis == :mole &&
                2 <= k <= B-1 && (j == 1 || j == N ||
                    f isa BurnerFlame && (j == 2 || j == N-1))
            steps[j] = cantera_fd ? _flame_cantera_species_perturbation(u,k,j) : _flame_perturbation(f,u,k,j)
            perturbed[k,j] = u[k,j]+steps[j]
            perturbed_any=true
            if local_property_updates && f isa Union{FreeFlame,BurnerFlame} &&
                    (update_transport || k != B) && perturbed[k,j] != u[k,j]
                node_count += 1
                w.property_nodes.indices[node_count] = j
            end
        end
        perturbed_any || continue
        # Transport stays frozen for Newton; adjoint sensitivities update it.
        nodes = if !local_property_updates || !(f isa Union{FreeFlame,BurnerFlame})
            color:3:N
        elseif !update_transport && k == B
            ()
        else
            w.property_nodes.count = node_count
            w.property_nodes
        end
        flame_residual!(rp,f,perturbed,w; previous,dt,previous_enthalpy,update_transport,nodes)
        for j in color:3:N
            auto_species_column = analytic && !update_transport && f isa Union{FreeFlame,BurnerFlame} &&
                f.transport_model == :mixture_averaged && !f.soret_enabled &&
                f.flux_gradient_basis == :mole && 2 <= k <= B-1 && 2 <= j < N
            if auto_species_column && !(f isa BurnerFlame && (j == 2 || j == N-1))
                continue
            end
            column = (j-1)*B+k
            for jj in max(1,j-1):min(N,j+1), kk in 1:B
                cantera_boundary_species_omission = cantera_omission && !update_transport &&
                    f isa BurnerFlame && f.transport_model == :mixture_averaged &&
                    !f.soret_enabled && f.flux_gradient_basis == :mole &&
                    2 <= k <= B-1 && ((j == 2 && jj == 1) || (j == N-1 && jj == N))
                cantera_boundary_species_omission && continue
                if auto_species_column &&
                        !((j == 2 && jj == 1) || (j == N-1 && jj == N))
                    continue
                end
                row = (jj-1)*B+kk
                band[kl+ku+1+row-column,column] = (rp[kk,jj]-r[kk,jj])/steps[j]
            end
            perturbed[k,j] = u[k,j]
        end
        copyto!(w.X,w.base_X)
        copyto!(w.rho,w.base_rho)
        copyto!(w.cp,w.base_cp)
        copyto!(w.h,w.base_h)
        copyto!(w.source,w.base_source)
    end
    if analytic && !update_transport && f isa Union{FreeFlame,BurnerFlame}
        _flame_auto_species_jacobian!(band,f,u,w,r;previous,dt) || begin
            # Unsupported chemistry derivative metadata falls back to full FD.
            fill!(band,0)
            fallback=_flame_jacobian(f,u,w,r;previous,dt,previous_enthalpy,
                analytic=false,update_transport,cantera_omission)
            return fallback
        end
    end
    return band
end

"""Populate node and adjacent face active regions for conservative flame derivatives."""
function _flame_active_region!(region::BitMatrix,u::AbstractMatrix)
    n_species = size(u, 1) - 2
    n_nodes = size(u, 2)
    size(region) == (n_species, 2n_nodes - 1) ||
        throw(DimensionMismatch("active region must have size ($(n_species), $(2n_nodes - 1))"))
    @inbounds for j in 1:n_nodes
        # Signed node properties have no branch; only clipped faces do.
        for k in 1:n_species
            region[k,j] = true
        end
        if j < n_nodes
            face = n_nodes + j
            for k in 1:n_species
                region[k,face] = _flame_species_forward(u[k+1,j] + u[k+1,j+1])
            end
        end
    end
    return region
end

# Frozen outer-state component weights; temperature is stored as T/1000.
"""
    NewtonPolicyOptions(; kwargs...)

Immutable controls for Newton damping, bounds, weighted correction norms,
convergence, iteration caps, and Jacobian refresh policy. The defaults preserve
Arrhenius' current behavior. Pass `cantera_compatible=true` to `solve!` to select
the corresponding Cantera policy values, or pass a `NewtonPolicyOptions` value
to select individual settings.
"""
Base.@kwdef struct NewtonPolicyOptions
    cantera_compatible::Bool = false
    damping_trials::Union{Nothing,Int} = nothing
    damping_factor::Union{Nothing,Float64} = nothing
    damping_acceptance::Union{Nothing,Symbol} = nothing
    step_bounds::Union{Nothing,Symbol} = nothing
    norm_denominator::Union{Nothing,Symbol} = nothing
    flow_weight::Union{Nothing,Symbol} = nothing
    steady_convergence::Union{Nothing,Symbol} = nothing
    premixed_steady_maxiters::Union{Nothing,Int} = nothing
    counterflow_steady_maxiters::Union{Nothing,Int} = nothing
    premixed_counterflow_steady_maxiters::Union{Nothing,Int} = nothing
    premixed_transient_maxiters::Union{Nothing,Int} = nothing
    premixed_policy_transient_maxiters::Union{Nothing,Int} = nothing
    counterflow_transient_maxiters::Union{Nothing,Int} = nothing
    jacobian_age::Union{Nothing,Int} = nothing
    jacobian_reevaluations::Union{Nothing,Int} = nothing
end

function _newton_policy(options, cantera_compatible)
    options === nothing || options isa NewtonPolicyOptions ||
        throw(ArgumentError("newton_options must be a NewtonPolicyOptions"))
    cantera_compatible = cantera_compatible || (options !== nothing && options.cantera_compatible)
    base = options === nothing ? NewtonPolicyOptions() : options
    defaults = (damping_trials=7,damping_factor=sqrt(2.0),
        damping_acceptance=:current,step_bounds=:current,
        norm_denominator=:current,flow_weight=:current,
        steady_convergence=:residual_inf,premixed_steady_maxiters=35,
        counterflow_steady_maxiters=45,premixed_counterflow_steady_maxiters=0,
        premixed_transient_maxiters=20,premixed_policy_transient_maxiters=200,
        counterflow_transient_maxiters=100,jacobian_age=20,
        jacobian_reevaluations=4)
    compatible = cantera_compatible || base.cantera_compatible
    cantera = (damping_trials=7,damping_factor=sqrt(2.0),
        damping_acceptance=:step_or_reduction,step_bounds=:cantera,
        norm_denominator=:cantera,flow_weight=:velocity,
        steady_convergence=:weighted_step,premixed_steady_maxiters=0,
        counterflow_steady_maxiters=0,premixed_counterflow_steady_maxiters=0,
        premixed_transient_maxiters=0,premixed_policy_transient_maxiters=0,
        counterflow_transient_maxiters=0,jacobian_age=20,
        jacobian_reevaluations=4)
    resolved = NamedTuple{keys(defaults)}(map(keys(defaults)) do name
        supplied = getfield(base,name)
        supplied === nothing ? (compatible ? getfield(cantera,name) : getfield(defaults,name)) : supplied
    end)
    base = NewtonPolicyOptions(cantera_compatible=compatible; resolved...)
    base.damping_trials>0 && isfinite(base.damping_factor) && base.damping_factor>1 ||
        throw(ArgumentError("Newton damping requires positive trials and a finite factor greater than one"))
    base.damping_acceptance in (:current,:step_or_reduction,:correction_reduction,:residual_merit) ||
        throw(ArgumentError("unsupported damping acceptance policy"))
    base.step_bounds in (:current,:cantera) || throw(ArgumentError("unsupported Newton step bounds"))
    base.norm_denominator in (:current,:cantera) || throw(ArgumentError("unsupported Newton norm denominator"))
    base.flow_weight in (:current,:velocity) || throw(ArgumentError("unsupported Newton flow weight"))
    base.steady_convergence in (:residual_inf,:weighted_step) ||
        throw(ArgumentError("unsupported Newton steady convergence policy"))
    all(>=(0),(base.premixed_steady_maxiters,base.counterflow_steady_maxiters,
        base.premixed_counterflow_steady_maxiters,base.premixed_transient_maxiters,
        base.premixed_policy_transient_maxiters,base.counterflow_transient_maxiters,
        base.jacobian_age,base.jacobian_reevaluations)) ||
        throw(ArgumentError("Newton iteration caps and Jacobian ages must be nonnegative"))
    return base
end

@inline _cantera_damped_step_converged(accepted, s1) = accepted && s1 <= 1
@inline _rejected_small_step(accepted, s1) = !accepted && s1 <= 1
@inline _damping_failure_fallback(accepted, s1, fallback) = fallback && !_rejected_small_step(accepted, s1)

@inline function _newton_merit_accepts(options,s1,s0,allow_converged)
    options.damping_acceptance === :step_or_reduction && return s1<1 || s1<s0
    options.damping_acceptance === :correction_reduction && return s1<s0
    return allow_converged ? s1<1 || s1<s0 : s1<s0
end

function _flame_correction_weights(f,u,transient,options=NewtonPolicyOptions())
    weights=_flame_correction_weights(u,transient)
    if options.flow_weight === :velocity
        absolute=transient ? 1e-11 : 1e-9
        velocity_sum=0.0
        for j in axes(u,2)
            velocity_sum+=abs(_flame_state_velocity(f,u,j))
        end
        weights[end]=1e-4*velocity_sum/size(u,2)+absolute
    end
    return weights
end

@inline function _flame_state_density(f,u,j)
    n=f.gas.n_species
    inv_molecular_weight=0.0
    @inbounds for k in 1:n
        inv_molecular_weight+=u[k+1,j]/f.gas.MW[k]
    end
    return f.pressure/(R*(1000*u[1,j])*inv_molecular_weight)
end
@inline _flame_state_velocity(f,u,j) = u[end,j]/_flame_state_density(f,u,j)
@inline function _flame_step_velocity(f,u,step,j)
    n=f.gas.n_species
    inv_molecular_weight=0.0
    @inbounds for k in 1:n
        inv_molecular_weight+=(u[k+1,j]+step[k+1,j])/f.gas.MW[k]
    end
    density=f.pressure/(R*(1000*(u[1,j]+step[1,j]))*inv_molecular_weight)
    return (u[end,j]+step[end,j])/density
end

# Frozen outer-state component weights; temperature is stored as T/1000.
function _flame_correction_weights(u,transient)
    absolute=transient ? 1e-11 : 1e-9
    weights=Vector{Float64}(undef,size(u,1))
    for k in axes(u,1)
        weights[k]=1e-4*sum(abs,@view(u[k,:]))/size(u,2)+(k==1 ? absolute/1000 : absolute)
    end
    return weights
end
function _flame_correction_norm(step,weights,denominator=length(step))
    total=0.0
    for j in axes(step,2),k in axes(step,1)
        total += (step[k,j]/weights[k])^2
    end
    return sqrt(total/denominator)
end
function _flame_policy_correction_norm(step,weights,f,options,state=f.state)
    denominator=length(step)
    if options.norm_denominator === :cantera
        denominator=(f.gas.n_species+6)*size(step,2)
    end
    if options.flow_weight === :velocity
        total=0.0
        B=size(step,1)
        for j in axes(step,2),k in 1:B-1
            total+=(step[k,j]/weights[k])^2
        end
        for j in axes(step,2)
            u0=_flame_state_velocity(f,state,j)
            u1=_flame_step_velocity(f,state,step,j)
            total+=((u1-u0)/weights[end])^2
        end
        return sqrt(total/denominator)
    end
    return _flame_correction_norm(step,weights,denominator)
end

function _flame_solve_correction!(correction,band,pivots,residual,bandwidth)
    @inbounds for i in eachindex(correction,residual)
        correction[i] = -residual[i]
    end
    LinearAlgebra.LAPACK.gbtrs!('N',bandwidth,bandwidth,length(residual),band,pivots,correction)
    return correction
end
# Prescribed-temperature systems retain residual merit; correction merit helps
# the coupled species and energy solve.
_flame_correction_enabled(f) = _conservative_flame(f) &&
    !(f isa BurnerFlame && !isempty(f.imposed_temperature))

function _flame_newton!(f,w; previous=nothing,dt=Inf,maxiters=35,tolerance=1e-8,
        require_positive=false,minimum_iterations=0,loglevel=0,analytic=true,iterations=Ref(0),
        iteration_callback=nothing,newton_options=NewtonPolicyOptions(),step_trace=nothing,cantera_omission=false)
    newton_options=_newton_policy(newton_options,newton_options.cantera_compatible)
    unlimited_iterations=maxiters<=0
    _premixed_stats(f).newton_calls+=1; PREMIXED_STATS.newton_calls+=1
    u = f.state
    r,trial,rt = similar(u),similar(u),similar(u)
    bandwidth = 2*size(u,1)-1
    pivots = LinearAlgebra.BlasInt[]
    policy = _flame_persistent_jacobian(f)
    # Cantera keeps the Jacobian age on the solver across calls (MultiNewton's
    # Jacobian is owned by SteadyStateSystem); other flames refresh per call.
    age = policy ? w.jacobian_age : Ref(newton_options.jacobian_age)
    force_new_jacobian = false
    jacobian_reevals = 0
    rdt = 0.0
    if policy
        rdt = previous === nothing ? 0.0 : 1.0/dt
        # Cantera initTimeInteg/setSteadyMode: a time-step change refactorizes
        # the stored Jacobian with the new transient diagonal instead of
        # re-evaluating it.
        if w.jacobian_valid[] && w.jacobian_rdt[] != rdt
            pivots = try
                _flame_factorize!(f,w,rdt)
            catch e
                _flame_recoverable_factorization_error(e) || rethrow()
                return false
            end
        end
        # The policy path always solves against the stored factorization; when
        # this call neither refactorized nor refreshed, copy the stored pivots
        # so the local variable names a valid ipiv buffer.
        if w.jacobian_valid[]
            pivots = w.jacobian_pivots[]
        end
    end
    last_contraction = Inf
    # `previous` is immutable throughout this Newton solve. Cache its enthalpy
    # once per solve, never across successive pseudo-time states.
    previous_enthalpy = previous !== nothing ? _flame_previous_enthalpy(f,w,previous) : nothing
    residual_valid = false
    region_enabled = _conservative_flame(f)
    # Cantera applies its weighted-step merit in every Newton solve
    # (MultiNewton::dampStep), independent of the spatial discretization and
    # of rdt: acceptance is s1 < 1 || s1 < s0 and s1 < 1 converges.
    # Conservative prescribed-temperature systems retain residual merit (4d09446).
    fixed_temperature = f isa BurnerFlame && !isempty(f.imposed_temperature)
    cantera_transient = previous !== nothing && !_flame_correction_enabled(f) &&
        (!fixed_temperature || !_conservative_flame(f))
    policy_merit = (policy && (!fixed_temperature || !_conservative_flame(f))) || newton_options.cantera_compatible
    correction_enabled = _flame_correction_enabled(f) || cantera_transient || policy_merit ||
        newton_options.steady_convergence===:weighted_step ||
        newton_options.damping_acceptance in (:step_or_reduction,:correction_reduction)
    correction_weights = correction_enabled ? _flame_correction_weights(f,u,previous !== nothing,newton_options) : Float64[]
    trial_correction = correction_enabled ? Vector{Float64}(undef,length(u)) : Float64[]
    current_region = falses(size(u,1)-2,2size(u,2)-1)
    factored_region = similar(current_region)
    region_valid = false
    iteration = 0
    while unlimited_iterations || iteration < maxiters
        iteration += 1
        iterations[] = iteration
        _premixed_stats(f).newton_iterations+=1; PREMIXED_STATS.newton_iterations+=1
        if !residual_valid
            flame_residual!(r,f,u,w; previous,dt,previous_enthalpy)
        end
        residual_valid = false
        residual_norm = norm(r,Inf)
        iteration_callback === nothing || iteration_callback()
        # Final-grid polishing also checks the physical species criterion.
        # Intermediate grids retain the existing equation-residual criterion.
        if newton_options.steady_convergence !== :weighted_step && !_flame_needs_correction_check(f) && iteration > minimum_iterations && residual_norm < tolerance && (!require_positive ||
                minimum(@view(u[2:end-1,:])) > -1e-12)
            return true
        end
        region_changed = false
        if region_enabled
            _flame_active_region!(current_region,u)
            region_changed = region_valid && current_region != factored_region
        end
        # Cantera MultiNewton::solve refreshes the Jacobian only past its
        # maximum age (m_maxAge = 20) or on demand after a damping failure.
        refresh = policy ? force_new_jacobian || !w.jacobian_valid[] || age[] > newton_options.jacobian_age :
            age[] >= newton_options.jacobian_age || region_changed
        step = try
            if refresh
                if policy
                    # Cantera OneDim::evalJacobian evaluates the steady
                    # Jacobian (rdt = 0); the transient diagonal enters at
                    # factorization (SystemJacobian::updateTransient).
                    flame_residual!(rt,f,u,w)
                    J = _flame_jacobian(f,u,w,rt; analytic,cantera_omission)
                    copyto!(w.steady_jacobian,J)
                    w.jacobian_evals[] += 1
                    _premixed_stats(f).jacobians+=1; PREMIXED_STATS.jacobians+=1
                    pivots = _flame_factorize!(f,w,rdt)
                else
                    J = _flame_jacobian(f,u,w,r; previous,dt,previous_enthalpy,analytic,cantera_omission)
                    _premixed_stats(f).jacobians+=1; PREMIXED_STATS.jacobians+=1
                    _,pivots = LinearAlgebra.LAPACK.gbtrf!(bandwidth,bandwidth,length(u),J)
                end
                if region_enabled
                    copyto!(factored_region,current_region)
                    region_valid = true
                end
                age[] = 0
                force_new_jacobian = false
            end
            correction = _flame_solve_correction!(w.correction,w.band,pivots,r,bandwidth)
            reshape(correction,size(u))
        catch e
            _flame_recoverable_factorization_error(e) || rethrow()
            return false
        end
        all(isfinite,step) || return false
        # MultiNewton::solve increments the Jacobian age after every step.
        policy && (age[] += 1)
        if newton_options.steady_convergence !== :weighted_step && _flame_needs_correction_check(f) && iteration > minimum_iterations &&
                residual_norm < tolerance && _flame_policy_correction_norm(step,correction_weights,f,newton_options) < 1
            return true
        end
        # The strict final-grid tolerance can lie below the residual's round-off
        # floor. If Newton stalls there, accept Cantera's criterion: a fresh
        # Newton correction below one in the steady weighted norm.
        at_floor = newton_options.steady_convergence !== :weighted_step && previous === nothing && refresh && tolerance <= 1e-8 && !require_positive &&
            iteration > minimum_iterations && residual_norm < 10tolerance &&
            _flame_policy_correction_norm(step,_flame_correction_weights(f,u,false,newton_options),f,newton_options) < 1
        alpha = 1.0
        for j in axes(u,2), k in axes(u,1)
            # Small negative species values are permitted in Newton iterates;
            # thermodynamic and kinetic evaluations use the nonnegative state.
            low,high = _flame_bounds(f,k,size(u,1))
            if newton_options.step_bounds === :cantera
                low,high = k==1 ? _cantera_flame_temperature_bounds(f) :
                    k==size(u,1) ? _cantera_flame_flow_state_bounds(f,u,j) :
                    _cantera_flame_species_bounds()
            end
            a = Inf
            if step[k,j] < 0
                a = (u[k,j]-low)/(-step[k,j])
            elseif step[k,j] > 0
                a = (high-u[k,j])/step[k,j]
            end
            alpha = min(alpha,a)
        end
        # Compute the current correction merit before the bound rejection path
        # so its trace and failure fallback have a defined s0 as well.
        step_merit = correction_enabled ? _flame_policy_correction_norm(step,correction_weights,f,newton_options) : 0.0
        # Cantera MultiNewton::boundStep lands exactly on its bounds, and
        # dampStep rejects a bound factor no greater than 1e-10.
        if alpha <= 1e-10
            step_trace === nothing || push!(step_trace,(false,step_merit,Inf,false))
            _damping_failure_fallback(false,step_merit,at_floor) && return true
            if policy
                # MultiNewton::solve: a damping failure re-evaluates an old
                # Jacobian (at most four times) and fails with a fresh one.
                if age[] > 1 && jacobian_reevals < newton_options.jacobian_reevaluations
                    jacobian_reevals += 1
                    force_new_jacobian = true
                    continue
                end
                return false
            end
            refresh && return false
            age[] = 20
            continue
        end
        # Trial correction merit uses the already factored (lagged) Newton J.
        accepted = false
        trial_merit = Inf
        # Cantera's MultiNewton::dampStep: at most seven trials, each sqrt(2)
        # shorter, before the Jacobian is refreshed or the step fails.
        for backtrack in 1:newton_options.damping_trials
            @. trial = u + alpha*step
            flame_residual!(rt,f,trial,w; previous,dt,previous_enthalpy)
            trial_contracts = false
            if all(isfinite,rt)
                if correction_enabled
                    copyto!(trial_correction,vec(rt))
                    LinearAlgebra.LAPACK.gbtrs!('N',bandwidth,bandwidth,length(u),w.band,pivots,trial_correction)
                    trial_merit = _flame_policy_correction_norm(reshape(trial_correction,size(u)),correction_weights,f,newton_options)
                    # Cantera dampStep accepts s1 < 1 or s1 < s0 and treats
                    # s1 < 1 as convergence of the Newton solve.
                    trial_contracts = if newton_options.damping_acceptance === :residual_merit
                        norm(rt) < norm(r)*(1-1e-4*alpha)
                    else
                        _newton_merit_accepts(newton_options,trial_merit,step_merit,
                            cantera_transient || policy_merit)
                    end
                else
                    trial_contracts = norm(rt) < norm(r)*(1-1e-4*alpha)
                end
            end
            if trial_contracts
                u .= trial
                accepted = true
                break
            end
            alpha /= newton_options.damping_factor
        end
        cantera_converged = _cantera_damped_step_converged(accepted,trial_merit)
        step_trace === nothing || push!(step_trace,(accepted,step_merit,trial_merit,cantera_converged))
        loglevel > 1 && println("Newton ",iteration," residual=",residual_norm," damping=",alpha)
        if accepted && ((newton_options.steady_convergence === :weighted_step && cantera_converged) ||
                (newton_options.steady_convergence !== :weighted_step && cantera_transient && trial_merit <= 1))
            return true
        end
        if !accepted || alpha <= 1e-10
            _damping_failure_fallback(false,step_merit,at_floor) && return true
            if policy
                if age[] > 1 && jacobian_reevals < newton_options.jacobian_reevaluations
                    jacobian_reevals += 1
                    force_new_jacobian = true
                    continue
                end
                return false
            end
            refresh && return false
            age[] = 20
            continue
        end
        at_floor && alpha < 1e-3 && return true
        last_contraction = norm(rt)/norm(r)
        # The accepted trial already evaluated the full residual and properties
        # at exactly the newly accepted state. A rejected search never reuses it.
        r,rt = rt,r
        residual_valid = _flame_reuse_trial(f)
        policy || (age[] += 1)
    end
    return false
end

const _timestep_growth_strategies = (:fixed_growth,:steady_norm,:transient_residual,:residual_ratio,:newton_iterations)

function _check_timestep_controls(growth,strategy,regrid,jacobian_mode)
    isfinite(growth) && growth >= 1 || throw(ArgumentError("time step growth factor must be finite and at least 1"))
    strategy in _timestep_growth_strategies || throw(ArgumentError("unknown time step growth strategy $strategy"))
    regrid isa Integer && regrid >= 0 || throw(ArgumentError("time_step_regrid must be a nonnegative integer"))
    jacobian_mode in (:auto,:finite_difference,:analytic) || throw(ArgumentError("unknown jacobian_mode $jacobian_mode"))
    return nothing
end

# Cantera's SteadyStateSystem::calculateTimeStepGrowthFactor after a successful
# step from `before` to `after`; `residual!(r,x,previous,dt)` evaluates in `r`.
function _timestep_growth(residual!,r,strategy,growth,before,after,dt,iterations)
    strategy === :fixed_growth && return growth
    strategy === :newton_iterations && return iterations <= 3 ? growth : 1.0
    steady = strategy === :steady_norm
    residual!(r,before,steady ? nothing : before,steady ? Inf : dt); a = norm(r,Inf)
    residual!(r,after,steady ? nothing : before,steady ? Inf : dt); b = norm(r,Inf)
    strategy === :residual_ratio || return b < a ? growth : 1.0
    b > 0 && a > b || return 1.0
    return min(growth,max(1.0,(a/b)^.2))
end

function _flame_steady!(f; loglevel=0,max_time_steps=500,timestep=Ref(1e-6),tolerance=1e-8,
        clip_species=false,growth_factor=1.5,growth_strategy=:fixed_growth,step_callback=nothing,
        analytic=true,steady_callback=nothing,cantera_omission=false,newton_options=NewtonPolicyOptions())
    newton_options=_newton_policy(newton_options,newton_options.cantera_compatible)
    w = _take_flame_workspace(f)
    try
        return _flame_steady_with_workspace!(f,w;loglevel,max_time_steps,timestep,tolerance,
            clip_species,growth_factor,growth_strategy,step_callback,analytic,steady_callback,cantera_omission,newton_options)

    finally
        _return_flame_workspace!(f,w)
    end
end

function _flame_steady_with_workspace!(f,w; loglevel=0,max_time_steps=500,timestep=Ref(1e-6),tolerance=1e-8,
        clip_species=false,growth_factor=1.5,growth_strategy=:fixed_growth,step_callback=nothing,
        analytic=true,steady_callback=nothing,cantera_omission=false,newton_options=NewtonPolicyOptions())
    newton_options=_newton_policy(newton_options,newton_options.cantera_compatible)
    if w isa FlameWorkspace
        w.stats[] = _premixed_stats(f)
        w.jacobian_valid[] = false
        w.jacobian_age[] = 100000
        w.jacobian_rdt[] = NaN
        w.jacobian_evals[] = 0
        empty!(w.jacobian_pivots[])
    end
    push!(_premixed_stats(f).grids,length(f.grid)); push!(PREMIXED_STATS.grids,length(f.grid))
    # Species reset and clipping apply to neutral premixed flames only; ionized
    # flames keep signed charged species.
    neutral = w isa FlameWorkspace
    clip_species && !neutral && return false
    neutral && (w.clip_species[] = clip_species)
    policy = _flame_persistent_jacobian(f)
    # Cantera time-steps from the state before the failed steady attempt
    # (MultiNewton::solve leaves its input untouched on failure).
    steady_state = _flame_checkpoint_steady(f) || policy ? copy(f.state) : nothing
    steady_ok = _flame_newton!(f,w; loglevel,tolerance,analytic,cantera_omission,newton_options,
        maxiters=newton_options.premixed_steady_maxiters,iteration_callback=steady_callback)
    steady_ok && return true
    _premixed_stats(f).steady_failures+=1; PREMIXED_STATS.steady_failures+=1
    # A failed steady trial must not replace the accepted pseudo-time state.
    steady_state === nothing || copyto!(f.state,steady_state)
    policy && return _flame_time_stepping!(f,w,steady_state; loglevel,max_time_steps,
        timestep,tolerance,growth_factor,growth_strategy,step_callback,analytic,steady_callback,cantera_omission,newton_options)
    dt = timestep[]
    previous = similar(f.state)
    failures = 0
    iterations = Ref(0)
    growth_residual = similar(f.state)
    residual_at!(r,x,previous,dt) = flame_residual!(r,f,x,w; previous,dt)
    for step in 1:max_time_steps
        copyto!(previous,f.state)
        if _flame_newton!(f,w; previous,dt,maxiters=newton_options.premixed_transient_maxiters,
                loglevel=0,analytic,cantera_omission,iterations,newton_options)
            _premixed_stats(f).time_steps+=1; PREMIXED_STATS.time_steps+=1
            failures = 0
            dt = min(dt*_timestep_growth(residual_at!,growth_residual,growth_strategy,growth_factor,
                previous,f.state,dt,iterations[]),1.0)
            timestep[] = dt
            step_callback === nothing || step_callback(dt)
            if step % 10 == 0
                loglevel > 0 && println("Transient step ",step," dt=",dt)
                steady_state === nothing || copyto!(steady_state,f.state)
                trial_ok = _flame_newton!(f,w; loglevel,tolerance,analytic,cantera_omission,newton_options,
                    iteration_callback=steady_callback)
                trial_ok && return true
                steady_state === nothing || copyto!(f.state,steady_state)
            end
        else
            copyto!(f.state,previous)
            failures += 1
            _premixed_stats(f).time_step_failures+=1; PREMIXED_STATS.time_step_failures+=1
            if neutral && failures > 2 && _flame_reset_species!(f.state)
                # As Cantera's resetBadValues: signed species pinned at their
                # lower bound stall every step; clip and renormalize instead.
                failures = 0
                continue
            end
            dt *= .25
            timestep[] = dt
            dt > 1e-12 || return false
        end
    end
    return false
end

# Cantera SteadyStateSystem::solve/timeStep for _flame_persistent_jacobian
# flames: batches of ten successful time steps between steady attempts
# (m_steps = {10}), dt grown only after a step that evaluated no Jacobian
# (SteadyStateSystem.cpp:309-311), halved on failure (m_tfactor = 0.5, minimum
# 1e-16, maximum 1e8), with clipped-species resets after three successive
# failures (resetBadValues). A failed steady attempt restarts time stepping
# from the state before the attempt.
function _flame_time_stepping!(f,w,steady_state; loglevel,max_time_steps,timestep,
        tolerance,growth_factor,growth_strategy,step_callback,analytic,steady_callback=nothing,
        cantera_omission=false,newton_options=NewtonPolicyOptions())
    dt = timestep[]
    previous = similar(f.state)
    iterations = Ref(0)
    growth_residual = similar(f.state)
    residual_at!(r,x,previous,dt) = flame_residual!(r,f,x,w; previous,dt)
    successive_failures = 0
    successes = 0
    batch = 0
    while successes < max_time_steps
        copyto!(previous,f.state)
        jacobians_before = _flame_jacobian_count(w,f)
        # Cantera MultiNewton::solve has no iteration cap: a lagged Jacobian
        # stall is broken by the age-20 refresh, not by giving up the step.
        # 200 is a safety valve only; age refreshes every 21 iterations exit
        # via convergence or damping failure long before it.
        if _flame_newton!(f,w; previous,dt,maxiters=newton_options.premixed_policy_transient_maxiters,
                loglevel=0,analytic,cantera_omission,iterations,newton_options,iteration_callback=steady_callback)
            successive_failures = 0
            successes += 1
            _premixed_stats(f).time_steps+=1; PREMIXED_STATS.time_steps+=1
            batch += 1
            # A step that refreshed the Jacobian does not grow dt.
            if _flame_jacobian_count(w,f) == jacobians_before
                dt *= _timestep_growth(residual_at!,growth_residual,growth_strategy,growth_factor,
                    previous,f.state,dt,iterations[])
            end
            step_callback === nothing || step_callback(dt)
            dt = min(dt,1e8)
            timestep[] = dt
            if batch >= 10
                batch = 0
                loglevel > 0 && println("Transient step ",successes," dt=",dt)
                copyto!(steady_state,f.state)
                trial_ok = _flame_newton!(f,w; loglevel,tolerance,analytic,cantera_omission,newton_options,
                    iteration_callback=steady_callback)
                trial_ok && return true
                copyto!(f.state,steady_state)
            end
        else
            copyto!(f.state,previous)
            successive_failures += 1
            if successive_failures > 2
                # Cantera SteadyStateSystem::timeStep: every third successive
                # failure resets bad values (clip negative species and
                # renormalize) instead of reducing dt.
                _flame_reset_species!(f.state)
                successive_failures = 0
                continue
            end
            dt *= .5
            timestep[] = dt
            dt > 1e-16 || return false
        end
    end
    return false
end

# Clip negative mass fractions and renormalize each node; true if any changed.
function _flame_reset_species!(u)
    changed = false
    @inbounds for j in axes(u,2)
        total = 0.0
        for k in 2:size(u,1)-1
            if u[k,j] < 0
                u[k,j] = 0.0; changed = true
            end
            total += u[k,j]
        end
        total > 0 && (@views u[2:end-1,j] ./= total)
    end
    return changed
end

# Edge temperature gradients relative to the mean gradient (Cantera's check_width).
function _flame_domain_too_narrow(f)
    T = temperature(f); x = f.grid; N = length(x)
    N < 4 && return false
    reference = (T[end]-T[1])/(x[end]-x[1])
    reference == 0 && return false
    left = (T[2]-T[1])/(x[2]-x[1])/reference
    right = (T[N-2]-T[N])/(x[N-2]-x[N])/reference
    return left > 0.02 || right > 0.02
end

_flame_refine_active(f,k) = true
_flame_refine_protect!(keep,f) = nothing

function _refine_flame!(f; ratio=3.,slope=.06,curve=.12,prune=0.,max_points=1000)
    isfinite(prune) && prune <= min(slope,curve) ||
        throw(ArgumentError("prune must be finite and no greater than slope or curve"))
    u,z = f.state,f.grid
    B,N = size(u)
    insert = falses(N-1)
    pruning = prune > 0
    # UNSET points survive unless an active component requests removal. KEEP
    # wins across components; a constant component must not force removal.
    keep = zeros(Int8,N)
    keep[1] = keep[end] = keep[f.anchor] = 1
    spacing = diff(z)
    for k in 1:B
        pruning && !_flame_refine_active(f,k) && continue
        # An imposed piecewise-linear profile has slope discontinuities that
        # cannot be removed by refinement. Its shape is supplied, not solved.
        if k in (1,B) && f isa BurnerFlame && !isempty(f.imposed_temperature)
            continue
        end
        # Refine physical velocity rather than the mass-flux algebraic variable.
        values = k == B ? velocity(f) : @view(u[k,:])
        low,high = extrema(values)
        span = high-low
        threshold = _flame_refine_threshold(f,k)
        if span > .01*max(abs(low),abs(high))
            for j in 1:N-1
                abs(values[j+1]-values[j]) > slope*span+threshold &&
                    spacing[j]>=2e-10 && (insert[j] = true)
                if pruning
                    change = abs(values[j+1]-values[j])/(slope*span+threshold)
                    if change >= prune
                        keep[j] = keep[j+1] = 1
                    elseif keep[j] == 0
                        keep[j] = -1
                    end
                end
            end
        end
        gradients = diff(values) ./ spacing
        glow,ghigh = extrema(gradients)
        gspan = ghigh-glow
        if gspan > .01*max(abs(glow),abs(ghigh))
            for j in 1:N-2
                if _conservative_flame(f) && f isa BurnerFlame && f.soret_enabled &&
                        !isempty(f.imposed_temperature)
                    knot = searchsortedfirst(f.profile_positions,z[j+1])
                    if knot <= length(f.profile_positions) && f.profile_positions[knot] == z[j+1]
                        # Soret couples an imposed T-gradient jump to a real
                        # species-gradient jump. Refinement cannot smooth it.
                        continue
                    end
                end
                if abs(gradients[j+1]-gradients[j]) > curve*gspan+threshold/spacing[j] &&
                        min(spacing[j],spacing[j+1])>=2e-10
                    insert[j] = true
                    insert[j+1] = true
                end
                if pruning
                    change = abs(gradients[j+1]-gradients[j])/(curve*gspan+threshold/spacing[j])
                    if change >= prune
                        keep[j+1] = 1
                    elseif keep[j+1] == 0
                        keep[j+1] = -1
                    end
                end
            end
        end
    end
    for j in 1:N-2
        left,right = z[j+1]-z[j],z[j+2]-z[j+1]
        left > ratio*right && (insert[j] = true)
        right > ratio*left && (insert[j+1] = true)
        if pruning
            left > ratio*right && (keep[max(1,j-1):min(N,j+2)] .= 1)
            right > ratio*left && (keep[j:min(N,j+3)] .= 1)
            # Removing the middle point must not violate either outer ratio.
            j > 1 && z[j+2]-z[j] > ratio*spacing[j-1] && (keep[j+1] = 1)
            j < N-2 && z[j+2]-z[j] > ratio*spacing[j+2] && (keep[j+1] = 1)
        end
    end
    _mark_profile_chord_defects!(insert,f,spacing)
    if pruning
        _flame_refine_protect!(keep,f)
        for j in 3:N-1
            keep[j] == -1 && keep[j-1] == -1 && (keep[j] = 1)
        end
    end
    removed = count(==(-1),keep)
    added = count(j->insert[j] && keep[j] != -1,1:N-1)
    added+removed == 0 && return false
    N-removed+added <= max_points || error("flame refinement exceeds max_points=$max_points")
    # Only free/burner states store mdot in the last row. Ionized and
    # counterflow states use different last-row components and keep direct
    # component interpolation here.
    planar_mdot_state = f isa Union{FreeFlame,BurnerFlame}
    axial_velocity = planar_mdot_state ? velocity(f) : nothing
    newz = Float64[]
    columns = Vector{Float64}[]
    anchor_z = z[f.anchor]
    for j in 1:N
        keep[j] == -1 && continue
        push!(newz,z[j]); push!(columns,u[:,j])
        if j < N && insert[j]
            push!(newz,(z[j]+z[j+1])/2)
            midpoint=(u[:,j]+u[:,j+1])/2
            # Sim1D::refine interpolates the physical velocity component U.
            # The native state stores mdot=rho*u, so rebuild it from the
            # interpolated T/Y and midpoint velocity after averaging.
            if planar_mdot_state
                rho=_flame_density(f,1000*midpoint[1],
                    @view(midpoint[2:f.gas.n_species+1])).density
                midpoint[end]=rho*(axial_velocity[j]+axial_velocity[j+1])/2
            end
            push!(columns,midpoint)
        end
    end
    f.grid = newz
    f.state = reduce(hcat,columns)
    f.anchor = findfirst(==(anchor_z),newz)
    if f isa BurnerFlame && !isempty(f.imposed_temperature)
        f.imposed_temperature = [_interpolate_profile(f.profile_positions,f.profile_temperatures,x) for x in newz]
    end
    f isa AbstractPremixedFlame && (_premixed_stats(f).regrids+=1; PREMIXED_STATS.regrids+=1)
    return true
end

# Cantera's FreeFlame.solve(auto=True) width check (onedim.py check_width):
# thrown after a steady solve whose edge temperature gradients exceed 2% of
# the mean gradient; solve! doubles the domain and restarts the staged ladder.
struct FlameDomainTooNarrow <: Exception end

# FreeFlame's staged initial guess uses the Cantera 30%-50% temperature and
# composition ramp while retaining the inlet mass flux from the preceding
# stage. The compatibility option applies Cantera's fixed-temperature snap or
# Tmid insertion rule on the target grid.
function _fresh_analytic_guess!(f::FreeFlame, z = f.grid; cantera_compatible=true)
    x = f.inlet_Y ./ f.gas.MW
    x ./= sum(x)
    eq = equilibrate(f.gas; T=f.inlet_temperature, P=f.pressure, X=x, mode=:HP)
    T0, y = f.inlet_temperature, f.inlet_Y
    L = z[end] - z[1]
    mdot = f.state[end, 1]
    state = Matrix{Float64}(undef, size(f.state, 1), length(z))
    for j in eachindex(z)
        fraction = clamp(((z[j] - z[1]) / L - .3) / .2, 0, 1)
        state[1, j] = (T0 + fraction * (eq.T - T0)) / 1000
        state[2:size(state, 1) - 1, j] .= y .+ fraction .* (eq.Y .- y)
        state[end, j] = mdot
    end
    grid = collect(Float64, z)
    target = (0.75 * T0 + 0.25 * eq.T) / 1000
    anchor = clamp(argmin(abs.(@view(state[1, :]) .- target)), 2, length(grid) - 1)
    fixed_temperature = 1000 * state[1, anchor]
    if cantera_compatible
        # Cantera FreeFlame.set_initial_guess keeps a nearby grid temperature
        # only within 20% of the corresponding inlet-to-midpoint / midpoint-
        # to-equilibrium interval. Otherwise it inserts the Tmid crossing.
        i = findlast(<(target), @view(state[1, :]))
        if i !== nothing && i < length(grid)
            if target - state[1, i] < 0.2 * (target - T0 / 1000)
                anchor, fixed_temperature = i, 1000 * state[1, i]
            elseif state[1, i + 1] - target < 0.2 * (eq.T / 1000 - target)
                anchor, fixed_temperature = i + 1, 1000 * state[1, i + 1]
            else
                # Sim1D::setFixedTemperature has its own 10%-of-local-gap
                # snap before inserting the crossing point.
                t1, t2 = 1000 * state[1, i], 1000 * state[1, i + 1]
                local_threshold = min(1.0, 0.1 * (t2 - t1))
                if abs(1000 * target - t1) <= local_threshold
                    anchor, fixed_temperature = i, t1
                elseif abs(t2 - 1000 * target) <= local_threshold
                    anchor, fixed_temperature = i + 1, t2
                else
                    fraction = (target - state[1, i]) / (state[1, i + 1] - state[1, i])
                    zfixed = grid[i] + fraction * (grid[i + 1] - grid[i])
                    grid = vcat(grid[1:i], zfixed, grid[i + 1:end])
                    state = hcat(state[:, 1:i],
                        (1 - fraction) .* state[:, i] .+ fraction .* state[:, i + 1],
                        state[:, i + 1:end])
                    anchor, fixed_temperature = i + 1, 1000 * target
                end
            end
        end
    end
    f.grid = grid
    f.state = state
    f.anchor = clamp(anchor, 2, length(grid) - 1)
    f.fixed_temperature = fixed_temperature
    f.converged = false
    return f
end

# Cantera BurnerFlame.set_initial_guess: ramp temperature and composition
# from the burner state to the HP-equilibrium state over the first 20% of the
# domain, while preserving the imposed inlet mass flux.
function _fresh_analytic_guess!(f::BurnerFlame, z = f.grid; cantera_compatible=true)
    x = f.inlet_Y ./ f.gas.MW
    x ./= sum(x)
    eq = equilibrate(f.gas; T=f.inlet_temperature, P=f.pressure, X=x, mode=:HP)
    znew = collect(Float64, z)
    L = znew[end] - znew[1]
    state = Matrix{Float64}(undef, size(f.state, 1), length(znew))
    for j in eachindex(znew)
        fraction = clamp((znew[j] - znew[1]) / (0.2 * L), 0, 1)
        state[1, j] = (f.inlet_temperature + fraction * (eq.T - f.inlet_temperature)) / 1000
        state[2:end-1, j] .= f.inlet_Y .+ fraction .* (eq.Y .- f.inlet_Y)
        state[end, j] = f.mass_flux
    end
    f.grid, f.state = znew, state
    f.anchor = clamp(argmin(abs.(@view(state[1, :]) .- (0.75 * f.inlet_temperature + 0.25 * eq.T) / 1000)),
        2, length(znew) - 1)
    f.fixed_temperature = 1000 * state[1, f.anchor]
    f.converged = false
    return f
end

# Sim1D.cpp's solve(refine_grid=true) loop: steady solve, width check, one
# refinement pass, repeat while refinement changes the grid (`while
# new_points > 0`, entered with `new_points = 1`, so the entry state is always
# solved. An enabled expensive-options final stage re-solves this grid after
# transport is restored. Grids that still guide refinement use the loose
# tolerance. The accepted grid re-solves to the full tolerance (solve!'s sequencing). A
# failed steady solve regrids and retries. The non-compatibility path retains
# its clipped-species retry; the Cantera-compatible path skips that attempt.
function _auto_refine_loop!(f; refine_grid,ratio,slope,curve,prune,max_points,
        max_time_steps,loglevel,domain_check,timestep,controls,initial_time_step,
        ts_regrid=3,cantera_compatible=false)
    strict = !refine_grid || !_flame_step_sequencing(f)
    regrids = 0
    for pass in 1:40
        loglevel > 0 && println("Solving on ", length(f.grid), " points")
        tolerance = strict ? 1e-8 : 1e-6
        # Cantera SteadyStateSystem::solve restarts time stepping from the
        # configured initial step on every grid.
        _flame_persistent_jacobian(f) && (timestep[] = initial_time_step)
        steady_ok = cantera_compatible ?
            _flame_steady!(f; loglevel, max_time_steps, timestep, tolerance, controls...) :
            _flame_steady!(f; loglevel, max_time_steps, timestep, tolerance)
        if !steady_ok && !cantera_compatible
            # The legacy Arrhenius path retries with clipped species.
            loglevel > 0 && println("Retrying ", length(f.grid), " points with clipped species")
            timestep[] = initial_time_step
            steady_ok = _flame_steady!(f; loglevel, max_time_steps, timestep, tolerance,
                clip_species = true, controls...)
        end
        if !steady_ok
            # Cantera's time_step_regrid refines a failed state and retries.
            if !(refine_grid && regrids < ts_regrid &&
                    _refine_flame!(f; ratio, slope, curve, prune, max_points))
                error("flame solve failed on $(length(f.grid)) points")
            end
            regrids += 1
            loglevel > 0 && println("Time stepping failed; regridding to ",
                length(f.grid), " points")
            timestep[] = initial_time_step
            strict = !_flame_step_sequencing(f)
            continue
        end
        domain_check && _flame_domain_too_narrow(f) && throw(FlameDomainTooNarrow())
        if refine_grid && _refine_flame!(f; ratio, slope, curve, prune, max_points)
            strict = !refine_grid || !_flame_step_sequencing(f)
        elseif !strict
            # The grid was accepted by refinement; re-solve to the full
            # tolerance and check refinement once more.
            strict = true
        else
            return f
        end
    end
    error("flame refinement did not converge after 40 passes")
end

struct FlameBlownOff <: Exception end

function _burner_blown_off(f::BurnerFlame)
    n = max(3, length(f.grid) ÷ 5) + 1
    T0, Tend = f.state[1, 1], f.state[1, end]
    rise = Tend - T0
    # Cantera tests the near-inlet rise against the total rise. Treat a
    # machine-flat profile as zero rise before evaluating that ratio.
    abs(rise) <= eps(Float64) * max(abs(T0), abs(Tend), 1.0) && return true
    return abs(f.state[1, n] - f.state[1, 1]) / rise < 1e-6
end

function _burner_blowoff_check!(f::BurnerFlame)
    _burner_blown_off(f) && throw(FlameBlownOff())
    return nothing
end

function _burner_blowoff_solve!(f::BurnerFlame; loglevel,max_time_steps,timestep,controls)
    f.state[1, :] .= f.inlet_temperature / 1000
    f.state[2:end-1, :] .= f.inlet_Y
    f.state[end, :] .= f.mass_flux
    fallback_controls = merge(controls, (; steady_callback=nothing))
    _flame_steady!(f; loglevel,max_time_steps,timestep,tolerance=1e-8,
        fallback_controls...) || error("burner blowoff fallback failed")
    f.converged = true
    return f
end

# Cantera's Sim1D.solve(auto=True) staged ladder: default tolerances and
# Soret disabled on uniform grids of len(grid), 12, 24 and 48 points. Each grid
# gets a fresh kind-specific guess; the first successful fixed-grid solve then
# enters refinement. The expensive-options solve is conditional after restore.
function _auto_staged!(f; refine_grid,ratio,slope,curve,prune,max_points,
        max_time_steps,loglevel,domain_check,timestep,controls,ts_regrid,
        cantera_compatible=true,steady_solve=_flame_steady!)
    initial_step = timestep[]
    soret0 = f.soret_enabled
    transport0 = f.transport_model
    f.soret_enabled = false
    cantera_compatible && transport0 == :multicomponent && (f.transport_model = :mixture_averaged)
    try
        zmin, zmax = f.grid[1], f.grid[end]
        solved = false
        for N in (length(f.grid), 12, 24, 48)
            N > max_points &&
                error("flame refinement exceeds max_points=$max_points")
            # A failed fixed-grid solve must not become the starting state of
            # the next ladder rung. Cantera initializes a fresh guess per N.
            grid_before = copy(f.grid)
            state_before = copy(f.state)
            anchor_before, fixed_temperature_before = f.anchor, f.fixed_temperature
            converged_before = f.converged
            z = N == length(f.grid) ? f.grid : range(zmin, zmax, length = N)
            _fresh_analytic_guess!(f, z; cantera_compatible)
            timestep[] = initial_step
            solved = steady_solve(f; loglevel, max_time_steps, timestep,
                tolerance = 1e-6, controls...)
            if !solved && cantera_compatible
                # Cantera retries a failed energy-enabled fixed grid with
                # temperature frozen, then re-enables energy.
                FLAME_TEMPERATURE_LOCKS[f] = copy(@view f.state[1, :])
                fixed_temperature_solved = try
                    timestep[] = initial_step
                    steady_solve(f; loglevel, max_time_steps, timestep,
                        tolerance = 1e-6, controls...)
                finally
                    delete!(FLAME_TEMPERATURE_LOCKS, f)
                end
                if fixed_temperature_solved
                    timestep[] = initial_step
                    solved = steady_solve(f; loglevel, max_time_steps, timestep,
                        tolerance = 1e-6, controls...)
                end
            end
            if !solved
                f.grid = grid_before
                f.state = state_before
                f.anchor = anchor_before
                f.fixed_temperature = fixed_temperature_before
                f.converged = converged_before
                refine_grid || break
                continue
            end
            domain_check && _flame_domain_too_narrow(f) && throw(FlameDomainTooNarrow())
            if refine_grid
                try
                    _auto_refine_loop!(f; refine_grid, ratio, slope, curve, prune,
                        max_points, max_time_steps, loglevel, domain_check,
                        timestep, controls, initial_time_step = initial_step,
                        ts_regrid = ts_regrid, cantera_compatible)
                catch err
                    # A width failure restarts the domain loop. Other failures
                    # discard this trial state and continue with the next fresh N.
                    (err isa FlameDomainTooNarrow || err isa FlameBlownOff) && rethrow()
                    solved = false
                    continue
                end
            end
            break
        end
        solved || error("Could not find a solution for the 1D problem")
    finally
        f.soret_enabled = soret0
        f.transport_model = transport0
    end
    # Cantera's final solve is needed only when expensive transport options
    # (or user tolerances, which Arrhenius does not expose yet) were requested.
    # The compatibility switch retains Arrhenius's historical unconditional
    # final stage when false.
    # Cantera _onedim.py:1587: `if have_user_tolerances or solve_multi or soret_doms:` with solve_multi and soret_doms
    # saved BEFORE the stages switched them off; transport0/soret0 are those saved flags.
    if !cantera_compatible || transport0 == :multicomponent || soret0
        _auto_refine_loop!(f; refine_grid, ratio, slope, curve, prune, max_points,
            max_time_steps, loglevel, domain_check, timestep, controls,
            initial_time_step = initial_step, ts_regrid = ts_regrid,
            cantera_compatible)
    end
    return f
end

# Cantera checks free-flame domain width after each successful steady solve.
# BurnerFlame uses the shared staged ladder without domain expansion.
function _solve_auto!(f; refine_grid,ratio,slope,curve,prune,max_points,
        max_time_steps,loglevel,domain_check,timestep,controls,ts_regrid,
        cantera_compatible=true)
    f isa Union{FreeFlame,BurnerFlame} || throw(ArgumentError(
        "Cantera-compatible auto staging supports FreeFlame and BurnerFlame"))
    f isa BurnerFlame && !isempty(f.imposed_temperature) && throw(ArgumentError(
        "cantera_compatible auto staging is for energy-enabled burner flames"))
    if f isa BurnerFlame
        domain_check = false
        if cantera_compatible
            controls = merge(controls, (; steady_callback=()->_burner_blowoff_check!(f)))
        end
    end
    initial_step = timestep[]
    expansions = f isa FreeFlame ? 12 : 1
    for _ in 1:expansions
        try
            timestep[] = initial_step
            _auto_staged!(f; refine_grid, ratio, slope, curve, prune,
                max_points, max_time_steps, loglevel, domain_check, timestep,
                controls, ts_regrid = ts_regrid, cantera_compatible)
        catch err
            if err isa FlameBlownOff
                return _burner_blowoff_solve!(f; loglevel,max_time_steps,timestep,controls)
            end
            err isa FlameDomainTooNarrow || rethrow()
            f.grid .*= 2
            loglevel > 0 && println("Expanding domain to ",
                f.grid[end] - f.grid[1], " m")
            refine_grid && _refine_flame!(f; ratio, slope, curve, prune, max_points)
            continue
        end
        while _conservative_flame(f) && minimum(@view(f.state[2:end-1, :])) <= -1e-12
            # A clipped negative trial can converge on an unphysical branch.
            # Build a nonnegative starting guess, preserving normalization,
            # then require an accepted Newton step and a fresh full solve.
            for j in axes(f.state, 2)
                correction = 0.0
                for k in 2:size(f.state, 1) - 1
                    if f.state[k, j] < 0
                        correction -= f.state[k, j]
                        f.state[k, j] = 0.0
                    end
                end
                dominant = argmax(@view(f.state[2:end-1, j])) + 1
                f.state[dominant, j] -= correction
            end
            _flame_newton!(f, FlameWorkspace(f); require_positive = true,
                minimum_iterations = 1, loglevel,
                analytic = controls.analytic,cantera_omission=controls.cantera_omission) ||
                error("flame final species polishing failed on $(length(f.grid)) points")
            refine_grid && _refine_flame!(f; ratio, slope, curve, prune, max_points) ||
                break
            _auto_refine_loop!(f; refine_grid, ratio, slope, curve, prune,
                max_points, max_time_steps, loglevel, domain_check, timestep,
                controls, initial_time_step = initial_step, ts_regrid = ts_regrid)
        end
        f.converged = true
        return f
    end
    error("Cantera-compatible auto solve did not converge")
end

# Cantera's Sim1D solve stats (jacobian/iteration/time-step counts), the
# premixed counterpart of COUNTERFLOW_STATS. Stats live per flame object
# (Cantera keeps solver_stats on the instance); the exported PREMIXED_STATS
# mirrors the accumulated family total for whole-run consumers.
const PREMIXED_STATS=PremixedStats()
# Per-flame stats storage; entries die with their flame (weak keys).
const FLAME_STATS=WeakKeyDict{Any,Any}()
function reset!(s::PremixedStats)
    for name in fieldnames(PremixedStats)
        v=getfield(s,name)
        v isa Vector ? empty!(v) : setfield!(s,name,0)
    end
    return s
end

"""
    solver_stats(flame) -> stats

Return the mutable solve-counter object owned by `flame`, mirroring Cantera's
per-instance `Sim1D.solver_stats`: a `PremixedStats` for `FreeFlame`/
`BurnerFlame` (and ionized premixed flames), a `CounterflowStats` for the
counterflow and impinging-jet kinds. These counters are a diagnostic only;
they never influence solver decisions.

Counted fields per kind (a field a kind never increments stays zero):
- premixed kinds (`FreeFlame`, `BurnerFlame`, `IonizedFlame`): `newton_calls`,
  `newton_iterations`, `jacobians`, `residual_evals`, `property_nodes`, `time_steps`,
  `time_step_failures`, `steady_failures`, `regrids`, `grids`; `property_nodes`
  counts node-property refreshes for `FreeFlame` and `BurnerFlame` and stays zero
  for `IonizedFlame`;
- counterflow kinds (`CounterflowDiffusionFlame`, the premixed-counterflow
  wrappers, impinging and catalytic jets): all of the above plus
  `fd_jacobians`, `chem_jac_failures`, `recoveries` and per-attempt `calls`
  tuples — except catalytic jets, which record `grids` but never `calls`.
`residual_evals` counts every call to the public residual evaluator
(`flame_residual!`/`counterflow_residual!`) at its entry, before argument
validation: damped-step trials, time-step-growth probes and finite-difference
Jacobian columns are all included. Cantera's counterpart excludes the
finite-difference Jacobian evaluations (`OneDim::eval` counts only `count!=0`
calls, and `evalJacobian` passes 0), so native `residual_evals` runs higher
than Cantera's on finite-difference solves.

Exact semantics:
- `solve!` resets the flame's counters at entry and they then report that
  `solve!` only. The reset happens after the keyword-argument validation, so
  an invalid `ratio`/`slope`/`prune`/... argument throws with the previous
  counts intact. Checks that run after the reset point throw with the counters
  already zeroed (counterflow two-point-control validation, catalytic
  transport requirements) or partially counted (runtime `error`s such as a
  `max_points` overrun or `auto` non-convergence).
- Direct calls to the internal helpers (`_flame_newton!`, `_flame_steady!`,
  `_refine_flame!`, `_counterflow_newton!`, `_counterflow_steady!`) accumulate
  on top of the current counters without resetting them.
- A steady attempt that throws is not rolled back: the scalar counters it
  already incremented stay, and its `grids` entry (pushed at attempt entry)
  stays, but counterflow `calls` summaries are appended only when the attempt
  returns, so a thrown attempt gets none.
- Immutable objects driven through the helpers (property adapters, test
  fixtures) count into the family global instead of a per-flame entry.
- The exported `PREMIXED_STATS`/`COUNTERFLOW_STATS` globals mirror the
  accumulated family totals across all flames and are never reset by `solve!`;
  call `Arrhenius.reset!` on them manually for per-run attribution.
"""
solver_stats(f::AbstractPremixedFlame)=_premixed_stats(f)
# WeakKeyDict keys must be mutable; immutable shim objects (property adapters,
# test fixtures driving the helpers directly) share the family global.
_premixed_stats(f)=Base.ismutable(f) ? get!(() -> PremixedStats(), FLAME_STATS, f) : PREMIXED_STATS

"""
    solve!(flame; refine_grid=true, ratio=3, slope=0.06, curve=0.12, prune=0,
           max_points=1000, max_time_steps=500, loglevel=0, auto=false)

Solve a native premixed flame using damped Newton iterations and implicit
pseudo-time stepping. Failure raises an error and leaves `converged=false`.
`cantera_compatible=true` selects the Cantera-compatible Newton policies and staged ladder
for free and burner flames; `auto=true` is an alias. It tries the input grid and uniform
12/24/48-point grids, retries failed-stage solves with temperature locked, then refines
and conditionally repeats after restoring expensive transport. Free flames also expand
a too-narrow domain; burner flames detect blowoff and return a nonreacting solution.
`newton_options` selects individual Newton policy values. Cantera's time-step controls:
`time_step_growth_factor`, `time_step_growth_strategy` (`:fixed_growth`,
`:steady_norm`, `:transient_residual`, `:residual_ratio`, `:newton_iterations`),
`time_step_callback(dt)`, `time_step_regrid` and `jacobian_mode` (`:finite_difference` by default,
`:analytic` to enable residual-faithful analytic species columns, or `:auto` to additionally match Cantera's omitted Burner boundary blocks). `cantera_compatible=true` selects `jacobian_mode=:auto` for this solve; it does not change other defaults. `prune` is Cantera's fourth refine criterion: grid
points whose slope/curve indicators fall below it are removed on each
refinement pass (0 disables pruning; it must not exceed `slope` or `curve`).
"""
function solve!(f::AbstractPremixedFlame; refine_grid=true,ratio=3.,slope=.06,curve=.12,prune=0.,
        max_points=1000,max_time_steps=500,loglevel=0,
        initial_time_step=1e-5,domain_check=f isa FreeFlame,
        time_step_growth_factor=1.5,time_step_growth_strategy=:fixed_growth,
        time_step_callback=nothing,time_step_regrid=3,jacobian_mode=:finite_difference,auto=false,
        newton_options=nothing,cantera_compatible::Bool=auto)
    isfinite(ratio) && ratio > 1 && 0 < slope <= 1 && 0 < curve <= 1 || throw(ArgumentError("invalid refinement criteria"))
    isfinite(prune) && prune <= min(slope,curve) || throw(ArgumentError("prune must be finite and no greater than slope or curve"))
    isfinite(initial_time_step) && initial_time_step > 0 || throw(ArgumentError("initial time step must be finite and positive"))
    cantera_compatible && (jacobian_mode = :auto)
    _check_timestep_controls(time_step_growth_factor,time_step_growth_strategy,time_step_regrid,jacobian_mode)
    newton_options = _newton_policy(newton_options,cantera_compatible)
    reset!(_premixed_stats(f))  # per-flame counters report this solve only
    analytic = jacobian_mode !== :finite_difference
    cantera_compatible = cantera_compatible || auto
    newton_options = _newton_policy(newton_options,cantera_compatible)
    cantera_compatible && (jacobian_mode = :auto)
    analytic = jacobian_mode !== :finite_difference
    controls = (;growth_factor=time_step_growth_factor,growth_strategy=time_step_growth_strategy,
        step_callback=time_step_callback,analytic,cantera_omission=(jacobian_mode === :auto),newton_options)
    timestep = Ref(Float64(initial_time_step))
    (auto || cantera_compatible) && return _solve_auto!(f; refine_grid,ratio,slope,curve,prune,max_points,
        max_time_steps,loglevel,domain_check,timestep,controls,ts_regrid=time_step_regrid,
        cantera_compatible)
    regrids = 0
    f.converged = false
    # Grids that refinement may still change only guide refinement; they use a
    # looser steady tolerance. A grid accepted by refinement is then solved to
    # the full residual tolerance and checked again.
    strict = !refine_grid || !_flame_step_sequencing(f)
    expansions = 0
    for pass in 1:40
        loglevel > 0 && println("Solving on ",length(f.grid)," points")
        tolerance = strict ? 1e-8 : 1e-6
        # Cantera SteadyStateSystem::solve restarts time stepping from the
        # configured initial step on every grid.
        _flame_persistent_jacobian(f) && (timestep[] = initial_time_step)
        if !_flame_steady!(f; loglevel,max_time_steps,timestep,tolerance,controls...)
            # Retry the grid with clipped species; at a nonnegative converged
            # state both forms share the same solution.
            loglevel > 0 && println("Retrying ",length(f.grid)," points with clipped species")
            timestep[] = initial_time_step
            if !_flame_steady!(f; loglevel,max_time_steps,timestep,tolerance,clip_species=true,controls...)
                # Cantera's time_step_regrid: refine the failed state and retry.
                refine_grid && regrids < time_step_regrid && _refine_flame!(f;ratio,slope,curve,prune,max_points) ||
                    error("flame solve failed on $(length(f.grid)) points")
                regrids += 1
                loglevel > 0 && println("Time stepping failed; regridding to ",length(f.grid)," points")
                timestep[] = initial_time_step
                strict = !_flame_step_sequencing(f)
                continue
            end
        end
        if domain_check && expansions < 12 && _flame_domain_too_narrow(f)
            # Cantera's FreeFlame.solve(auto=True): double a domain whose edge
            # temperature gradients are significant, refine once and re-solve.
            f.grid .*= 2
            expansions += 1
            loglevel > 0 && println("Expanding domain to ",f.grid[end]-f.grid[1]," m")
            refine_grid && _refine_flame!(f;ratio,slope,curve,prune,max_points)
            strict = !refine_grid || !_flame_step_sequencing(f)
            continue
        end
        refined = refine_grid && _refine_flame!(f;ratio,slope,curve,prune,max_points)
        if !refined && !strict
            strict = true
            continue
        end
        refined && (strict = !_flame_step_sequencing(f))
        if !refined
            if _conservative_flame(f) && minimum(@view(f.state[2:end-1,:])) <= -1e-12
                # A clipped negative trial can converge on an unphysical branch.
                # Build a nonnegative starting guess, preserving normalization,
                # then require an accepted Newton step and a fresh full solve.
                for j in axes(f.state,2)
                    correction = 0.0
                    for k in 2:size(f.state,1)-1
                        if f.state[k,j] < 0
                            correction -= f.state[k,j]
                            f.state[k,j] = 0.0
                        end
                    end
                    dominant = argmax(@view(f.state[2:end-1,j]))+1
                    f.state[dominant,j] -= correction
                end
                _flame_newton!(f,FlameWorkspace(f);require_positive=true,
                    minimum_iterations=1,loglevel,analytic,cantera_omission=(jacobian_mode === :auto),newton_options) ||
                    error("flame final species polishing failed on $(length(f.grid)) points")
                refine_grid && _refine_flame!(f;ratio,slope,curve,prune,max_points) && continue
            end
            f.converged = true
            return f
        end
    end
    error("flame refinement did not converge after 40 passes")
end
export FreeFlame, BurnerFlame, set_transport!, set_temperature_profile!, solve!, temperature, mass_fractions, velocity, flame_speed, flame_residual!, solver_stats, PremixedStats, PREMIXED_STATS, NewtonPolicyOptions
