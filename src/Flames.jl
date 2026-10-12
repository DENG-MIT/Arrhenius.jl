abstract type AbstractPremixedFlame end

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
end

"""
    FreeFlame(gas; T=300, P=one_atm, X, width=0.03, grid=nothing)

Construct a native Julia premixed flame with mixture-averaged diffusion using
mole-fraction gradients. The inlet mass flux (and flame speed) is an eigenvalue.
Initialization uses native constant-enthalpy chemical equilibrium. Call `solve!`
to solve the species and energy equations with adaptive grid refinement.
`discretization=:conservative` uses shared finite-volume species and total
enthalpy fluxes; `:finite_difference` retains the original discretization.
"""
function FreeFlame(gas::Solution; T=300.0, P=one_atm, X, width=0.03, grid=nothing,
        transport_model=gas.trans.model==:ionized_gas ? :ionized_gas : :mixture_averaged,
        multicomponent_data=nothing, soret=false,flux_gradient_basis=:mole,
        discretization=gas.trans.model==:ionized_gas ? :finite_difference : :conservative)
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
        1000*state[1,anchor],argmax(y),false,:mixture_averaged,false,nothing,:mole,discretization)
    return set_transport!(flame,transport_model; data=multicomponent_data,soret,flux_gradient_basis)
end

"""
    BurnerFlame(gas; mdot, T=300, P=one_atm, X, width=0.03, grid=nothing)

Construct a burner-stabilized premixed flame. `mdot` is the inlet mass flux in
kg/(m² s). The energy equation is enabled unless `set_temperature_profile!`
specifies a measured temperature profile.
Use `discretization=:conservative` for finite-volume species and total enthalpy
balances, or `:finite_difference` for the original discretization.
"""
function BurnerFlame(gas::Solution; mdot,T=300.0,P=one_atm,X,width=.03,grid=nothing,
        transport_model=gas.trans.model==:ionized_gas ? :ionized_gas : :mixture_averaged,
        multicomponent_data=nothing,soret=false,flux_gradient_basis=:mole,
        discretization=gas.trans.model==:ionized_gas ? :finite_difference : :conservative)
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
    if isnothing(grid)
        # The burner preheat length is set by conduction versus convection,
        # lambda/(mdot*cp), rather than by an arbitrary fraction of the domain.
        # Resolve this layer even when most of the domain contains products.
        Xeq = Yeq ./ gas.MW
        Xeq ./= sum(Xeq)
        cp = dot(Yeq,cal_cp_R(gas,1000Teq,P,Xeq).*R./gas.MW)
        _,lambda = mixture_transport!(TransportWorkspace(gas),gas,P,1000Teq,Xeq)
        ell = min(4lambda/(mdot*cp),width/16)
        base.grid = unique(vcat(ell .* [0,.5,1,2,4,8,16],Float64(width)))
        base.state = zeros(gas.n_species+2,length(base.grid))
        base.anchor = 2
        transition_length = 2ell
    end
    for j in eachindex(base.grid)
        fraction = clamp((base.grid[j]-base.grid[1])/transition_length,0,1)
        base.state[1,j] = T/1000+fraction*(Teq-T/1000)
        base.state[2:end-1,j] = base.inlet_Y+fraction*(Yeq-base.inlet_Y)
        base.state[end,j] = mdot
    end
    return BurnerFlame(base.gas,base.grid,base.pressure,base.inlet_temperature,
        base.inlet_Y,base.state,base.anchor,base.fixed_temperature,
        base.dependent_species,false,Float64(mdot),Float64[],Float64[],Float64[],:full_knots,
        base.transport_model,base.soret_enabled,base.multicomponent_data,base.flux_gradient_basis,base.discretization)
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

struct FlameWorkspace{K}
    kinetics::K
    transport::TransportWorkspace
    X::Matrix{Float64}
    rho::Vector{Float64}
    cp::Matrix{Float64}
    h::Matrix{Float64}
    source::Matrix{Float64}
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
    multi_transport::Union{Nothing,MultiTransportWorkspace}
    multi_prefactor::Array{Float64,3}
    thermal_diffusion::Matrix{Float64}
    mixture_thermal::Union{Nothing,MixtureThermalDiffusionWorkspace}
    rate_caches::Vector{_KineticsTemperatureCache}
    conservative::Union{Nothing,ConservativeFlameWorkspace}
end
function FlameWorkspace(f::AbstractPremixedFlame)
    n, N = f.gas.n_species, length(f.grid)
    return FlameWorkspace(KineticsWorkspace(f.gas.reaction),TransportWorkspace(f.gas),
        zeros(n,N),zeros(N),zeros(n,N),zeros(n,N),zeros(n,N),zeros(n,N-1),zeros(N-1),
        zeros(n,N-1),zeros(n),zeros(n),zeros(n),zeros(n),
        zeros(6*(n+2)-2,(n+2)*N),zeros(n+2,N),zeros(n+2,N),zeros(N),
        f.transport_model == :multicomponent ? MultiTransportWorkspace(f.multicomponent_data) : nothing,
        zeros(n,n,f.transport_model == :multicomponent ? N-1 : 0),zeros(n,N-1),
        f.transport_model == :mixture_averaged && f.soret_enabled ?
            MixtureThermalDiffusionWorkspace(f.multicomponent_data) : nothing,
        [_KineticsTemperatureCache(f.gas.reaction) for _ in 1:N],
        _conservative_flame(f) ? ConservativeFlameWorkspace(n,N,f.gas) : nothing)
end

function _flame_properties!(w,f,u; update_transport=true,nodes=eachindex(f.grid))
    gas, n, N = f.gas, f.gas.n_species, length(f.grid)
    conservative = _conservative_flame(f)
    @inbounds for j in nodes
        T = 1000*u[1,j]
        sumY = 0.0
        inverseMW = 0.0
        for k in 1:n
            y = max(u[k+1,j],0.0)
            sumY += y
            inverseMW += y/gas.MW[k]
        end
        meanMW = sumY/inverseMW
        w.rho[j] = f.pressure*meanMW/(R*T)
        for k in 1:n
            w.X[k,j] = max(u[k+1,j],0.0)/gas.MW[k]/inverseMW
            w.C[k] = f.pressure/(R*T)*w.X[k,j]
        end
        x = @view w.X[:,j]
        h = @view w.h[:,j]
        cp = @view w.cp[:,j]
        cal_h_RT!(h,gas,T,f.pressure,x)
        cal_cp_R!(cp,gas,T,f.pressure,x)
        cal_s0_R!(w.entropy,gas,T,f.pressure,x)
        h .*= R*T
        cp .*= R
        w.entropy .*= R
        wdot!(@view(w.source[:,j]),gas.reaction,T,w.C,w.entropy,h,w.kinetics;
            temperature_cache=w.rate_caches[j])
    end
    @inbounds for j in 1:N-1
        Tmid = 500*(u[1,j]+u[1,j+1])
        ysum, inverseMW = 0.0, 0.0
        for k in 1:n
            w.ymid[k] = max(0.0, .5*(u[k+1,j]+u[k+1,j+1]))
            ysum += w.ymid[k]
            inverseMW += w.ymid[k]/gas.MW[k]
        end
        meanMW = ysum/inverseMW
        for k in 1:n
            w.ymid[k] /= ysum
            w.xmid[k] = w.ymid[k]*meanMW/gas.MW[k]
        end
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
                    mixture_transport!(w.transport,gas,f.pressure,Tmid,w.xmid)
                    for k in 1:n
                        w.conservative.density_diffusion[k,j] =
                            f.pressure*meanMW/(R*Tmid)*w.transport.diffusion[k]
                    end
                end
            else
                _, w.conductivity[j] = mixture_transport!(w.transport,gas,f.pressure,Tmid,w.xmid;
                    basis=f.flux_gradient_basis)
                for k in 1:n
                    weight = f.flux_gradient_basis == :mass ? meanMW : gas.MW[k]
                    w.diffusion_prefactor[k,j] = f.pressure/(R*Tmid)*weight*w.transport.diffusion[k]
                    if conservative
                        w.conservative.density_diffusion[k,j] =
                            f.pressure*meanMW/(R*Tmid)*w.transport.diffusion[k]
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
        residual[f.dependent_species+1,j] = sum(@view(u[2:n+1,j]))-1
    end
    return residual
end

"Evaluate the discretized steady species, energy, and mass-flow residual."
function flame_residual!(residual, f::AbstractPremixedFlame, u=f.state, w=FlameWorkspace(f);
        previous=nothing, dt=Inf, update_transport=true,nodes=eachindex(f.grid),
        previous_enthalpy=nothing)
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
        residual[f.dependent_species+1,j] = sum(@view(u[2:n+1,j]))-1
    end
    return residual
end

# Model-specific state contracts share the same banded Newton engine.
_flame_previous_enthalpy(f,w,previous) = _conservative_flame(f) ?
    _flame_previous_enthalpy!(w.conservative,f,previous) : nothing
_flame_perturbation(f,u,k,j)=1e-7*max(abs(u[k,j]),k==1 ? .1 : 1e-5)
_flame_correction_weights(f,u,transient)=_flame_correction_weights(u,transient)
_flame_needs_correction_check(f)=false
_flame_checkpoint_steady(f)=_conservative_flame(f)
_flame_reuse_trial(f)=_conservative_flame(f)
_flame_refine_threshold(f,k)=k==1 ? sqrt(eps(Float64))/1000 : sqrt(eps(Float64))
function _flame_bounds(f,k,B)
    low=k==1 ? .2 : k==B ? 1e-6 : -1e-7
    high=k==1 ? 6.0 : k==B ? 100.0 : 1.00001
    return low,high
end

# Three grid colors exploit the nearest-neighbor block stencil. Each residual
# evaluation perturbs one component at every third point without overlap.
function _flame_jacobian(f,u,w,r; previous=nothing,dt=Inf,previous_enthalpy=nothing,analytic=true)
    if analytic && _conservative_flame(f) && _conservative_flame_jacobian!(f,u,w,r;previous,dt,previous_enthalpy)
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
    base_properties = (copy(w.X),copy(w.rho),copy(w.cp),copy(w.h),copy(w.source))
    if _conservative_flame(f) && previous !== nothing && previous_enthalpy === nothing
        previous_enthalpy = _flame_previous_enthalpy!(w.conservative,f,previous)
    end
    for k in 1:B, color in 1:3
        for j in color:3:N
            steps[j] = _flame_perturbation(f,u,k,j)
            perturbed[k,j] = u[k,j]+steps[j]
        end
        flame_residual!(rp,f,perturbed,w; previous,dt,previous_enthalpy,update_transport=false,nodes=color:3:N)
        for j in color:3:N
            column = (j-1)*B+k
            for jj in max(1,j-1):min(N,j+1), kk in 1:B
                row = (jj-1)*B+kk
                band[kl+ku+1+row-column,column] = (rp[kk,jj]-r[kk,jj])/steps[j]
            end
            perturbed[k,j] = u[k,j]
        end
        w.X .= base_properties[1]
        w.rho .= base_properties[2]
        w.cp .= base_properties[3]
        w.h .= base_properties[4]
        w.source .= base_properties[5]
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
        for k in 1:n_species
            region[k,j] = u[k+1,j] >= 0
        end
        if j < n_nodes
            face = n_nodes + j
            for k in 1:n_species
                region[k,face] = u[k+1,j] + u[k+1,j+1] >= 0
            end
        end
    end
    return region
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
function _flame_correction_norm(step,weights)
    total=0.0
    for j in axes(step,2),k in axes(step,1)
        total += (step[k,j]/weights[k])^2
    end
    return sqrt(total/length(step))
end
# Prescribed-temperature systems retain residual merit; correction merit helps
# the coupled species and energy solve.
_flame_correction_enabled(f) = _conservative_flame(f) &&
    !(f isa BurnerFlame && !isempty(f.imposed_temperature))

function _flame_newton!(f,w; previous=nothing,dt=Inf,maxiters=35,tolerance=1e-8,
        require_positive=false,minimum_iterations=0,loglevel=0)
    u = f.state
    r,trial,rt = similar(u),similar(u),similar(u)
    bandwidth = 2*size(u,1)-1
    pivots = LinearAlgebra.BlasInt[]
    age = 5
    last_contraction = Inf
    # `previous` is immutable throughout this Newton solve. Cache its enthalpy
    # once per solve, never across successive pseudo-time states.
    previous_enthalpy = previous !== nothing ? _flame_previous_enthalpy(f,w,previous) : nothing
    residual_valid = false
    region_enabled = _conservative_flame(f)
    correction_enabled = _flame_correction_enabled(f)
    correction_weights = correction_enabled ? _flame_correction_weights(f,u,previous !== nothing) : Float64[]
    trial_correction = correction_enabled ? Vector{Float64}(undef,length(u)) : Float64[]
    current_region = falses(size(u,1)-2,2size(u,2)-1)
    factored_region = similar(current_region)
    region_valid = false
    for iteration in 1:maxiters
        if !residual_valid
            flame_residual!(r,f,u,w; previous,dt,previous_enthalpy)
        end
        residual_valid = false
        residual_norm = norm(r,Inf)
        # Final-grid polishing also checks the physical species criterion.
        # Intermediate grids retain the existing equation-residual criterion.
        if !_flame_needs_correction_check(f) && iteration > minimum_iterations && residual_norm < tolerance && (!require_positive ||
                minimum(@view(u[2:end-1,:])) > -1e-12)
            return true
        end
        region_changed = false
        if region_enabled
            _flame_active_region!(current_region,u)
            region_changed = region_valid && current_region != factored_region
        end
        refresh = age >= 5 || last_contraction > .7 || region_changed
        step = try
            if refresh
                J = _flame_jacobian(f,u,w,r; previous,dt,previous_enthalpy)
                _,pivots = LinearAlgebra.LAPACK.gbtrf!(bandwidth,bandwidth,length(u),J)
                if region_enabled
                    copyto!(factored_region,current_region)
                    region_valid = true
                end
                age = 0
            end
            correction = -vec(copy(r))
            LinearAlgebra.LAPACK.gbtrs!('N',bandwidth,bandwidth,length(u),w.band,pivots,correction)
            reshape(correction,size(u))
        catch e
            e isa SingularException || rethrow()
            return false
        end
        all(isfinite,step) || return false
        if _flame_needs_correction_check(f) && iteration > minimum_iterations &&
                residual_norm < tolerance && _flame_correction_norm(step,correction_weights) < 1
            return true
        end
        alpha = 1.0
        for j in axes(u,2), k in axes(u,1)
            # Small negative species values are permitted in Newton iterates;
            # thermodynamic and kinetic evaluations use the nonnegative state.
            low,high = _flame_bounds(f,k,size(u,1))
            if step[k,j] < 0
                alpha = min(alpha,.99*(u[k,j]-low)/(-step[k,j]))
            elseif step[k,j] > 0
                alpha = min(alpha,.99*(high-u[k,j])/step[k,j])
            end
        end
        # Trial correction merit uses the already factored (lagged) Newton J.
        step_merit = correction_enabled ? _flame_correction_norm(step,correction_weights) : 0.0
        accepted = false
        for backtrack in 1:24
            @. trial = u + alpha*step
            flame_residual!(rt,f,trial,w; previous,dt,previous_enthalpy)
            trial_contracts = false
            if all(isfinite,rt)
                if correction_enabled
                    copyto!(trial_correction,vec(rt))
                    LinearAlgebra.LAPACK.gbtrs!('N',bandwidth,bandwidth,length(u),w.band,pivots,trial_correction)
                    trial_contracts = _flame_correction_norm(reshape(trial_correction,size(u)),correction_weights) < step_merit
                else
                    trial_contracts = norm(rt) < norm(r)*(1-1e-4*alpha)
                end
            end
            if trial_contracts
                u .= trial
                accepted = true
                break
            end
            alpha *= .5
        end
        loglevel > 1 && println("Newton ",iteration," residual=",residual_norm," damping=",alpha)
        if !accepted || alpha <= 1e-10
            refresh && return false
            age = 5
            continue
        end
        last_contraction = norm(rt)/norm(r)
        # The accepted trial already evaluated the full residual and properties
        # at exactly the newly accepted state. A rejected search never reuses it.
        r,rt = rt,r
        residual_valid = _flame_reuse_trial(f)
        age += 1
    end
    return false
end

function _flame_steady!(f; loglevel=0,max_time_steps=500,timestep=Ref(1e-6))
    w = FlameWorkspace(f)
    steady_state = _flame_checkpoint_steady(f) ? copy(f.state) : nothing
    _flame_newton!(f,w; loglevel) && return true
    # A failed steady trial must not replace the accepted pseudo-time state.
    steady_state === nothing || copyto!(f.state,steady_state)
    dt = timestep[]
    previous = similar(f.state)
    for step in 1:max_time_steps
        copyto!(previous,f.state)
        if _flame_newton!(f,w; previous,dt,maxiters=20,loglevel=0)
            dt = min(dt*1.5,1.0)
            timestep[] = dt
            if step % 10 == 0
                loglevel > 0 && println("Transient step ",step," dt=",dt)
                steady_state === nothing || copyto!(steady_state,f.state)
                _flame_newton!(f,w; loglevel) && return true
                steady_state === nothing || copyto!(f.state,steady_state)
            end
        else
            copyto!(f.state,previous)
            dt *= .25
            timestep[] = dt
            dt > 1e-12 || return false
        end
    end
    return false
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
    newz = Float64[]
    columns = Vector{Float64}[]
    anchor_z = z[f.anchor]
    for j in 1:N
        keep[j] == -1 && continue
        push!(newz,z[j]); push!(columns,u[:,j])
        if j < N && insert[j]
            push!(newz,(z[j]+z[j+1])/2)
            push!(columns,(u[:,j]+u[:,j+1])/2)
        end
    end
    f.grid = newz
    f.state = reduce(hcat,columns)
    f.anchor = findfirst(==(anchor_z),newz)
    if f isa BurnerFlame && !isempty(f.imposed_temperature)
        f.imposed_temperature = [_interpolate_profile(f.profile_positions,f.profile_temperatures,x) for x in newz]
        f.state[1,:] .= f.imposed_temperature ./ 1000
    end
    return true
end

"""
    solve!(flame; refine_grid=true, ratio=3, slope=0.06, curve=0.12,
           max_points=1000, max_time_steps=500, loglevel=0)

Solve a native premixed flame using damped Newton iterations and implicit
pseudo-time stepping. Failure raises an error and leaves `converged=false`.
"""
function solve!(f::AbstractPremixedFlame; refine_grid=true,ratio=3.,slope=.06,curve=.12,
        max_points=1000,max_time_steps=500,loglevel=0,
        initial_time_step=f isa FreeFlame ? 1e-6 : 1e-5)
    isfinite(ratio) && ratio > 1 && 0 < slope <= 1 && 0 < curve <= 1 || throw(ArgumentError("invalid refinement criteria"))
    isfinite(initial_time_step) && initial_time_step > 0 || throw(ArgumentError("initial time step must be finite and positive"))
    f.converged = false
    timestep = Ref(Float64(initial_time_step))
    for pass in 1:40
        loglevel > 0 && println("Solving on ",length(f.grid)," points")
        _flame_steady!(f; loglevel,max_time_steps,timestep) || error("flame solve failed on $(length(f.grid)) points")
        if !refine_grid || !_refine_flame!(f;ratio,slope,curve,max_points)
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
                    minimum_iterations=1,loglevel) ||
                    error("flame final species polishing failed on $(length(f.grid)) points")
                refine_grid && _refine_flame!(f;ratio,slope,curve,max_points) && continue
            end
            f.converged = true
            return f
        end
    end
    error("flame refinement did not converge after 40 passes")
end
export FreeFlame, BurnerFlame, set_transport!, set_temperature_profile!, solve!, temperature, mass_fractions, velocity, flame_speed, flame_residual!
