"""Planar premixed flame with signed charged species and an algebraic electric field.

Construct through FreeFlame or BurnerFlame with an ionized-gas mechanism.
The initial stage freezes charged diffusion. Enable it with
set_electric_field!(flame, true) before the second solve!.
Conservative flames first initialize with a frozen finite-difference solve.
"""
mutable struct IonizedFlame{G<:Solution} <: AbstractPremixedFlame
    gas::G
    grid::Vector{Float64}
    pressure::Float64
    inlet_temperature::Float64
    inlet_Y::Vector{Float64}
    state::Matrix{Float64}
    anchor::Int
    fixed_temperature::Float64
    converged::Bool
    kind::Symbol
    mass_flux::Float64
    field_enabled::Bool
    ion_data::IonTransportData
    discretization::Symbol
    initialized::Bool
end

const _ION_FLAME_FARADAY = _ION_QE*6.02214076e26
const _ION_FLAME_EPS0 = inv(299792458.0^2 *
    (2*7.2973525693e-3*6.62607015e-34/(_ION_QE^2*299792458.0)))

function IonizedFlame(gas::Solution;kind=:free,T=300.0,P=one_atm,X,
        width=.05,grid=nothing,mdot=nothing,discretization=:finite_difference)
    discretization=Symbol(discretization)
    discretization in (:finite_difference,:conservative) || throw(ArgumentError("invalid ionized discretization"))
    kind in (:free,:burner) || throw(ArgumentError("kind must be :free or :burner"))
    isfinite(T) && 200<=T<=6000 && isfinite(P) && P>0 ||
        throw(ArgumentError("finite temperature in 200-6000 K and positive pressure required"))
    isfinite(width) && width>0 || throw(ArgumentError("positive finite width required"))
    kind==:burner && !(mdot isa Real && isfinite(mdot) && mdot>0) &&
        throw(ArgumentError("burner mass flux must be finite and positive"))
    data=IonTransportData(gas)
    x=mole_fractions(gas,X);y=x.*gas.MW./dot(x,gas.MW)
    eq=equilibrate(gas;T,P,X=x,mode=:HP)
    rho=P*dot(x,gas.MW)/(R*T)
    massflux=kind==:burner ? Float64(mdot) : rho
    z=isnothing(grid) ? width.*[0,.2,.3,.35,.4,.5,.6,.8,1] : Float64.(grid)
    transition=.2*width
    if isnothing(grid) && kind==:burner
        _,lambda,_=ionized_transport!(IonTransportWorkspace(data),data,P,eq.T,eq.X)
        cp=dot(eq.Y,cal_cp_R(gas,eq.T,P,eq.X).*R./gas.MW)
        ell=min(4lambda/(massflux*cp),width/16)
        z=unique(vcat(ell.*[0,.5,1,2,4,8,16],Float64(width)))
        transition=2ell
    end
    length(z)>=5 && all(isfinite,z) && all(>(0),diff(z)) ||
        throw(ArgumentError("at least five increasing finite grid coordinates required"))
    span=z[end]-z[1];n=gas.n_species;u=zeros(n+3,length(z))
    for j in eachindex(z)
        fraction=kind==:free ? clamp(((z[j]-z[1])/span-.3)/.2,0,1) :
            clamp((z[j]-z[1])/transition,0,1)
        temp=T+fraction*(eq.T-T)
        u[1,j]=temp/1000
        for k in 1:n
            u[k+1,j]=y[k]+fraction*(eq.Y[k]-y[k])
        end
        inverseMW=sum(u[k+1,j]/gas.MW[k] for k in 1:n)
        u[end,j]=massflux*R*temp*inverseMW/P
    end
    anchor=clamp(argmin(abs.(u[1,:].-(.75*T+.25*eq.T)/1000)),2,length(z)-1)
    IonizedFlame(gas,z,Float64(P),Float64(T),y,u,anchor,1000*u[1,anchor],
        false,kind,massflux,false,data,discretization,false)
end

"""Enable or freeze charged diffusion and the electric-field equation for the next solve."""
function set_electric_field!(f::IonizedFlame,enabled::Bool)
    f.field_enabled=enabled
    f.converged=false
    return f
end
electric_field(f::IonizedFlame)=1000 .* vec(f.state[end-1,:])
mass_fractions(f::IonizedFlame)=copy(f.state[2:f.gas.n_species+1,:])
velocity(f::IonizedFlame)=copy(vec(f.state[end,:]))

include("IonConservativeFlames.jl")
include("IonFlameProperties.jl")
include("IonFlameResidual.jl")

# Use the same signed, unnormalized thermodynamic convention as the residual.
function density(f::IonizedFlame)
    return f.pressure ./ (R .* temperature(f) .* vec(sum(mass_fractions(f)./f.gas.MW;dims=1)))
end
function heat_release_rate(f::IonizedFlame)
    w=FlameWorkspace(f)
    _ion_flame_properties!(w,f,f.state)
    return -vec(sum(w.h.*w.source;dims=1))
end
function set_transport!(f::IonizedFlame,model;data=nothing,soret=false,flux_gradient_basis=:mole)
    Symbol(replace(String(model),"-"=>"_"))==:ionized_gas &&
        data===nothing && soret===false && flux_gradient_basis==:mole ||
        throw(ArgumentError("ionized flames require ionized-gas mole-gradient transport without Soret"))
    f.converged=false
    return f
end


_flame_node_mass_flux(f::IonizedFlame,u,w,j)=w.rho[j]*u[end,j]
_flame_previous_enthalpy(f::IonizedFlame,w,previous)=w.conservative===nothing ? nothing :
    _flame_previous_enthalpy!(w.conservative,f,previous)
_flame_correction_enabled(f::IonizedFlame)=true
_flame_needs_correction_check(f::IonizedFlame)=true
_flame_checkpoint_steady(f::IonizedFlame)=true
_flame_reuse_trial(f::IonizedFlame)=true
_flame_refine_threshold(f::IonizedFlame,k)=
    sqrt(eps(Float64))/(k in (1,f.gas.n_species+2) ? 1000 : 1)
function _flame_perturbation(f::IonizedFlame,u,k,j)
    scale=k in (1,f.gas.n_species+2) ? 1000 : 1
    return copysign(1e-5*abs(u[k,j])+1e-10/scale,u[k,j])
end
function _flame_correction_weights(f::IonizedFlame,u,transient)
    n=f.gas.n_species;weights=zeros(n+3)
    for k in 1:n+3
        rtol=1e-4;atol=transient ? 1e-11 : 1e-9
        if 2<=k<=n+1 && f.ion_data.charges[k-1]!=0
            rtol=1e-5;atol=k-1==f.ion_data.electron ? 1e-20 : 1e-16
        end
        k in (1,n+2) && (atol/=1000)
        weights[k]=rtol*sum(abs,@view(u[k,:]))/size(u,2)+atol
    end
    return weights
end
function _flame_bounds(f::IonizedFlame,k,B)
    k==1 && return (.2,2*minimum(@view(f.gas.thermo.Trange[:,end]))/1000)
    k==B && return (-1e20,1e20)
    k==B-1 && return (-1e17,1e17)
    k-1==f.ion_data.electron && return (-1e-14,1.0)
    f.ion_data.charges[k-1]!=0 && return (-1e-10,1.0)
    return (-1e-7,1e5)
end

function solve!(f::IonizedFlame;kwargs...)
    if f.discretization==:conservative && !f.initialized
        requested_discretization=f.discretization
        requested_field=f.field_enabled
        try
            f.discretization=:finite_difference
            f.field_enabled=false
            invoke(solve!,Tuple{AbstractPremixedFlame},f;kwargs...)
            # Preserve an accepted seed if the subsequent conservative solve fails.
            f.initialized=true
        finally
            f.discretization=requested_discretization
            f.field_enabled=requested_field
            f.converged=false
        end
    end
    invoke(solve!,Tuple{AbstractPremixedFlame},f;kwargs...)
    f.initialized=true
    return f
end

export IonizedFlame, set_electric_field!, electric_field
