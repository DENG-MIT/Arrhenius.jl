# Axisymmetric boundary configurations. Flow1D / Boundary1D and onedim.py at
# Cantera commit 726522be4e2a13454d8415b7ef799d621f665cf3 define the equations
# and native equilibrium initial profiles used here.
# Copyright (c) 2001-2009, California Institute of Technology. All rights reserved.
# Copyright (c) 2009 Sandia Corporation. Under the terms of Contract
# AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain
# rights in this software.
# Copyright (c) 2011-2026, Cantera Developers. All rights reserved.
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# - Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# - Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# - Neither the name of the California Institute of Technology, Sandia
#   Corporation nor the names of other contributors may be used to endorse or
#   promote products derived from this software without specific prior written
#   permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
abstract type AbstractPremixedCounterflow end

"Premixed reactants opposed by a prescribed product stream."
struct CounterflowPremixedFlame{F,C} <: AbstractPremixedCounterflow
    flow::F
    configuration::C
end

"Half of a symmetric pair of premixed counterflow flames."
struct CounterflowTwinPremixedFlame{F,C} <: AbstractPremixedCounterflow
    flow::F
    configuration::C
end

"An axisymmetric jet impinging on an inert, isothermal surface."
struct ImpingingJet{F,C} <: AbstractPremixedCounterflow
    flow::F
    configuration::C
end

@inline Base.getproperty(f::AbstractPremixedCounterflow,name::Symbol) =
    name in (:flow,:configuration) ? getfield(f,name) : getproperty(getfield(f,:flow),name)
@inline Base.setproperty!(f::AbstractPremixedCounterflow,name::Symbol,value) =
    setproperty!(getfield(f,:flow),name,value)
Base.propertynames(f::AbstractPremixedCounterflow,private::Bool=false) =
    (propertynames(getfield(f,:flow),private)...,:flow,:configuration)

function _premixed_counterflow_storage(gas;reactants,mdot_reactants,mdot_products,
        T_reactants,P,width,grid,boundary,products=nothing,T_products=nothing,
        initial_products=:equilibrium)
    all(x->isfinite(x)&&x>0,(P,mdot_reactants,width)) &&
        isfinite(mdot_products)&&mdot_products>=0 ||
        throw(ArgumentError("positive pressure, width and reactant mass flux required; product mass flux must be nonnegative"))
    isfinite(T_reactants)&&200<=T_reactants<=6000 || throw(ArgumentError("inlet temperature must lie in 200–6000 K"))
    gas.thermo isa IdealGasThermo || throw(ArgumentError("ideal gas thermodynamics required"))
    gas.trans.poly_order==5 || throw(ArgumentError("native transport sidecar required"))
    initial_products in (:equilibrium,:inlet) || throw(ArgumentError("initial_products must be :equilibrium or :inlet"))
    fractions=boundary==:products ? [0,.3,.5,.7,1] : boundary==:symmetry ? [0,.2,.4,.5,.6,.8,1] : [0,.2,.4,.6,.8,1]
    z=isnothing(grid) ? Float64(width).*fractions : Float64.(grid)
    length(z)>=5 && all(isfinite,z) && all(>(0),diff(z)) ||
        throw(ArgumentError("at least five increasing grid points required"))
    X=mole_fractions(gas,reactants); Y=X.*gas.MW./dot(X,gas.MW)
    equilibrium=initial_products==:equilibrium ? equilibrate(gas;T=T_reactants,P,X,mode=:HP) : nothing
    Teq=isnothing(equilibrium) ? T_reactants : equilibrium.T
    Yeq=isnothing(equilibrium) ? Y : equilibrium.Y
    Yright=isnothing(products) ? copy(Yeq) : begin
        xr=mole_fractions(gas,products); xr.*gas.MW./dot(xr,gas.MW)
    end
    Tright=isnothing(T_products) ? Teq : Float64(T_products)
    isfinite(Tright)&&200<=Tright<=6000 || throw(ArgumentError("product/surface temperature must lie in 200–6000 K"))
    rhou=P*dot(X,gas.MW)/(R*T_reactants)
    rhob=P/(R*Tright*sum(Yright./gas.MW))
    uu=mdot_reactants/rhou; ub=mdot_products/rhob
    width=z[end]-z[1]
    strain=boundary==:symmetry ? 2uu/width : (uu+ub)/width
    curvature=boundary==:symmetry ? -rhou*strain^2 : -.5*(rhou+rhob)*strain^2
    stagnation=mdot_reactants/(mdot_reactants+mdot_products)
    n=gas.n_species
    state=zeros(n+4,length(z))
    for j in eachindex(z)
        fraction=(z[j]-z[1])/width
        if boundary==:surface
            if initial_products==:equilibrium
                blend=clamp(fraction/.3,0,1)
                temperature= fraction<.7 ? T_reactants+blend*(Teq-T_reactants) : Teq+(fraction-.7)/.3*(Tright-Teq)
                y=Y.+blend.*(Yeq.-Y)
            else
                temperature=T_reactants+fraction*(Tright-T_reactants)
                y=Y
            end
        else
            blend=clamp((fraction-.4)/.2,0,1)
            tail=max(0.,(fraction-.6)/.4)
            temperature=T_reactants+blend*(Teq-T_reactants)+tail*(Tright-Teq)
            y=Y.+blend.*(Yeq.-Y).+tail.*(Yright.-Yeq)
        end
        state[1,j]=temperature/1000
        state[2:n+1,j]=y
        V=boundary==:symmetry ? fraction*strain : boundary==:surface ? 0. :
            fraction<=stagnation ? strain*fraction/stagnation : strain*(1-fraction)/(1-stagnation)
        state[n+2,j]=V/_counterflow_Vscale
        state[n+3,j]=(boundary==:surface ? 0. : curvature)/_counterflow_Lscale
        rho=P/(R*temperature*sum(y./gas.MW))
        state[end,j]=rho*((1-fraction)*uu-fraction*ub)
    end
    state[end,1]=mdot_reactants; state[end,end]=-mdot_products
    return CounterflowDiffusionFlame(gas,z,Float64(P),Float64(T_reactants),Tright,
        Y,Yright,Float64(mdot_reactants),Float64(mdot_products),state,1,argmax(Y),false,
        :mixture_averaged,false,nothing,:mole,Float64[],false,(0.,0.),nothing)
end

"""
    CounterflowPremixedFlame(gas; reactants, mdot_reactants, mdot_products,
                            T_reactants=300, P=one_atm, width=.02)

Opposed premixed reactants and products, with positive inlet mass-flux magnitudes
in kg/(m² s). By default, the product temperature and composition are calculated
by native HP equilibrium. Set both `products` and `T_products` for another
prescribed product stream. No reference flame profiles are used to initialize.
"""
function CounterflowPremixedFlame(gas::Solution;reactants,mdot_reactants,mdot_products,
        T_reactants=300.,P=one_atm,width=.02,grid=nothing,products=nothing,T_products=nothing)
    xor(isnothing(products),isnothing(T_products)) &&
        throw(ArgumentError("products and T_products must be specified together"))
    config=(;reactants,mdot_reactants,mdot_products,T_reactants,P,width,products,T_products)
    flow=_premixed_counterflow_storage(gas;config...,grid,boundary=:products)
    return CounterflowPremixedFlame(flow,config)
end

"""
    CounterflowTwinPremixedFlame(gas; reactants, mdot, T=300, P=one_atm, width=.025)

Solve one half of a symmetric pair of opposed premixed jets. `width` is the
distance from the inlet to the symmetry plane. The plane has zero axial velocity,
zero temperature and spread-rate gradients, and zero diffusive species fluxes.
"""
function CounterflowTwinPremixedFlame(gas::Solution;reactants,mdot,T=300.,P=one_atm,width=.025,grid=nothing)
    config=(;reactants,mdot,T,P,width)
    flow=_premixed_counterflow_storage(gas;reactants,mdot_reactants=mdot,mdot_products=0.,
        T_reactants=T,P,width,grid,boundary=:symmetry)
    return CounterflowTwinPremixedFlame(flow,config)
end

"""
    ImpingingJet(gas; reactants, mdot, T_inlet=300, T_surface=300,
                 P=one_atm, width=.02, initial_products=:equilibrium)

An axisymmetric jet against an inert isothermal wall. The wall has zero axial
velocity and spread rate and zero net species flux. `initial_products=:inlet`
initializes a nonburning flow; `:equilibrium` initializes a detached flame using
native HP equilibrium. Reactive surfaces are not represented by this type.
"""
function ImpingingJet(gas::Solution;reactants,mdot,T_inlet=300.,T_surface=300.,
        P=one_atm,width=.02,grid=nothing,initial_products=:equilibrium)
    config=(;reactants,mdot,T_inlet,T_surface,P,width,initial_products)
    flow=_premixed_counterflow_storage(gas;reactants,mdot_reactants=mdot,mdot_products=0.,
        T_reactants=T_inlet,T_products=T_surface,P,width,grid,boundary=:surface,initial_products)
    return ImpingingJet(flow,config)
end

temperature(f::AbstractPremixedCounterflow)=temperature(f.flow)
mass_fractions(f::AbstractPremixedCounterflow)=mass_fractions(f.flow)
velocity(f::AbstractPremixedCounterflow)=velocity(f.flow)
spread_rate(f::AbstractPremixedCounterflow)=spread_rate(f.flow)
pressure_curvature(f::AbstractPremixedCounterflow)=pressure_curvature(f.flow)
function density(f::AbstractPremixedCounterflow)
    return [f.pressure/(R*1000*f.state[1,j]*sum(f.state[k+1,j]/f.gas.MW[k]
        for k in 1:f.gas.n_species)) for j in eachindex(f.grid)]
end
"""
    heat_release_rate(f)

Volumetric chemical heat release rate [W/m^3] at each grid point.
"""
function heat_release_rate(f::Union{AbstractPremixedCounterflow,CounterflowDiffusionFlame})
    properties=CounterflowWorkspace(f).properties
    _flame_properties!(properties,f,f.state)
    return [-dot(@view(properties.h[:,j]),@view(properties.source[:,j])) for j in eachindex(f.grid)]
end
CounterflowWorkspace(f::AbstractPremixedCounterflow)=CounterflowWorkspace(f.flow)
extinct(f::AbstractPremixedCounterflow)=maximum(temperature(f))-f.fuel_temperature<10
extinct(f::ImpingingJet)=maximum(temperature(f))-max(f.fuel_temperature,f.oxidizer_temperature)<10

"""
Return consumption speed [m/s], characteristic upstream strain rate [1/s], and
its grid index for a twin flame. Definitions follow the published Cantera twin
flame example: integrate heat release divided by mass heat capacity, and use the
largest forward velocity gradient upstream of the preheat velocity minimum.
"""
function twin_flame_diagnostics(f::CounterflowTwinPremixedFlame)
    p=CounterflowWorkspace(f).properties
    _flame_properties!(p,f,f.state)
    N=length(f.grid)
    integrand=zeros(N)
    for j in 1:N
        cp=sum(f.state[k+1,j]*p.cp[k,j]/f.gas.MW[k] for k in 1:f.gas.n_species)
        integrand[j]=-dot(@view(p.h[:,j]),@view(p.source[:,j]))/cp
    end
    T=temperature(f)
    total=sum(.5*(integrand[j]+integrand[j+1])*(f.grid[j+1]-f.grid[j]) for j in 1:N-1)
    consumption=total/((maximum(T)-minimum(T))*maximum(p.rho))
    speed=velocity(f)
    rates=vcat(diff(speed)./diff(f.grid),(speed[end]-speed[end-1])/(f.grid[end]-f.grid[end-1]))
    maximum_location=argmax(abs.(rates))
    maximum_location>2 || throw(ArgumentError("no resolved upstream preheat region"))
    minimum_velocity=argmin(@view speed[1:maximum_location-1])
    minimum_velocity>1 || throw(ArgumentError("no resolved upstream velocity minimum"))
    point=argmax(abs.(rates[1:minimum_velocity-1]))
    return (consumption_speed=consumption,characteristic_strain_rate=abs(rates[point]),
        strain_rate_point=point,strain_rate_profile=rates)
end

function counterflow_residual!(residual,f::AbstractPremixedCounterflow,u=f.state,
        w=CounterflowWorkspace(f);kwargs...)
    counterflow_residual!(residual,f.flow,u,w;kwargs...)
    if (f isa CounterflowTwinPremixedFlame || f isa ImpingingJet) && !f.soret_enabled
        # At an impermeable boundary, the independent zero-diffusion-flux
        # equations and normalization are equivalent to equal adjacent Y.
        # This form stays differentiable for negative Newton trial traces.
        for k in 1:f.gas.n_species
            k==f.dependent_species && continue
            residual[k+1,end]=u[k+1,end]-u[k+1,end-1]
        end
    end
    if f isa CounterflowTwinPremixedFlame
        residual[end-2,end]=u[end-2,end]-u[end-2,end-1]
        if isempty(f.fixed_temperature)
            residual[1,end]=u[1,end]-u[1,end-1]
        end
    end
    return residual
end

function _premixed_reinitialize(f::CounterflowPremixedFlame,grid)
    config=merge(f.configuration,(mdot_reactants=f.fuel_mass_flux,mdot_products=f.oxidizer_mass_flux))
    return CounterflowPremixedFlame(f.gas;config...,grid)
end
function _premixed_reinitialize(f::CounterflowTwinPremixedFlame,grid)
    return CounterflowTwinPremixedFlame(f.gas;merge(f.configuration,(mdot=f.fuel_mass_flux,))...,grid)
end
function _premixed_reinitialize(f::ImpingingJet,grid)
    return ImpingingJet(f.gas;merge(f.configuration,(mdot=f.fuel_mass_flux,))...,grid)
end

"Change positive inlet mass-flux magnitudes [kg/(m² s)] while retaining the current profile."
function set_mass_flux!(f::AbstractPremixedCounterflow;reactants=f.fuel_mass_flux,products=f.oxidizer_mass_flux)
    isfinite(reactants)&&reactants>0 && isfinite(products)&&products>=0 ||
        throw(ArgumentError("positive reactant and nonnegative product mass flux required"))
    !(f isa CounterflowPremixedFlame) && products!=0 &&
        throw(ArgumentError("an impermeable wall or symmetry plane has zero mass flux"))
    f.fuel_mass_flux=reactants
    f.oxidizer_mass_flux=products
    f.converged=false
    return f
end
function _counterflow_reset!(f::AbstractPremixedCounterflow,grid)
    initialized=_premixed_reinitialize(f,grid)
    f.grid=initialized.grid; f.state=initialized.state; f.anchor=1
    empty!(f.fixed_temperature)
    return f
end

function _premixed_counterflow_steady!(f;kwargs...)
    try
        return _counterflow_steady!(f;kwargs...)
    catch error
        # A singular coarse-grid factorization is a failed initialization.
        error isa LinearAlgebra.LAPACKException || rethrow()
        return false
    end
end

# Cantera Refiner::analyze insertion/pruning criteria, evaluated in physical
# units (including velocity instead of the axial mass-flux unknown).
function _refine_axisymmetric!(f;ratio,slope,curve,prune,grid_min,max_points)
    u,z=f.state,f.grid
    B,N=size(u)
    insert=falses(N-1)
    keep=zeros(Int8,N) # 0 unset, -1 remove, 1 keep
    keep[1]=keep[end]=1
    spacing=diff(z)
    threshold=sqrt(eps(Float64))
    for k in 1:B
        values=k==B ? velocity(f) : k==1 ? temperature(f) :
            k==B-2 ? spread_rate(f) : k==B-1 ? pressure_curvature(f) : @view(u[k,:])
        lo,hi=extrema(values)
        if hi-lo>.01*max(abs(lo),abs(hi))
            change=slope*(hi-lo)+threshold
            for j in 1:N-1
                quotient=abs(values[j+1]-values[j])/change
                quotient>1 && spacing[j]>=2grid_min && (insert[j]=true)
                if quotient>=prune
                    keep[j]=keep[j+1]=1
                elseif keep[j]==0
                    keep[j]=-1
                end
            end
        end
        gradients=diff(values)./spacing
        lo,hi=extrema(gradients)
        if hi-lo>.01*max(abs(lo),abs(hi))
            change=curve*(hi-lo)
            for j in 1:N-2
                quotient=abs(gradients[j+1]-gradients[j])/(change+threshold/spacing[j])
                if quotient>1 && min(spacing[j],spacing[j+1])>=2grid_min
                    insert[j]=insert[j+1]=true
                end
                if quotient>=prune
                    keep[j+1]=1
                elseif keep[j+1]==0
                    keep[j+1]=-1
                end
            end
        end
    end
    for j in 2:N-1
        if spacing[j]>ratio*spacing[j-1]
            insert[j]=true
            keep[max(1,j-1):min(N,j+2)] .= 1
        end
        if spacing[j-1]>ratio*spacing[j]
            insert[j-1]=true
            keep[max(1,j-2):min(N,j+1)] .= 1
        end
        j>2 && z[j+1]-z[j-1]>ratio*spacing[j-2] && (keep[j]=1)
        j<N-1 && z[j+1]-z[j-1]>ratio*spacing[j+1] && (keep[j]=1)
    end
    for j in 3:N-1
        keep[j]==-1 && keep[j-1]==-1 && (keep[j]=1)
    end
    !any(insert) && all(x->x>=0,keep) && return false
    newgrid=Float64[]; columns=Vector{Float64}[]
    for j in 1:N
        keep[j]==-1 && continue
        push!(newgrid,z[j]); push!(columns,u[:,j])
        if j<N && insert[j]
            push!(newgrid,(z[j]+z[j+1])/2)
            push!(columns,(u[:,j]+u[:,j+1])/2)
        end
    end
    length(newgrid)<=max_points || error("axisymmetric refinement exceeds max_points=$max_points")
    f.grid=newgrid; f.state=reduce(hcat,columns); f.anchor=1
    return true
end

"Solve the axisymmetric flow with native equilibrium initialization and adaptive refinement."
function solve!(f::AbstractPremixedCounterflow;refine_grid=true,ratio=3.,slope=.1,curve=.2,
        prune=0.,grid_min=1e-10,max_points=1200,max_time_steps=800,loglevel=0,
        initial_time_step=1e-6,auto=true)
    isfinite(ratio)&&ratio>1 && 0<slope<=1 && 0<curve<=1 || throw(ArgumentError("invalid refinement criteria"))
    isfinite(initial_time_step)&&initial_time_step>0 || throw(ArgumentError("positive initial time step required"))
    isfinite(prune)&&prune<=min(slope,curve) && isfinite(grid_min)&&grid_min>0 ||
        throw(ArgumentError("finite prune <= min(slope,curve) and positive grid_min required"))
    f.converged=false
    timestep=Ref(Float64(initial_time_step))
    initial_N=length(f.grid)
    attempts=auto&&refine_grid ? unique(vcat(initial_N,[N for N in (12,24,48) if N>initial_N])) : [initial_N]
    require_burning=!(f isa ImpingingJet && f.configuration.initial_products==:inlet)
    acceptable()=!require_burning || !extinct(f)
    success=false
    for (attempt,N) in enumerate(attempts)
        N<=max_points || error("counterflow initialization exceeds max_points=$max_points")
        if attempt>1
            _counterflow_reset!(f,collect(range(f.grid[1],f.grid[end];length=N)))
            timestep[]=initial_time_step
        end
        loglevel>0 && println("Premixed counterflow solve on ",length(f.grid)," points")
        success=_premixed_counterflow_steady!(f;max_time_steps,timestep,loglevel)
        if auto && (!success || !acceptable())
            loglevel>0 && println("Recovering with native initial temperature profile")
            _counterflow_reset!(f,copy(f.grid))
            f.fixed_temperature=temperature(f)
            success=_premixed_counterflow_steady!(f;max_time_steps,timestep,loglevel)
            empty!(f.fixed_temperature)
            success && (success=_premixed_counterflow_steady!(f;max_time_steps,timestep,loglevel))
        end
        success && acceptable() && break
    end
    success && acceptable() || error("native premixed counterflow solver failed")
    for pass in 1:40
        if !refine_grid || !_refine_axisymmetric!(f;ratio,slope,curve,prune,grid_min,max_points)
            f.converged=true
            return f
        end
        loglevel>0 && println("Premixed counterflow solve on ",length(f.grid)," points")
        _premixed_counterflow_steady!(f;max_time_steps,timestep,loglevel) || error("native premixed counterflow failed on $(length(f.grid)) points")
        acceptable() || error("premixed counterflow extinguished during refinement")
    end
    error("premixed counterflow refinement did not finish")
end

export CounterflowPremixedFlame, CounterflowTwinPremixedFlame, ImpingingJet, set_mass_flux!
export twin_flame_diagnostics
