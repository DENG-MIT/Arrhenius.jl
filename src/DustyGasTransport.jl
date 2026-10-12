# Dusty-gas equations follow Cantera DustyGasTransport.cpp at
# 726522be4e2a13454d8415b7ef799d621f665cf3.
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

"Stationary porous-medium geometry; lengths are in m and permeability in m²."
struct DustyGasParameters
    porosity::Float64
    tortuosity::Float64
    mean_pore_radius::Float64
    mean_particle_diameter::Float64
    permeability::Union{Nothing,Float64}
end

function DustyGasParameters(;porosity,tortuosity=1.,mean_pore_radius,
        mean_particle_diameter=0.,permeability=nothing)
    isfinite(porosity)&&0<porosity<1 || throw(ArgumentError("porosity must lie strictly between zero and one"))
    all(x->isfinite(x)&&x>0,(tortuosity,mean_pore_radius)) ||
        throw(ArgumentError("positive finite tortuosity and pore radius required"))
    isfinite(mean_particle_diameter)&&mean_particle_diameter>=0 ||
        throw(ArgumentError("nonnegative finite particle diameter required"))
    isnothing(permeability) || (isfinite(permeability)&&permeability>=0) ||
        throw(ArgumentError("permeability must be nonnegative or nothing"))
    return DustyGasParameters(Float64(porosity),Float64(tortuosity),Float64(mean_pore_radius),
        Float64(mean_particle_diameter),isnothing(permeability) ? nothing : Float64(permeability))
end

"""
    DustyGasTransport(gas; porosity, mean_pore_radius, tortuosity=1,
                      mean_particle_diameter=0, permeability=nothing,
                      multicomponent_data=nothing)

Native ideal-gas transport through a stationary porous medium. Molecular and
Knudsen diffusion both include porosity/tortuosity. With `permeability=nothing`,
the Darcy term uses the close-packed-sphere estimate from particle diameter.
Provide `MultiTransportData` to evaluate gas-phase thermal conductivity.
The reusable workspace holds `diffusion` [m²/s], effective `binary` [m²/s],
`knudsen` [m²/s], and the gas viscosity [Pa s].
"""
mutable struct DustyGasTransport{G}
    gas::G
    parameters::DustyGasParameters
    gas_transport::TransportWorkspace
    multicomponent_data::Union{Nothing,MultiTransportData}
    multicomponent::Union{Nothing,MultiTransportWorkspace}
    diffusion::Matrix{Float64}
    binary::Matrix{Float64}
    knudsen::Vector{Float64}
    resistance::Matrix{Float64}
    factor::Matrix{Float64}
    X::Vector{Float64}
    concentration::Vector{Float64}
    gradient::Vector{Float64}
    rhs::Vector{Float64}
    gas_viscosity::Float64
    temperature::Float64
    pressure::Float64
    coefficients_valid::Bool
    thermal_temperature::Float64
end

function DustyGasTransport(gas::Solution;multicomponent_data=nothing,kwargs...)
    gas.thermo isa IdealGasThermo || throw(ArgumentError("dusty-gas transport requires ideal-gas thermodynamics"))
    any(e->lowercase(e)=="e",gas.elements) && throw(ArgumentError("charged-species transport is not supported"))
    gas.trans.poly_order==5 || throw(ArgumentError("native transport fits required"))
    p=DustyGasParameters(;kwargs...)
    if !isnothing(multicomponent_data)
        multicomponent_data.species_names==gas.species_names &&
            multicomponent_data.molecular_weights≈gas.MW ||
            throw(ArgumentError("multicomponent data must match gas species and order"))
    end
    multi=isnothing(multicomponent_data) ? nothing : MultiTransportWorkspace(multicomponent_data)
    n=gas.n_species
    return DustyGasTransport(gas,p,TransportWorkspace(gas),multicomponent_data,multi,
        zeros(n,n),zeros(n,n),zeros(n),zeros(n,n),zeros(n,n),zeros(n),zeros(n),
        zeros(n),zeros(n),NaN,NaN,NaN,false,NaN)
end

"Update medium properties and invalidate cached coefficients."
function set_porous_medium!(w::DustyGasTransport;porosity=w.parameters.porosity,
        tortuosity=w.parameters.tortuosity,mean_pore_radius=w.parameters.mean_pore_radius,
        mean_particle_diameter=w.parameters.mean_particle_diameter,permeability=w.parameters.permeability)
    w.parameters=DustyGasParameters(;porosity,tortuosity,mean_pore_radius,mean_particle_diameter,permeability)
    w.coefficients_valid=false
    return w
end

"The prescribed or close-packed-sphere Darcy permeability [m²]."
function dusty_gas_permeability(w::DustyGasTransport)
    p=w.parameters
    return isnothing(p.permeability) ? p.porosity^3*p.mean_particle_diameter^2/
        (72*p.tortuosity*(1-p.porosity)^2) : p.permeability
end

function _dusty_composition(v,n)
    length(v)==n || throw(DimensionMismatch("one fraction per species required"))
    all(x->isfinite(x)&&x>=0,v) && isapprox(sum(v),1;atol=1e-10,rtol=1e-10) ||
        throw(ArgumentError("normalized nonnegative species fractions required"))
end

"""
    dusty_gas_diffusion!(workspace, P, T, X)

Return the full dusty-gas diffusion matrix [m²/s]. `P` is Pa, `T` is K, and
`X` contains normalized mole fractions. The returned matrix belongs to the
workspace and is replaced by the next coefficient or flux evaluation.
"""
function dusty_gas_diffusion!(w::DustyGasTransport,P,T,X)
    n=w.gas.n_species
    _dusty_composition(X,n)
    all(x->isfinite(x)&&x>0,(P,T)) || throw(ArgumentError("positive finite pressure and temperature required"))
    if w.coefficients_valid && P==w.pressure && T==w.temperature &&
            all(w.X[k]==max(X[k],1e-20) for k in 1:n)
        return w.diffusion
    end
    w.coefficients_valid=false
    w.gas_viscosity,_=mixture_transport!(w.gas_transport,w.gas,P,T,X)
    p=w.parameters
    effective=p.porosity/p.tortuosity
    @inbounds for k in 1:n
        w.X[k]=max(X[k],1e-20) # Cantera's Tiny; deliberately not renormalized
        w.knudsen[k]=(2/3)*p.mean_pore_radius*effective*sqrt(8*R*T/(pi*w.gas.MW[k]))
    end
    @inbounds for j in 1:n,k in 1:n
        w.binary[k,j]=effective*w.gas_transport.binary[k,j]
        w.resistance[k,j]=-w.X[k]/w.binary[k,j]
    end
    @inbounds for k in 1:n
        diagonal=1/w.knudsen[k]
        for j in 1:n
            j==k && continue
            diagonal+=w.X[j]/w.binary[k,j]
        end
        w.resistance[k,k]=diagonal
    end
    copyto!(w.factor,w.resistance)
    fill!(w.diffusion,0)
    @inbounds for k in 1:n
        w.diffusion[k,k]=1
    end
    ldiv!(lu!(w.factor),w.diffusion)
    w.pressure=P; w.temperature=T; w.coefficients_valid=true
    return w.diffusion
end

"Gas-phase multicomponent thermal conductivity [W/(m K)]; no porous-solid correction is included."
function dusty_gas_thermal_conductivity(w::DustyGasTransport,P,T,X)
    isnothing(w.multicomponent_data) && throw(ArgumentError("MultiTransportData is required for gas-phase thermal conductivity"))
    _dusty_composition(X,w.gas.n_species)
    all(x->isfinite(x)&&x>0,(P,T)) || throw(ArgumentError("positive finite pressure and temperature required"))
    # Ideal-gas conductivity depends on T and composition, not medium geometry
    # or pressure. Keep its cache separate from the molecular-diffusion cache.
    if T==w.thermal_temperature && all(w.multicomponent.X[k]==max(X[k],1e-20) for k in eachindex(X))
        return w.multicomponent.conductivity
    end
    conductivity=multicomponent_thermal_conductivity!(w.multicomponent,w.multicomponent_data,w.gas,P,T,X)
    w.thermal_temperature=T
    return conductivity
end

"""
    dusty_gas_molar_fluxes!(flux, workspace, T1, T2, rho1, rho2, Y1, Y2, delta)

Write species molar fluxes [kmol/(m² s)] from point 1 toward point 2, separated
by positive `delta` [m]. Temperatures are K, densities kg/m³, and `Y1`, `Y2`
are normalized mass fractions. The calculation includes molecular diffusion,
Knudsen diffusion and pressure-driven Darcy flow. No zero-total-flux correction
is applied: momentum can be transferred to the stationary porous medium.

Coefficients are evaluated at arithmetic mean T/P and mole fractions formed
from mean endpoint concentrations. Darcy's term uses those mean concentrations
directly, including when endpoint temperatures differ.
"""
function dusty_gas_molar_fluxes!(flux,w::DustyGasTransport,T1,T2,rho1,rho2,Y1,Y2,delta)
    n=w.gas.n_species
    length(flux)==n || throw(DimensionMismatch("one molar flux per species required"))
    _dusty_composition(Y1,n); _dusty_composition(Y2,n)
    all(x->isfinite(x)&&x>0,(T1,T2,rho1,rho2,delta)) ||
        throw(ArgumentError("positive finite temperatures, densities and separation required"))
    c1sum=0.; c2sum=0.
    @inbounds for k in 1:n
        c1=rho1*Y1[k]/w.gas.MW[k]
        c2=rho2*Y2[k]/w.gas.MW[k]
        w.concentration[k]=.5*(c1+c2)
        w.gradient[k]=(c2-c1)/delta
        c1sum+=c1; c2sum+=c2
    end
    p1=c1sum*R*T1; p2=c2sum*R*T2
    concentration_sum=sum(w.concentration)
    @inbounds for k in 1:n
        w.rhs[k]=w.concentration[k]/concentration_sum
    end
    dusty_gas_diffusion!(w,.5*(p1+p2),.5*(T1+T2),w.rhs)
    darcy=dusty_gas_permeability(w)*(p2-p1)/(delta*w.gas_viscosity)
    @inbounds for k in 1:n
        w.rhs[k]=w.gradient[k]+w.concentration[k]/w.knudsen[k]*darcy
    end
    mul!(flux,w.diffusion,w.rhs,-1.,0.)
    return flux
end

"Allocate and return species molar fluxes [kmol/(m² s)] between two endpoint states."
function dusty_gas_molar_fluxes(w::DustyGasTransport,args...)
    return dusty_gas_molar_fluxes!(zeros(w.gas.n_species),w,args...)
end

export DustyGasTransport, DustyGasParameters, set_porous_medium!, dusty_gas_permeability
export dusty_gas_diffusion!, dusty_gas_thermal_conductivity, dusty_gas_molar_fluxes!, dusty_gas_molar_fluxes
