# Dixon–Lewis multicomponent transport and mixture Soret, following Cantera
# MultiTransport.cpp and MixTransport.cpp
# at 726522be4e2a13454d8415b7ef799d621f665cf3. The translated equations retain
# Cantera's license:
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

"Temperature-independent neutral ideal-gas collision fits and species constants."
struct MultiTransportData
    species_names::Vector{String}
    molecular_weights::Vector{Float64}
    epsilon_over_k::Vector{Float64}
    log_epsilon_ij_over_k::Matrix{Float64}
    rotational_heat_capacity::Vector{Float64}
    rotational_relaxation::Vector{Float64}
    viscosity_poly::Matrix{Float64}
    binary_poly::Array{Float64,3}
    astar_poly::Array{Float64,3}
    bstar_poly::Array{Float64,3}
    cstar_poly::Array{Float64,3}
    gas_constant::Float64
end

function MultiTransportData(path::AbstractString; mechanism=nothing)
    a = npzread(path)
    text(key) = String(vec(UInt8.(a[key])))
    text("format_utf8") == "arrhenius-multicomponent-v1" ||
        throw(ArgumentError("unsupported multicomponent sidecar format"))
    if mechanism !== nothing && haskey(a, "source_sha256_utf8")
        _matches_mechanism_hash(text("source_sha256_utf8"), mechanism) ||
            throw(ArgumentError("multicomponent sidecar does not match mechanism"))
    end
    names = String.(split(text("species_names_utf8"), '\n'))
    n = length(names)
    data = MultiTransportData(names, vec(a["molecular_weights"]),
        vec(a["epsilon_over_k"]), a["log_epsilon_ij_over_k"],
        vec(a["rotational_heat_capacity"]), vec(a["rotational_relaxation"]),
        a["viscosity_poly"], a["binary_poly"], a["astar_poly"],
        a["bstar_poly"], a["cstar_poly"], only(a["gas_constant"]))
    n > 1 || throw(ArgumentError("multicomponent transport requires at least two species"))
    size(data.viscosity_poly) == (5,n) && size(data.binary_poly) == (5,n,n) &&
        all(size(p) == (9,n,n) for p in (data.astar_poly,data.bstar_poly,data.cstar_poly)) ||
        throw(DimensionMismatch("invalid multicomponent polynomial dimensions"))
    all(length(v) == n for v in (data.molecular_weights,data.epsilon_over_k,
        data.rotational_heat_capacity,data.rotational_relaxation)) &&
        size(data.log_epsilon_ij_over_k) == (n,n) ||
        throw(DimensionMismatch("invalid multicomponent species dimensions"))
    all(all(isfinite,a) for a in (data.molecular_weights,data.epsilon_over_k,
        data.rotational_heat_capacity,data.rotational_relaxation,
        data.log_epsilon_ij_over_k,data.viscosity_poly,data.binary_poly,
        data.astar_poly,data.bstar_poly,data.cstar_poly)) &&
        all(>(0),data.molecular_weights) && all(>(0),data.epsilon_over_k) &&
        isfinite(data.gas_constant) && data.gas_constant > 0 ||
        throw(ArgumentError("invalid multicomponent transport constants"))
    return data
end

function MultiTransportData(path::AbstractString, gas::Solution; kwargs...)
    data = MultiTransportData(path; kwargs...)
    data.species_names == gas.species_names || throw(ArgumentError("transport species order mismatch"))
    data.molecular_weights ≈ gas.MW || throw(ArgumentError("transport molecular weights mismatch"))
    return data
end

"Reusable storage. `diffusion` is m²/s; `thermal_diffusion` is kg/(m s); `binary` stores pressure times diffusivity."
mutable struct MultiTransportWorkspace
    X::Vector{Float64}
    viscosity::Vector{Float64}
    binary::Matrix{Float64}
    astar::Matrix{Float64}
    bstar::Matrix{Float64}
    cstar::Matrix{Float64}
    rotational_relaxation::Vector{Float64}
    internal_heat_capacity::Vector{Float64}
    cp_R::Vector{Float64}
    L::Matrix{Float64}
    rhs::Vector{Float64}
    L00::Matrix{Float64}
    diffusion::Matrix{Float64}
    thermal_diffusion::Vector{Float64}
    conductivity::Float64
end

function MultiTransportWorkspace(data::MultiTransportData)
    n = length(data.molecular_weights)
    return MultiTransportWorkspace(zeros(n),zeros(n),zeros(n,n),zeros(n,n),zeros(n,n),
        zeros(n,n),zeros(n),zeros(n),zeros(n),zeros(3n,3n),zeros(3n),
        zeros(n,n),zeros(n,n),zeros(n),NaN)
end

@inline function _multi_poly(a, i, j, z)
    p = a[end,i,j]
    @inbounds for k in size(a,1)-1:-1:1
        p = p*z + a[k,i,j]
    end
    return p
end
@inline _multi_frot(tr) = 1 + (0.5*pi*sqrt(pi))*sqrt(tr) +
    (0.25*pi*pi+2)*tr + (pi*sqrt(pi))*sqrt(tr)*tr

function _multi_temperature!(w, data, T, cp_R)
    n = length(w.X)
    logT = log(T)
    sqrtT = sqrt(T)
    @inbounds for i in 1:n
        p = data.viscosity_poly[5,i]
        for k in 4:-1:1
            p = p*logT + data.viscosity_poly[k,i]
        end
        w.viscosity[i] = (sqrt(sqrtT)*p)^2
        w.rotational_relaxation[i] = max(1.,data.rotational_relaxation[i]) *
            _multi_frot(data.epsilon_over_k[i]/298) / _multi_frot(data.epsilon_over_k[i]/T)
        w.internal_heat_capacity[i] = cp_R[i] - 2.5
    end
    @inbounds for j in 1:n, i in 1:j
        z = logT - data.log_epsilon_ij_over_k[i,j]
        w.astar[i,j] = w.astar[j,i] = _multi_poly(data.astar_poly,i,j,z)
        w.bstar[i,j] = w.bstar[j,i] = _multi_poly(data.bstar_poly,i,j,z)
        w.cstar[i,j] = w.cstar[j,i] = _multi_poly(data.cstar_poly,i,j,z)
        w.binary[i,j] = w.binary[j,i] = T*sqrtT*_multi_poly(data.binary_poly,i,j,logT)
    end
    @inbounds for i in 1:n
        w.binary[i,i] = 1.2*data.gas_constant*T*w.viscosity[i]*w.astar[i,i]/data.molecular_weights[i]
    end
end

function _multi_L00!(L, w, data, T)
    n = length(w.X)
    x, mw, D = w.X, data.molecular_weights, w.binary
    prefactor = 16*T/25
    @inbounds for i in 1:n
        s = -x[i]/D[i,i]
        for k in 1:n
            s += x[k]/D[i,k]
        end
        s /= mw[i]
        for j in 1:n
            L[i,j] = prefactor*x[j]*(mw[j]*s + x[i]/D[i,j])
        end
        L[i,i] = 0
    end
end

function _multi_thermal!(w, data, T)
    n = length(w.X)
    x, mw, D = w.X, data.molecular_weights, w.binary
    A, B, C, L = w.astar,w.bstar,w.cstar,w.L
    cr, ci, zr = data.rotational_heat_capacity,w.internal_heat_capacity,w.rotational_relaxation
    fill!(L,0)
    _multi_L00!(L,w,data,T)
    @inbounds for j in 1:n
        s = 0.
        for i in 1:n
            L[i,j+n] = -1.6*T*x[i]*x[j]*mw[i]*(1.2*C[j,i]-1)/((mw[j]+mw[i])*D[j,i])
            s -= L[i,j+n]
        end
        L[j,j+n] += s
    end
    @inbounds for j in 1:n, i in 1:n
        L[i+n,j] = L[j,i+n]
    end
    @inbounds for j in 1:n
        c1 = 16*T/25*x[j]
        c2 = 13.75*mw[j]^2
        c3 = cr[j]/zr[j]
        c4 = 7.5*mw[j]^2
        s = 0.
        for i in 1:n
            t1 = D[i,j]*(mw[i]+mw[j])^2
            t2 = 4*mw[j]*A[i,j]*(1+5/(3*pi)*(c3+cr[i]/zr[i]))
            L[i+n,j+n] = c1*x[i]*mw[i]/(mw[j]*t1)*(c2-3*mw[j]^2*B[i,j]-t2*mw[j])
            s += x[i]/t1*(c4+mw[i]^2*(6.25-3*B[i,j])+t2*mw[i])
        end
        L[j+n,j+n] -= s*c1
        if ci[j] > .001
            c = 32*T/(5*pi)*mw[j]*x[j]*cr[j]/(ci[j]*zr[j])
            s = 0.
            for i in 1:n
                L[i+n,j+2n] = c*A[j,i]*x[i]/((mw[j]+mw[i])*D[j,i])
                s += L[i+n,j+2n]
            end
            L[j+n,j+2n] += s
        end
    end
    @inbounds for j in 1:n, i in 1:n
        L[i+2n,j+n] = L[j+n,i+2n]
    end
    @inbounds for i in 1:n
        if ci[i] > .001
            c1 = 4*T*x[i]/ci[i]
            c2 = 12*mw[i]*cr[i]/(5*pi*ci[i]*zr[i])
            s = 0.
            for k in 1:n
                s += x[k]/D[i,k]
                if k != i
                    s += x[k]*A[i,k]*c2/(mw[k]*D[i,k])
                end
            end
            L[i+2n,i+2n] = -8/pi*mw[i]*x[i]^2*cr[i]/
                (ci[i]^2*data.gas_constant*w.viscosity[i]*zr[i])-c1*s
        else
            L[i+2n,i+2n] = 1
        end
        w.rhs[i] = 0
        w.rhs[i+n] = x[i]
        w.rhs[i+2n] = ci[i] > .001 ? x[i] : 0
    end
    ldiv!(lu!(L),w.rhs)
    conductivity = 0.
    @inbounds for i in 1:n
        conductivity += x[i]*w.rhs[i+n]
        if ci[i] > .001
            conductivity += x[i]*w.rhs[i+2n]
        end
        w.thermal_diffusion[i] = 1.6/data.gas_constant*mw[i]*x[i]*w.rhs[i]
    end
    w.conductivity = -4*conductivity
end

function _multicomponent_state!(w::MultiTransportWorkspace,data::MultiTransportData,P,T,X,cp_R)
    n = length(data.molecular_weights)
    length(X) == length(cp_R) == length(w.X) == n || throw(DimensionMismatch("transport species size mismatch"))
    isfinite(P) && P > 0 && isfinite(T) && T > 0 || throw(ArgumentError("positive finite P and T required"))
    all(isfinite,X) && all(x -> x >= 0,X) && isapprox(sum(X),1;atol=1e-10,rtol=1e-10) ||
        throw(ArgumentError("normalized nonnegative mole fractions required"))
    all(isfinite,cp_R) || throw(ArgumentError("finite species heat capacities required"))
    @inbounds for i in 1:n
        w.X[i] = max(X[i],1e-20)
    end
    _multi_temperature!(w,data,T,cp_R)
    return n
end

"""
    multicomponent_transport!(workspace, data, P, T, X, cp_R)

Compute the Dixon–Lewis multicomponent matrix [m²/s], thermal diffusion
coefficients [kg/(m s)], and conductivity [W/(m K)]. Return conductivity and
store the other results in `workspace.diffusion` and `workspace.thermal_diffusion`.
`cp_R` contains species reference ideal-gas heat capacities divided by R.
Mole fractions must be normalized and nonnegative; zero entries use the same
1e-20 transport floor as Cantera. This floor is not renormalized.
"""
function multicomponent_transport!(w::MultiTransportWorkspace,data::MultiTransportData,P,T,X,cp_R)
    n = _multicomponent_state!(w,data,P,T,X,cp_R)
    _multi_L00!(w.L00,w,data,T)
    fill!(w.diffusion,0)
    @inbounds for i in 1:n
        w.diffusion[i,i] = 1
    end
    ldiv!(lu!(w.L00),w.diffusion)
    prefactor = 16*T*dot(X,data.molecular_weights)/(25*P)
    @inbounds for i in 1:n
        diagonal = w.diffusion[i,i]
        for j in 1:n
            w.diffusion[i,j] = prefactor/data.molecular_weights[j]*w.X[i]*(w.diffusion[i,j]-diagonal)
        end
    end
    _multi_thermal!(w,data,T)
    return w.conductivity
end

function multicomponent_transport!(w::MultiTransportWorkspace,data::MultiTransportData,
                                   gas::Solution,P,T,X)
    gas.species_names == data.species_names || throw(ArgumentError("transport species order mismatch"))
    cal_cp_R!(w.cp_R,gas,T,P,X)
    return multicomponent_transport!(w,data,P,T,X,w.cp_R)
end

"""
    multicomponent_thermal_conductivity!(workspace, data, gas, P, T, X)
    multicomponent_thermal_conductivity!(workspace, data, P, T, X, cp_R)

Return Dixon–Lewis thermal conductivity [W/(m K)] without computing the
multicomponent diffusion matrix. This updates thermal diffusion coefficients
and the thermal workspace; `workspace.diffusion` retains its previous value.
Input requirements are the same as for `multicomponent_transport!`.
"""
function multicomponent_thermal_conductivity!(w::MultiTransportWorkspace,data::MultiTransportData,P,T,X,cp_R)
    _multicomponent_state!(w,data,P,T,X,cp_R)
    _multi_thermal!(w,data,T)
    return w.conductivity
end

function multicomponent_thermal_conductivity!(w::MultiTransportWorkspace,data::MultiTransportData,
                                            gas::Solution,P,T,X)
    gas.species_names == data.species_names || throw(ArgumentError("transport species order mismatch"))
    cal_cp_R!(w.cp_R,gas,T,P,X)
    return multicomponent_thermal_conductivity!(w,data,P,T,X,w.cp_R)
end

"""
    multicomponent_fluxes!(flux, workspace, data, P, T, X, grad_X, grad_T)

Return mass fluxes [kg/(m² s)] relative to the mass-averaged velocity, using
already updated transport coefficients. `grad_X` is in 1/m and `grad_T` in K/m.
The Soret contribution is `-thermal_diffusion * grad_T / T`.
"""
function multicomponent_fluxes!(flux,w::MultiTransportWorkspace,data::MultiTransportData,
                                P,T,X,grad_X,grad_T)
    n = length(w.X)
    length(flux) == length(X) == length(grad_X) == n || throw(DimensionMismatch("flux species size mismatch"))
    isfinite(P) && P > 0 && isfinite(T) && T > 0 && isfinite(grad_T) &&
        all(isfinite,grad_X) || throw(ArgumentError("invalid flux state or gradients"))
    mw = data.molecular_weights
    mean_mw = dot(X,mw)
    rho = P*mean_mw/(data.gas_constant*T)
    @inbounds for i in 1:n
        s = 0.
        for j in 1:n
            s += mw[j]*w.diffusion[i,j]*grad_X[j]
        end
        flux[i] = rho*mw[i]/mean_mw^2*s - w.thermal_diffusion[i]*grad_T/T
    end
    return flux
end

export MultiTransportData, MultiTransportWorkspace, multicomponent_transport!, multicomponent_fluxes!
export multicomponent_thermal_conductivity!

"Reusable storage for the separate mixture-averaged Soret model."
struct MixtureThermalDiffusionWorkspace
    X::Vector{Float64}
    Y::Vector{Float64}
    viscosity::Vector{Float64}
    binary::Matrix{Float64}
    phi::Matrix{Float64}
    a::Vector{Float64}
    diffusion::Vector{Float64}
    thermal_diffusion::Vector{Float64}
end

function MixtureThermalDiffusionWorkspace(data::MultiTransportData)
    n=length(data.molecular_weights)
    return MixtureThermalDiffusionWorkspace(zeros(n),zeros(n),zeros(n),zeros(n,n),
        zeros(n,n),zeros(n),zeros(n),zeros(n))
end

"""
    mixture_thermal_diffusion!(workspace, data, P, T, X)

Compute the mixture-averaged Soret coefficients [kg/(m s)], following
Cantera 4 `MixTransport::getThermalDiffCoeffs`. Return the updated vector
`workspace.thermal_diffusion`. These are distinct from multicomponent Soret
coefficients; their mass-flux contribution is `-DT * grad_T / T` for either
mole- or mass-gradient mixture diffusion. Requires the collision sidecar,
but no species heat capacities or coupled multicomponent matrix solve.
"""
function mixture_thermal_diffusion!(w::MixtureThermalDiffusionWorkspace,
                                    data::MultiTransportData,P,T,X)
    n=length(data.molecular_weights)
    length(X) == length(w.X) == n || throw(DimensionMismatch("transport species size mismatch"))
    isfinite(P) && P>0 && isfinite(T) && T>0 || throw(ArgumentError("positive finite P and T required"))
    all(isfinite,X) && all(x -> x>=0,X) && isapprox(sum(X),1;atol=1e-10,rtol=1e-10) ||
        throw(ArgumentError("normalized nonnegative mole fractions required"))
    mw=data.molecular_weights
    mmw=dot(X,mw)
    logT=log(T); sqrtT=sqrt(T)
    @inbounds for k in 1:n
        w.X[k]=max(1e-20,X[k])
        w.Y[k]=X[k]*mw[k]/mmw
        p=data.viscosity_poly[5,k]
        for i in 4:-1:1
            p=p*logT+data.viscosity_poly[i,k]
        end
        w.viscosity[k]=(sqrt(sqrtT)*p)^2
        w.thermal_diffusion[k]=0
    end
    @inbounds for j in 1:n, k in j:n
        vr=w.viscosity[k]/w.viscosity[j]
        wr=mw[j]/mw[k]
        factor=1+sqrt(vr)*sqrt(sqrt(wr))
        w.phi[k,j]=factor^2/(sqrt(8)*sqrt(1+mw[k]/mw[j]))
        w.phi[j,k]=w.phi[k,j]/(vr*wr)
        w.binary[k,j]=w.binary[j,k]=T*sqrtT*_multi_poly(data.binary_poly,k,j,logT)
    end
    @inbounds for k in 1:n
        if w.Y[k]<1e-20
            w.a[k]=0
            continue
        end
        s=0.
        for j in 1:n
            j==k && continue
            s+=w.X[j]*w.phi[k,j]
        end
        w.a[k]=(15/4)*w.viscosity[k]/mw[k]/(1+1.065*s/w.X[k])
    end
    @inbounds for k in 1:n-1, j in k+1:n
        z=logT-data.log_epsilon_ij_over_k[k,j]
        Cstar=_multi_poly(data.cstar_poly,k,j,z)
        dt_T=(1.2*Cstar-1)/(w.binary[k,j]/P)/(mw[k]+mw[j])
        delta=dt_T*(w.Y[k]*w.a[j]-w.Y[j]*w.a[k])
        w.thermal_diffusion[k]+=delta
        w.thermal_diffusion[j]-=delta
    end
    norm=0.
    @inbounds for k in 1:n
        s=0.
        for j in 1:n
            j==k && continue
            s+=w.X[j]/w.binary[j,k]
        end
        w.diffusion[k]=s<=0 ? w.binary[k,k]/P : (mmw-w.X[k]*mw[k])/(P*mmw*s)
        w.thermal_diffusion[k]*=w.diffusion[k]*mw[k]*mmw
        norm+=w.thermal_diffusion[k]
    end
    @inbounds for k in 1:n
        w.thermal_diffusion[k]-=w.Y[k]*norm
    end
    return w.thermal_diffusion
end

export MixtureThermalDiffusionWorkspace, mixture_thermal_diffusion!
