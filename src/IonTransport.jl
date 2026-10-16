# Ionized-gas transport and charged diffusion fluxes, following Cantera
# IonGasTransport.cpp, GasTransport.cpp, MixTransport.cpp and IonFlow.cpp
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


const _ION_KB = 1.380649e-23
const _ION_QE = 1.602176634e-19
const _ION_ELECTRON_MW = 0.0005485799088728283 # Cantera's standard E atomic weight [kg/kmol]

"Model-specific collision fits and charge information for ionized-gas transport."
struct IonTransportData
    species_names::Vector{String}
    molecular_weights::Vector{Float64}
    charges::Vector{Float64}
    electron::Int
    neutrals::Vector{Int}
    ions::Vector{Int}
    viscosity_poly::Matrix{Float64}
    conductivity_poly::Matrix{Float64}
    binary_poly::Matrix{Float64}
    weight_quarter_ratio::Matrix{Float64}
    inverse_phi_denominator::Matrix{Float64}
end

"""
    IonTransportData(gas)

Prepare ionized-gas transport from a mechanism preprocessed with
mechanism/export_sidecar.py. The collision fits must come from the
ionized-gas model. Electron mobility follows that model's fixed 0.4 m²/(V s).
"""
function IonTransportData(gas::Solution)
    gas.trans.model == :ionized_gas || throw(ArgumentError("ionized-gas transport fits are required"))
    n = gas.n_species
    n > 0 || throw(ArgumentError("at least one species is required"))
    t = gas.trans
    t.poly_order == 5 && size(t.species_viscosities_poly) == size(t.thermal_conductivity_poly) == (5,n) &&
        size(t.binary_diff_coeffs_poly) == (5,n*n) || throw(DimensionMismatch("invalid ionized transport fits; regenerate the sidecar"))
    all(all(isfinite,p) for p in (t.species_viscosities_poly,t.thermal_conductivity_poly,t.binary_diff_coeffs_poly)) &&
        all(isfinite,gas.MW) && all(>(0),gas.MW) || throw(ArgumentError("invalid ionized transport data"))
    e = findfirst(==("E"),gas.elements)
    charges = e === nothing ? zeros(n) : -Float64.(gas.ele_matrix[e,:])
    all(isfinite,charges) || throw(ArgumentError("species charges must be finite"))
    electrons = findall(k -> gas.MW[k] == eltype(gas.MW)(_ION_ELECTRON_MW) && charges[k] == -1,1:n)
    length(electrons) <= 1 || throw(ArgumentError("multiple electron species are not supported"))
    electron = isempty(electrons) ? 0 : only(electrons)
    mw = Float64.(gas.MW)
    IonTransportData(copy(gas.species_names),mw,charges,electron,findall(iszero,charges),
        findall(k -> charges[k] != 0 && k != electron,1:n),
        Matrix{Float64}(t.species_viscosities_poly),Matrix{Float64}(t.thermal_conductivity_poly),
        Matrix{Float64}(t.binary_diff_coeffs_poly),
        [sqrt(sqrt(mw[j]/mw[k])) for k in 1:n,j in 1:n],
        [inv(sqrt(8*(1+mw[k]/mw[j]))) for k in 1:n,j in 1:n])
end

"Reusable ionized transport storage; binary stores D*P [m² Pa/s], diffusion D [m²/s], mobility [m²/(V s)]."
mutable struct IonTransportWorkspace
    X::Vector{Float64}
    viscosity::Vector{Float64}
    sqrt_viscosity::Vector{Float64}
    conductivity::Vector{Float64}
    diffusion::Vector{Float64}
    mobility::Vector{Float64}
    binary::Matrix{Float64}
    pressure::Float64
    temperature::Float64
    mean_MW::Float64
end
function IonTransportWorkspace(data::IonTransportData)
    n = length(data.species_names)
    IonTransportWorkspace(zeros(n),zeros(n),zeros(n),zeros(n),zeros(n),zeros(n),zeros(n,n),NaN,NaN,NaN)
end

"""
    ionized_transport!(workspace, data, P, T, X)

Return mixture viscosity [Pa s], thermal conductivity [W/(m K)], and electrical
conductivity [S/m], and update diffusion coefficients and mobilities. Input mole
fractions sum to one; signed solver traces are floored only in the transport
sums, following Cantera. All properties use gas temperature. The optional
mean_molecular_weight preserves a solver state's molecular weight when its
underlying mass fractions are not normalized.
"""
function ionized_transport!(w::IonTransportWorkspace, data::IonTransportData, P, T, X; mean_molecular_weight=nothing)
    n = length(data.species_names)

    length(data.molecular_weights) == n || throw(DimensionMismatch("molecular_weights length mismatch"))
    length(data.charges) == n || throw(DimensionMismatch("charges length mismatch"))
    size(data.viscosity_poly) == (5, n) || throw(DimensionMismatch("viscosity_poly must be (5, n)"))
    size(data.conductivity_poly) == (5, n) || throw(DimensionMismatch("conductivity_poly must be (5, n)"))
    size(data.binary_poly) == (5, n * n) || throw(DimensionMismatch("binary_poly must be (5, n*n)"))
    size(data.weight_quarter_ratio) == (n, n) || throw(DimensionMismatch("weight_quarter_ratio must be (n, n)"))
    size(data.inverse_phi_denominator) == (n, n) || throw(DimensionMismatch("inverse_phi_denominator must be (n, n)"))
    0 <= data.electron <= n || throw(DimensionMismatch("electron index out of range"))
    for k in data.neutrals
        1 <= k <= n || throw(DimensionMismatch("neutral index out of range"))
    end
    for k in data.ions
        1 <= k <= n || throw(DimensionMismatch("ion index out of range"))
    end
    length(X) == n || throw(DimensionMismatch("mole fraction vector length mismatch"))
    length(w.X) == n || throw(DimensionMismatch("workspace X length mismatch"))
    length(w.viscosity) == n || throw(DimensionMismatch("workspace viscosity length mismatch"))
    length(w.sqrt_viscosity) == n || throw(DimensionMismatch("workspace sqrt_viscosity length mismatch"))
    length(w.conductivity) == n || throw(DimensionMismatch("workspace conductivity length mismatch"))
    length(w.diffusion) == n || throw(DimensionMismatch("workspace diffusion length mismatch"))
    length(w.mobility) == n || throw(DimensionMismatch("workspace mobility length mismatch"))
    size(w.binary) == (n, n) || throw(DimensionMismatch("workspace binary must be (n, n)"))

    (isfinite(P) && P > 0) || throw(ArgumentError("pressure must be positive and finite"))
    (isfinite(T) && T > 0) || throw(ArgumentError("temperature must be positive and finite"))

    xsum = 0.0
    @inbounds for k in 1:n
        xk = X[k]
        isfinite(xk) || throw(ArgumentError("mole fractions must be finite"))
        xsum += xk
    end
    isapprox(xsum, 1.0; atol = 1e-10, rtol = 1e-10) ||
        throw(ArgumentError("mole fractions must be normalized"))

    MW = data.molecular_weights
    Xw = w.X
    mean_MW = 0.0
    @inbounds for k in 1:n
        xk = X[k]
        mean_MW += xk * MW[k]
        Xw[k] = max(xk, 1.0e-20)
    end

    if mean_molecular_weight !== nothing
        mean_MW = Float64(mean_molecular_weight)
    end
    isfinite(mean_MW) && mean_MW > 0 || throw(ArgumentError("mean molecular weight must be finite and positive"))
    logT = log(T)
    sqrtT = sqrt(T)
    T32 = T * sqrtT
    fourthT = sqrt(sqrtT)

    sqv = w.sqrt_viscosity
    visc = w.viscosity
    cond = w.conductivity
    @inbounds for k in 1:n
        sv = fourthT * _transport_polynomial(data.viscosity_poly, k, logT)
        sqv[k] = sv
        visc[k] = sv * sv
        cond[k] = sqrtT * _transport_polynomial(data.conductivity_poly, k, logT)
    end

    bin = w.binary
    @inbounds for j in 1:n
        for i in 1:n
            bin[i, j] = T32 * _transport_polynomial(data.binary_poly, i + (j - 1) * n, logT)
        end
    end

    wqr = data.weight_quarter_ratio
    ipd = data.inverse_phi_denominator
    vismix = 0.0
    @inbounds for k in data.neutrals
        svk = sqv[k]
        denom = 0.0
        for j in 1:n
            r = 1.0 + svk / sqv[j] * wqr[k, j]
            denom += r * r * ipd[k, j] * Xw[j]
        end
        vismix += Xw[k] * visc[k] / denom
    end

    sum1 = 0.0
    sum2 = 0.0
    @inbounds for k in data.neutrals
        sum1 += Xw[k] * cond[k]
        sum2 += Xw[k] / cond[k]
    end
    lambda = 0.5 * (sum1 + 1.0 / sum2)

    kBT = _ION_KB * T
    e = data.electron
    diff = w.diffusion
    if n == 1
        @inbounds diff[1] = bin[1, 1] / P
    else
        @inbounds for k in 1:n
            if k == e
                diff[k] = 0.4 * kBT / _ION_QE
            else
                s2 = 0.0
                for j in data.neutrals
                    if j != k
                        s2 += Xw[j] / bin[j, k]
                    end
                end
                if s2 <= 0.0
                    diff[k] = bin[k, k] / P
                else
                    diff[k] = (mean_MW - Xw[k] * MW[k]) / (P * mean_MW * s2)
                end
            end
        end
    end

    mob = w.mobility
    @inbounds for k in 1:n
        mob[k] = k == e ? 0.4 : 0.0
    end
    qe_over_kBT = _ION_QE / kBT
    @inbounds for k in data.ions
        sm = 0.0
        for j in data.neutrals
            sm += Xw[j] / (bin[k, j] * qe_over_kBT)
        end
        mob[k] = 1.0 / sm / P
    end

    p_over_kBT = P / kBT
    sigma = 0.0
    @inbounds for k in data.ions
        sigma += Xw[k] * p_over_kBT * abs(data.charges[k]) * _ION_QE * mob[k]
    end
    if e != 0
        @inbounds sigma += Xw[e] * p_over_kBT * _ION_QE * 0.4
    end

    w.pressure = P
    w.temperature = T
    w.mean_MW = mean_MW

    return (vismix, lambda, sigma)
end

"""
    ionized_flux!(flux, workspace, data, gradX, Yleft, Yright;
                  density, electric_field=0.0, frozen=false)

Fill mass fluxes [kg/(m² s)] with midpoint transport and mole-fraction gradients.
The density is the left-node density; drift uses the average of endpoint mass
fractions. The electric-field stage corrects neutral fluxes to conserve total
mass. The frozen stage sets charged fluxes to zero and applies the source
neutral correction without renormalizing neutral mass fractions. Endpoint
mass fractions may contain signed, unnormalized Newton trial values; they are
not clipped or normalized by this algebraic flux function.
"""
function ionized_flux!(flux::AbstractVector, w::IonTransportWorkspace, data::IonTransportData,
                       gradX, Yleft, Yright; density, electric_field=0.0, frozen=false)
    (isfinite(w.pressure) && w.pressure > 0.0 &&
     isfinite(w.temperature) && w.temperature > 0.0 &&
     isfinite(w.mean_MW) && w.mean_MW > 0.0) ||
        throw(DomainError((w.pressure, w.temperature, w.mean_MW),
            "ionized_flux! requires finite positive workspace pressure, temperature and mean_MW from a completed midpoint transport call"))
    n = length(data.species_names)
    all(length(v) == n for v in (flux,gradX,Yleft,Yright,w.diffusion,w.mobility)) ||
        throw(DimensionMismatch("one flux, gradient and mass fraction per species is required"))
    all(isfinite,gradX) && isfinite(electric_field) && isfinite(density) && density > 0 ||
        throw(ArgumentError("finite gradients, field and positive density are required"))
    all(isfinite,Yleft) && all(isfinite,Yright) ||
        throw(ArgumentError("finite endpoint mass fractions are required"))
    MW = data.molecular_weights
    q = data.charges
    RT = Arrhenius.R * w.temperature
    if frozen
        sumflux = 0.0
        @inbounds for k in data.neutrals
            f = -(MW[k] * w.pressure / RT) * w.diffusion[k] * gradX[k]
            flux[k] = f
            sumflux -= f
        end
        @inbounds for k in data.neutrals
            flux[k] += sumflux * Yleft[k]
        end
        @inbounds for k in data.ions
            flux[k] = 0.0
        end
        if data.electron != 0
            @inbounds flux[data.electron] = 0.0
        end
    else
        sumflux = 0.0
        @inbounds for k in 1:n
            f = -(MW[k] * w.pressure / RT) * w.diffusion[k] * gradX[k]
            flux[k] = f
            sumflux -= f
        end
        sumion = 0.0
        @inbounds for k in 1:n
            if q[k] != 0.0
                drift = density * 0.5 * (Yleft[k] + Yright[k]) * electric_field * q[k] * w.mobility[k]
                flux[k] += drift
                sumflux -= drift
                sumion += Yleft[k]
            end
        end
        corr = sumflux / (1.0 - sumion)
        @inbounds for k in data.neutrals
            flux[k] += Yleft[k] * corr
        end
    end
    return flux
end

export IonTransportData, IonTransportWorkspace, ionized_transport!, ionized_flux!
