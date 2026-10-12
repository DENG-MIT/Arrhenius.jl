# Native Julia adaptation of Cantera's coverage-dependent surface thermo
# (CoverageDependentSurfPhase.cpp/.h), applying adsorbate lateral-interaction
# correction terms to ideal surface species standard-state properties. The
# four enthalpy/entropy dependency laws (linear, polynomial, piecewise-linear,
# interpolative) and the log-quadratic coverage heat-capacity law with its
# 298.15 K-anchored integrals follow the Cantera 4.0 source.
# https://github.com/Cantera/cantera/blob/726522be4e2a13454d8415b7ef799d621f665cf3/src/thermo/CoverageDependentSurfPhase.cpp
#
# Copyright (c) 2001-2009, California Institute of Technology
# All rights reserved.
# Copyright (c) 2009 Sandia Corporation. Under the terms of Contract
# AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain
# rights in this software.
# Copyright (c) 2011-2026, Cantera Developers. All rights reserved.
#
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

# All quantities use the Arrhenius.jl SI convention: K, kmol, J. Dependency
# enthalpy parameters are in J/kmol, entropy and heat-capacity parameters in
# J/kmol/K, and coverages are dimensionless site fractions. This component is
# pure standard-state thermodynamics; it is deliberately decoupled from the
# ideal-surface kinetics in SurfaceKinetics.jl and asserts no chemistry
# coupling. The `+log(reference_coverage)` shift applied to the standard-state
# entropy (and the opposite shift to Gibbs energy) is the standard-state
# reference-coverage correction; it is distinct from the partial-molar mixing
# entropy `-R*log(theta_k/theta_ref)` that Cantera adds only in
# getPartialMolarEntropies/getChemPotentials, which this component does not
# compute.

"Supertype for composable coverage-dependency laws targeting species standard states."
abstract type AbstractCoverageDependency end

"""
    LinearDependency(target, influencing; enthalpy=0.0, entropy=0.0)

Linear lateral-interaction law. Adds `enthalpy*theta_j` [J/kmol] and
`entropy*theta_j` [J/kmol/K] to the standard state of species `target`,
where `j = influencing`.
"""
struct LinearDependency <: AbstractCoverageDependency
    target::Int
    influencing::Int
    enthalpy::Float64
    entropy::Float64
end
function LinearDependency(target, influencing; enthalpy=0.0, entropy=0.0)
    all(isfinite, (enthalpy, entropy)) ||
        throw(ArgumentError("finite linear dependency parameters required"))
    return LinearDependency(target, influencing, enthalpy, entropy)
end

"""
    PolynomialDependency(target, influencing; enthalpy_coefficients=zeros(4),
                         entropy_coefficients=zeros(4))

Quartic polynomial law `sum_i c_i*theta_j^i` for orders one through four.
Enthalpy coefficients are in J/kmol, entropy coefficients in J/kmol/K.
"""
struct PolynomialDependency <: AbstractCoverageDependency
    target::Int
    influencing::Int
    enthalpy_coefficients::NTuple{4,Float64}
    entropy_coefficients::NTuple{4,Float64}
end
function PolynomialDependency(target, influencing; enthalpy_coefficients=zeros(4),
                              entropy_coefficients=zeros(4))
    length(enthalpy_coefficients) == 4 && length(entropy_coefficients) == 4 ||
        throw(ArgumentError("four polynomial coefficients per property required"))
    all(isfinite, enthalpy_coefficients) && all(isfinite, entropy_coefficients) ||
        throw(ArgumentError("finite polynomial coefficients required"))
    return PolynomialDependency(target, influencing,
        NTuple{4,Float64}(enthalpy_coefficients), NTuple{4,Float64}(entropy_coefficients))
end

"""
    PiecewiseLinearDependency(target, influencing; enthalpy_low, enthalpy_high,
                              enthalpy_change, entropy_low, entropy_high,
                              entropy_change)

Two-segment law with slopes `*_low` below the change coverage and `*_high`
above it, continuous at `*_change`. Slopes are in J/kmol (enthalpy) or
J/kmol/K (entropy); change coverages are dimensionless in (0, 1).
"""
struct PiecewiseLinearDependency <: AbstractCoverageDependency
    target::Int
    influencing::Int
    enthalpy_low::Float64
    enthalpy_high::Float64
    enthalpy_change::Float64
    entropy_low::Float64
    entropy_high::Float64
    entropy_change::Float64
end
function PiecewiseLinearDependency(target, influencing; enthalpy_low, enthalpy_high,
                                   enthalpy_change, entropy_low, entropy_high,
                                   entropy_change)
    all(isfinite, (enthalpy_low, enthalpy_high, entropy_low, entropy_high)) ||
        throw(ArgumentError("finite piecewise-linear slopes required"))
    0 < enthalpy_change < 1 && 0 < entropy_change < 1 ||
        throw(ArgumentError("piecewise-linear change coverages must lie in (0, 1)"))
    return PiecewiseLinearDependency(target, influencing, enthalpy_low,
        enthalpy_high, enthalpy_change, entropy_low, entropy_high, entropy_change)
end

"""
    InterpolativeDependency(target, influencing; enthalpy_coverages, enthalpies,
                            entropy_coverages, entropies)

Tabulated law linearly interpolated between knots. Coverage knots must be
strictly increasing, start at 0.0, end at 1.0, and match the value lengths.
Enthalpies are in J/kmol, entropies in J/kmol/K.
"""
struct InterpolativeDependency <: AbstractCoverageDependency
    target::Int
    influencing::Int
    enthalpy_coverages::Vector{Float64}
    enthalpies::Vector{Float64}
    entropy_coverages::Vector{Float64}
    entropies::Vector{Float64}
end
function InterpolativeDependency(target, influencing; enthalpy_coverages,
                                 enthalpies, entropy_coverages, entropies)
    _valid_knots(enthalpy_coverages, enthalpies, "enthalpy")
    _valid_knots(entropy_coverages, entropies, "entropy")
    return InterpolativeDependency(target, influencing,
        collect(Float64, enthalpy_coverages), collect(Float64, enthalpies),
        collect(Float64, entropy_coverages), collect(Float64, entropies))
end
function _valid_knots(coverages, values, label)
    length(coverages) == length(values) && length(coverages) >= 2 ||
        throw(ArgumentError("$label coverage/value tables need matching lengths >= 2"))
    all(isfinite, coverages) && all(isfinite, values) ||
        throw(ArgumentError("finite $label interpolation tables required"))
    coverages[1] == 0.0 && coverages[end] == 1.0 ||
        throw(ArgumentError("$label coverages must start at 0.0 and end at 1.0"))
    all(diff(collect(coverages)) .> 0) ||
        throw(ArgumentError("$label coverages must be strictly increasing"))
    return nothing
end

"""
    HeatCapacityDependency(target, influencing; a, b)

Log-quadratic coverage heat capacity `(a*log(T) + b)*theta_j^2` [J/kmol/K].
Its contributions to enthalpy and entropy are integrated analytically and
anchored to zero at 298.15 K, following Cantera.
"""
struct HeatCapacityDependency <: AbstractCoverageDependency
    target::Int
    influencing::Int
    a::Float64
    b::Float64
end
function HeatCapacityDependency(target, influencing; a, b)
    all(isfinite, (a, b)) ||
        throw(ArgumentError("finite heat-capacity coefficients required"))
    return HeatCapacityDependency(target, influencing, a, b)
end

"""
    CoverageThermoModel(; species_names, thermo_type, thermo_coefficients,
                        reference_coverage=1.0,
                        dependencies=AbstractCoverageDependency[])

Coverage-dependent surface standard-state thermodynamics. `thermo_type` uses
the SurfaceMechanism convention (1 = NASA7, 2 = constant-cp) and
`thermo_coefficients` uses its per-species row layout, so the base
low-coverage states are evaluated with the existing `_surface_species_thermo`.
`reference_coverage` is the standard-state reference coverage theta_ref in
(0, 1]; `dependencies` lists self- and cross-interaction objects whose target
and influencing indices refer to `species_names`.
"""
struct CoverageThermoModel
    species_names::Vector{String}
    thermo_type::Vector{Int}
    thermo_coefficients::Matrix{Float64}
    reference_coverage::Float64
    dependencies::Vector{AbstractCoverageDependency}
end
function CoverageThermoModel(; species_names, thermo_type, thermo_coefficients,
                             reference_coverage=1.0,
                             dependencies=AbstractCoverageDependency[])
    names = String.(collect(species_names))
    n = length(names)
    n >= 1 && all(!isempty, names) ||
        throw(ArgumentError("at least one nonempty species name required"))
    length(unique(names)) == n || throw(ArgumentError("duplicate species names"))
    types = Int.(collect(thermo_type))
    length(types) == n && all(t -> t in (1, 2), types) ||
        throw(ArgumentError("one thermo_type (1=NASA7, 2=constant-cp) per species required"))
    coeffs = Matrix{Float64}(thermo_coefficients)
    size(coeffs, 1) == n || throw(DimensionMismatch("one thermo coefficient row per species"))
    ncols = any(==(1), types) ? 15 : 4
    size(coeffs, 2) >= ncols ||
        throw(DimensionMismatch("thermo coefficient layout requires >= $ncols columns"))
    all(isfinite, coeffs) || throw(ArgumentError("finite thermo coefficients required"))
    all(>(0), @view(coeffs[:, 1])) ||
        throw(ArgumentError("positive NASA7 split or constant-cp reference temperature required"))
    isfinite(reference_coverage) && 0 < reference_coverage <= 1 ||
        throw(ArgumentError("reference coverage must lie in (0, 1]"))
    deps = AbstractCoverageDependency[dependencies...]
    for d in deps
        1 <= d.target <= n && 1 <= d.influencing <= n ||
            throw(ArgumentError("dependency species index out of range"))
        _validate_coverage_dependency(d)
    end
    return CoverageThermoModel(names, types, coeffs, reference_coverage, deps)
end
_coverage_species_count(m::CoverageThermoModel) = length(m.species_names)

# Validate again at model construction so positional struct constructors cannot
# introduce an invalid physical law through the typed dependency collection.
_validate_coverage_dependency(d::LinearDependency) =
    LinearDependency(d.target, d.influencing; enthalpy=d.enthalpy, entropy=d.entropy)
_validate_coverage_dependency(d::PolynomialDependency) =
    PolynomialDependency(d.target, d.influencing;
        enthalpy_coefficients=d.enthalpy_coefficients, entropy_coefficients=d.entropy_coefficients)
_validate_coverage_dependency(d::PiecewiseLinearDependency) =
    PiecewiseLinearDependency(d.target, d.influencing; enthalpy_low=d.enthalpy_low,
        enthalpy_high=d.enthalpy_high, enthalpy_change=d.enthalpy_change,
        entropy_low=d.entropy_low, entropy_high=d.entropy_high, entropy_change=d.entropy_change)
function _validate_coverage_dependency(d::InterpolativeDependency)
    _valid_knots(d.enthalpy_coverages, d.enthalpies, "enthalpy")
    _valid_knots(d.entropy_coverages, d.entropies, "entropy")
end
_validate_coverage_dependency(d::HeatCapacityDependency) =
    HeatCapacityDependency(d.target, d.influencing; a=d.a, b=d.b)

"""
    CoverageThermoModel(path::AbstractString)

Load SI model parameters prepared by `mechanism/export_coverage_thermo.py`.
The parameter archive is JSON; Cantera is needed only for preprocessing.
"""
function CoverageThermoModel(path::AbstractString)
    data = YAML.load_file(path)
    data["format"] == "arrhenius-coverage-thermo-v1" ||
        throw(ArgumentError("unsupported coverage thermodynamics parameter format"))
    data["units"] == "K, J/kmol, J/kmol/K" ||
        throw(ArgumentError("coverage thermodynamics parameters must use SI units"))
    dependencies = AbstractCoverageDependency[]
    for entry in data["dependencies"]
        target, influencing = entry["target"], entry["influencing"]
        if entry["kind"] == "polynomial"
            push!(dependencies, PolynomialDependency(target, influencing;
                enthalpy_coefficients=entry["enthalpy"], entropy_coefficients=entry["entropy"]))
        elseif entry["kind"] == "interpolative"
            h, s = entry["enthalpy"], entry["entropy"]
            push!(dependencies, InterpolativeDependency(target, influencing;
                enthalpy_coverages=h["coverages"], enthalpies=h["values"],
                entropy_coverages=s["coverages"], entropies=s["values"]))
        else
            throw(ArgumentError("unsupported coverage dependency kind: $(entry["kind"])"))
        end
        a, b = entry["heat_capacity_a"], entry["heat_capacity_b"]
        if a != 0 || b != 0
            push!(dependencies, HeatCapacityDependency(target, influencing; a, b))
        end
    end
    coefficients = permutedims(hcat(data["thermo_coefficients"]...))
    return CoverageThermoModel(; species_names=data["species_names"], thermo_type=data["thermo_type"],
        thermo_coefficients=coefficients, reference_coverage=data["reference_state_coverage"], dependencies)
end

"""
    CoverageThermoWorkspace(model)

Reusable per-state output buffers. The arrays returned by
`coverage_thermo!` alias this workspace, so full coverage sweeps evaluate
without per-state allocation.
"""
struct CoverageThermoWorkspace
    enthalpy_RT::Vector{Float64}
    entropy_R::Vector{Float64}
    cp_R::Vector{Float64}
    gibbs_RT::Vector{Float64}
    h_cov::Vector{Float64}
    s_cov::Vector{Float64}
    cp_cov::Vector{Float64}
end
function CoverageThermoWorkspace(m::CoverageThermoModel)
    n = _coverage_species_count(m)
    return CoverageThermoWorkspace(zeros(n), zeros(n), zeros(n), zeros(n),
        zeros(n), zeros(n), zeros(n))
end

# Base constant-pressure heat capacity [J/kmol/K] matching the
# `_surface_species_thermo` coefficient layout (NASA7 / constant-cp).
@inline function _coverage_base_cp(m, k, T)
    c = m.thermo_coefficients
    if m.thermo_type[k] == 1
        j = T <= c[k, 1] ? 9 : 2
        a1, a2, a3, a4, a5 = (c[k, j+i] for i in 0:4)
        return R * (a1 + T*(a2 + T*(a3 + T*(a4 + T*a5))))
    end
    return c[k, 4]
end

@inline _poly4(theta, c) =
    theta*(c[1] + theta*(c[2] + theta*(c[3] + theta*c[4])))

@inline function _interp(coverages, values, theta)
    theta <= coverages[1] && return values[1]
    theta >= coverages[end] && return values[end]
    i = min(searchsortedlast(coverages, theta), length(coverages) - 1)
    return (values[i+1] - values[i]) / (coverages[i+1] - coverages[i]) *
           (theta - coverages[i]) + values[i]
end

function _coverage_correction!(m::CoverageThermoModel, h_cov, s_cov, cp_cov, T, theta)
    fill!(h_cov, 0.0)
    fill!(s_cov, 0.0)
    fill!(cp_cov, 0.0)
    for d in m.dependencies
        _add_dependency!(d, h_cov, s_cov, cp_cov, T, theta)
    end
    return nothing
end

function _add_dependency!(d::LinearDependency, h_cov, s_cov, cp_cov, T, theta)
    q = theta[d.influencing]
    h_cov[d.target] += d.enthalpy * q
    s_cov[d.target] += d.entropy * q
    return nothing
end

function _add_dependency!(d::PolynomialDependency, h_cov, s_cov, cp_cov, T, theta)
    q = theta[d.influencing]
    h_cov[d.target] += _poly4(q, d.enthalpy_coefficients)
    s_cov[d.target] += _poly4(q, d.entropy_coefficients)
    return nothing
end

function _add_dependency!(d::PiecewiseLinearDependency, h_cov, s_cov, cp_cov, T, theta)
    q = theta[d.influencing]
    h_cov[d.target] += q <= d.enthalpy_change ? d.enthalpy_low * q :
        d.enthalpy_low * d.enthalpy_change + d.enthalpy_high * (q - d.enthalpy_change)
    s_cov[d.target] += q <= d.entropy_change ? d.entropy_low * q :
        d.entropy_low * d.entropy_change + d.entropy_high * (q - d.entropy_change)
    return nothing
end

function _add_dependency!(d::InterpolativeDependency, h_cov, s_cov, cp_cov, T, theta)
    q = theta[d.influencing]
    h_cov[d.target] += _interp(d.enthalpy_coverages, d.enthalpies, q)
    s_cov[d.target] += _interp(d.entropy_coverages, d.entropies, q)
    return nothing
end

function _add_dependency!(d::HeatCapacityDependency, h_cov, s_cov, cp_cov, T, theta)
    q2 = theta[d.influencing]^2
    a, b = d.a, d.b
    cp_cov[d.target] += (a*log(T) + b) * q2
    h_cov[d.target] += (T*(a*log(T) - a + b) - 298.15*(a*log(298.15) - a + b)) * q2
    s_cov[d.target] += 0.5 * (log(T)*(a*log(T) + 2b) - log(298.15)*(a*log(298.15) + 2b)) * q2
    return nothing
end

function _check_coverage_state(m::CoverageThermoModel, T, theta)
    n = _coverage_species_count(m)
    length(theta) == n || throw(DimensionMismatch("coverage length"))
    isfinite(T) && T > 0 || throw(DomainError(T, "positive finite temperature required"))
    all(q -> isfinite(q) && 0 <= q <= 1, theta) ||
        throw(DomainError(theta, "coverages must be finite and lie in [0, 1]"))
    abs(sum(theta) - 1) <= 1e-8 ||
        throw(ArgumentError("coverages must sum to one"))
    return nothing
end

"""
    coverage_thermo!(workspace, model, T, theta)

Evaluate coverage-dependent standard-state properties at temperature `T` [K]
and coverages `theta` (must lie in [0, 1] and sum to one). Fills
`workspace.enthalpy_RT` (h°/RT), `workspace.entropy_R` (s°/R, including the
`+log(theta_ref)` standard-state reference-coverage correction),
`workspace.cp_R` (cp°/R), and `workspace.gibbs_RT` (g°/RT = h°/RT - s°/R).
The returned arrays alias the workspace.
"""
function coverage_thermo!(w::CoverageThermoWorkspace, m::CoverageThermoModel, T, theta)
    _check_coverage_state(m, T, theta)
    n = _coverage_species_count(m)
    all(v -> length(v) == n, (w.enthalpy_RT, w.entropy_R, w.cp_R, w.gibbs_RT,
                              w.h_cov, w.s_cov, w.cp_cov)) ||
        throw(DimensionMismatch("workspace and coverage model species counts differ"))
    _coverage_correction!(m, w.h_cov, w.s_cov, w.cp_cov, T, theta)
    RT = R * T
    log_theta_ref = log(m.reference_coverage)
    for k in 1:n
        h0, s0 = _surface_species_thermo(m, k, T)
        cp0 = _coverage_base_cp(m, k, T)
        w.enthalpy_RT[k] = (h0 + w.h_cov[k]) / RT
        w.entropy_R[k] = (s0 + w.s_cov[k]) / R + log_theta_ref
        w.cp_R[k] = (cp0 + w.cp_cov[k]) / R
        w.gibbs_RT[k] = w.enthalpy_RT[k] - w.entropy_R[k]
    end
    return w
end

"""
    coverage_thermo(model, T, theta)

Allocating convenience form of `coverage_thermo!`; returns a NamedTuple of
independently owned arrays.
"""
function coverage_thermo(m::CoverageThermoModel, T, theta)
    w = coverage_thermo!(CoverageThermoWorkspace(m), m, T, theta)
    return (enthalpy_RT=w.enthalpy_RT, entropy_R=w.entropy_R,
        cp_R=w.cp_R, gibbs_RT=w.gibbs_RT)
end

export AbstractCoverageDependency, LinearDependency, PolynomialDependency
export PiecewiseLinearDependency, InterpolativeDependency, HeatCapacityDependency
export CoverageThermoModel, CoverageThermoWorkspace
export coverage_thermo!, coverage_thermo
