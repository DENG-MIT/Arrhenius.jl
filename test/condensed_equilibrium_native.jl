module CondensedEquilibriumTests
using Arrhenius, ForwardDiff, LinearAlgebra, Test, YAML

# Graphite NASA7 polynomials for C(gr) from the NASA thermodynamic database
# (Gordon & McBride NASA SP-273; McBride, Gordon & Reno NASA TM-4513), with
# density 2.16 g/cm^3. These are solver input fixtures, not reference states.
const GRAPHITE_LOW = [-0.310872072, 4.40353686e-03, 1.90394118e-06, -6.38546966e-09,
                      2.98964248e-12, -108.650794, 1.11382953]
const GRAPHITE_HIGH = [1.45571829, 1.71702216e-03, -6.97562786e-07, 1.35277032e-10,
                       -9.67590652e-15, -695.138814, -8.52583033]
const GRAPHITE_COEFFICIENTS = vcat(1000.0, GRAPHITE_HIGH, GRAPHITE_LOW)

function graphite_phase()
    return StoichiometricCondensedPhase(;
        name="graphite", species="C(gr)", elements=Dict("C" => 1.0),
        molecular_weight=12.011, molar_volume=12.011 / 2160.0,
        thermo_model="NASA7", thermo_coefficients=GRAPHITE_COEFFICIENTS,
        temperature_range=[200.0, 5000.0], reference_pressure=one_atm)
end

"Constant-cp carbon matching the NASA7 high-temperature region at 1000 K,
with optional uniform h0/s0 shifts to move the formation affinity."
function constant_cp_phase(; h0_shift=0.0, s0_shift=0.0)
    a = GRAPHITE_HIGH
    T0 = 1000.0
    cp0 = R * (a[1] + T0 * (a[2] + T0 * (a[3] + T0 * (a[4] + T0 * a[5]))))
    h0 = R * (T0 * (a[1] + T0 * (a[2] / 2 + T0 * (a[3] / 3 + T0 * (a[4] / 4 + T0 * a[5] / 5)))) + a[6])
    s0 = R * (a[1] * log(T0) + T0 * (a[2] + T0 * (a[3] / 2 + T0 * (a[4] / 3 + T0 * a[5] / 4))) + a[7])
    return StoichiometricCondensedPhase(;
        name="carbon-constant-cp", species="C(s)", elements=Dict("C" => 1.0),
        molecular_weight=12.011, molar_volume=12.011 / 2160.0,
        thermo_model="constant-cp",
        thermo_coefficients=[T0, h0 + h0_shift, s0 + s0_shift, cp0],
        temperature_range=[200.0, 5000.0], reference_pressure=one_atm)
end

const CONDENSED_JSON_UNITS = Dict(
    "temperature" => "K", "pressure" => "Pa", "molecular_weight" => "kg/kmol",
    "density" => "kg/m^3", "molar_volume" => "m^3/kmol", "enthalpy" => "J/kmol",
    "entropy" => "J/kmol/K", "heat_capacity" => "J/kmol/K")

function condensed_json_dict(; format="arrhenius-condensed-phase-v1",
                             elements=Dict("C" => 1.0), model="NASA7",
                             coefficients=GRAPHITE_COEFFICIENTS,
                             temperature_range=[200.0, 5000.0],
                             molecular_weight=12.011, density=2160.0,
                             units=CONDENSED_JSON_UNITS)
    return Dict{String,Any}(
        "format" => format, "source_sha256" => repeat("0", 64),
        "phase_name" => "graphite", "species_name" => "C(gr)",
        "elements" => elements, "charge" => 0.0,
        "molecular_weight_kg_per_kmol" => molecular_weight,
        "density_kg_per_m3" => density,
        "molar_volume_m3_per_kmol" => molecular_weight / density,
        "thermo" => Dict{String,Any}(
            "model" => model, "coefficients" => coefficients,
            "temperature_range_K" => temperature_range,
            "reference_pressure_Pa" => 101325.0),
        "units" => units)
end

function write_json(path, payload)
    open(path, "w") do io
        YAML.write(io, payload, "")
    end
    return path
end

@testset "StoichiometricCondensedPhase keyword validation" begin
    phase = graphite_phase()
    @test phase.thermo_type == [1]
    @test size(phase.thermo_coefficients) == (1, 15)
    @test phase.thermo_coefficients[1, 1] == 1000.0
    @test vec(phase.thermo_coefficients[1, 2:8]) == GRAPHITE_HIGH
    @test vec(phase.thermo_coefficients[1, 9:15]) == GRAPHITE_LOW
    @test phase.elements == Dict("C" => 1.0)
    ccp = constant_cp_phase()
    @test ccp.thermo_type == [2]
    @test vec(ccp.thermo_coefficients[1, 1:4]) ≈
        [1000.0, ccp.thermo_coefficients[1, 2], ccp.thermo_coefficients[1, 3],
         ccp.thermo_coefficients[1, 4]]
    unbounded = StoichiometricCondensedPhase(;
        name="constant-cp", species="C(s)", elements=Dict("C"=>1),
        molecular_weight=12.011, molar_volume=12.011/2160.,
        thermo_model="constant-cp", thermo_coefficients=vec(ccp.thermo_coefficients[1,1:4]),
        temperature_range=[0., Inf])
    @test unbounded.temperature_range == [0., Inf]
    @test Arrhenius._condensed_standard_hs(unbounded,700.,one_atm) ==
        Arrhenius._condensed_standard_hs(ccp,700.,one_atm)
    base = (name="graphite", species="C(gr)", elements=Dict("C" => 1.0),
            molecular_weight=12.011, molar_volume=12.011 / 2160.0,
            thermo_model="NASA7", thermo_coefficients=GRAPHITE_COEFFICIENTS,
            temperature_range=[200.0, 5000.0])
    @test_throws ArgumentError StoichiometricCondensedPhase(; base..., elements=Dict("C" => -1.0))
    @test_throws ArgumentError StoichiometricCondensedPhase(; base..., elements=Dict("C" => 0.0))
    @test_throws ArgumentError StoichiometricCondensedPhase(; base..., elements=Dict("C" => "one"))
    @test_throws ArgumentError StoichiometricCondensedPhase(; base..., charge=1.0)
    @test_throws ArgumentError StoichiometricCondensedPhase(; base..., molecular_weight=0.0)
    @test_throws ArgumentError StoichiometricCondensedPhase(; base..., molecular_weight=NaN)
    @test_throws ArgumentError StoichiometricCondensedPhase(; base..., molar_volume=-1e-3)
    @test_throws ArgumentError StoichiometricCondensedPhase(; base..., reference_pressure=0.0)
    @test_throws ArgumentError StoichiometricCondensedPhase(; base..., thermo_model="Shomate")
    @test_throws ArgumentError StoichiometricCondensedPhase(; base..., thermo_coefficients=ones(14))
    @test_throws ArgumentError StoichiometricCondensedPhase(; base..., thermo_coefficients=ones(16))
    @test_throws ArgumentError StoichiometricCondensedPhase(; base...,
        thermo_coefficients=vcat(6000.0, GRAPHITE_HIGH, GRAPHITE_LOW))
    @test_throws ArgumentError StoichiometricCondensedPhase(; base..., temperature_range=[5000.0, 200.0])
    @test_throws ArgumentError StoichiometricCondensedPhase(; base..., temperature_range=[200.0])
    ccp_coeffs = constant_cp_phase().thermo_coefficients
    @test_throws ArgumentError StoichiometricCondensedPhase(; base..., thermo_model="constant-cp",
        thermo_coefficients=[1000.0, 1e5, 1e4, -2e4])
    @test_throws ArgumentError StoichiometricCondensedPhase(; base..., thermo_model="constant-cp",
        thermo_coefficients=[1000.0, 1e5, 1e4])
end

@testset "condensed-phase JSON round trip and guards" begin
    mktempdir() do directory
        direct = graphite_phase()
        loaded = StoichiometricCondensedPhase(
            write_json(joinpath(directory, "graphite.json"), condensed_json_dict()))
        @test loaded.name == direct.name
        @test loaded.species == direct.species
        @test loaded.elements == direct.elements
        @test loaded.molecular_weight == direct.molecular_weight
        @test loaded.molar_volume == direct.molar_volume
        @test loaded.reference_pressure == direct.reference_pressure
        @test loaded.temperature_range == direct.temperature_range
        open_range = condensed_json_dict(;model="constant-cp",
            coefficients=[298.15,0.,0.,2e4],temperature_range=[0.,nothing])
        unbounded = StoichiometricCondensedPhase(write_json(joinpath(directory,"unbounded.json"),open_range))
        @test unbounded.temperature_range == [0.,Inf]
        @test loaded.thermo_type == direct.thermo_type
        @test loaded.thermo_coefficients == direct.thermo_coefficients
        case(name, payload) = write_json(joinpath(directory, name), payload)
        @test_throws ArgumentError StoichiometricCondensedPhase(
            case("format.json", condensed_json_dict(; format="arrhenius-condensed-phase-v0")))
        missing_elements = condensed_json_dict()
        delete!(missing_elements, "elements")
        @test_throws ArgumentError StoichiometricCondensedPhase(case("missing.json", missing_elements))
        bad_weight = condensed_json_dict()
        bad_weight["molecular_weight_kg_per_kmol"] = "twelve"
        @test_throws ArgumentError StoichiometricCondensedPhase(case("string.json", bad_weight))
        @test_throws ArgumentError StoichiometricCondensedPhase(
            case("negative.json", condensed_json_dict(; elements=Dict("C" => -1.0))))
        @test_throws ArgumentError StoichiometricCondensedPhase(
            case("short.json", condensed_json_dict(; coefficients=ones(14))))
        @test_throws ArgumentError StoichiometricCondensedPhase(
            case("model.json", condensed_json_dict(; model="Shomate", coefficients=ones(15))))
        bad_units = merge(Dict(CONDENSED_JSON_UNITS), Dict("pressure" => "bar"))
        @test_throws ArgumentError StoichiometricCondensedPhase(
            case("units.json", condensed_json_dict(; units=bad_units)))
        @test_throws ArgumentError StoichiometricCondensedPhase(
            case("range.json", condensed_json_dict(; temperature_range=[5000.0, 200.0])))
        inconsistent = condensed_json_dict()
        inconsistent["molar_volume_m3_per_kmol"] = 2 * inconsistent["molar_volume_m3_per_kmol"]
        @test_throws ArgumentError StoichiometricCondensedPhase(case("volume.json", inconsistent))
        notanobject = write_json(joinpath(directory, "text.json"), "not a phase")
        @test_throws ArgumentError StoichiometricCondensedPhase(notanobject)
        broken = joinpath(directory, "broken.json")
        write(broken, "{unclosed")
        @test_throws ArgumentError StoichiometricCondensedPhase(broken)
    end
end

# Independent reference residual, written directly from the defining equations
# with explicit loops so it shares no code path with the workspace kernel:
#   v = Aᵀλ − g − logP;  x = softmax(v);  N = exp(ν)
#   t = N·A·x + As·ξ;   F = [log(t) − log(b); logsumexp(v); Asᵀλ − gp]
function reference_residual(A, As, g, gp, b, logpressure, u)
    ne, ns = size(A)
    np = size(As, 2)
    T = eltype(u)
    λ = u[1:ne]
    ν = u[ne+1]
    ξ = u[ne+2:ne+1+np]
    v = Vector{T}(undef, ns)
    for k in 1:ns
        acc = zero(T)
        for j in 1:ne
            acc += A[j, k] * λ[j]
        end
        v[k] = acc - g[k] - logpressure
    end
    vmax = maximum(v)
    denom = zero(T)
    for k in 1:ns
        denom += exp(v[k] - vmax)
    end
    x = [exp(v[k] - vmax) / denom for k in 1:ns]
    N = exp(ν)
    t = Vector{T}(undef, ne)
    for i in 1:ne
        acc = zero(T)
        for k in 1:ns
            acc += A[i, k] * x[k]
        end
        acc *= N
        for p in 1:np
            acc += As[i, p] * ξ[p]
        end
        t[i] = acc
    end
    F = Vector{T}(undef, ne + 1 + np)
    for i in 1:ne
        F[i] = log(t[i]) - log(b[i])
    end
    F[ne+1] = vmax + log(denom)
    for p in 1:np
        acc = zero(T)
        for j in 1:ne
            acc += As[j, p] * λ[j]
        end
        F[ne+1+p] = acc - gp[p]
    end
    return F
end

function reference_pieces(A, g, logpressure, λ)
    v = transpose(A) * λ .- g .- logpressure
    vmax = maximum(v)
    w = exp.(v .- vmax)
    x = w ./ sum(w)
    return x, A * x
end

@testset "condensed residual and corrected Jacobian" begin
    A1 = [1.0 0.0 1.0 2.0 0.0 1.0;
          0.0 1.0 0.0 1.0 2.0 0.0;
          1.0 0.0 2.0 0.0 1.0 1.0]
    As1 = [1.0 0.0;
           0.0 1.0;
           0.0 2.0]
    g1 = [2.0, -5.0, 8.0, -12.0, 3.5, -7.0]
    gp1 = [-1.0, 2.5]
    b1 = [1.3, 2.1, 0.8]
    logP1 = log(1.5)
    u1 = [-1.2, 0.7, -2.0, log(2.5), 0.4, 0.7]
    Fref1(uu) = reference_residual(A1, As1, g1, gp1, b1, logP1, uu)
    ws = Arrhenius.CondensedWorkspace(A1, As1, b1, g1, gp1, logP1)
    F, J = Arrhenius._condensed_residual_jacobian!(ws, u1)
    @test isapprox(F, Fref1(u1); atol=1e-10)
    @test isapprox(J, ForwardDiff.jacobian(Fref1, u1); atol=1e-10, rtol=1e-8)
    @test Arrhenius._condensed_residual!(ws, u1) ≈ F atol=1e-14
    @test cond(J) < 1e10
    # Normalization row and ν column: ā without an N factor, zeros elsewhere;
    # ∂F_i/∂ν = G_i/t_i strictly below one once solids carry elements.
    x, abar = reference_pieces(A1, g1, logP1, u1[1:3])
    t = exp(u1[4]) .* abar .+ As1 * u1[5:6]
    @test J[4, 1:3] ≈ abar atol=1e-12
    @test J[4, 4] == 0.0
    @test all(==(0.0), J[4, 5:6])
    @test J[1:3, 4] ≈ (exp(u1[4]) .* abar) ./ t atol=1e-12
    @test all(<(1), J[1:3, 4])

    # Gas-only limit (np = 0) reduces exactly to the moment form of the
    # gas-only equilibrium Jacobian.
    A2 = [1.0 0.0 2.0 0.0 1.0;
          0.0 1.0 0.0 1.0 0.0;
          1.0 0.0 0.0 2.0 1.0]
    g2 = [1.0, -3.0, 4.0, -8.0, 2.0]
    b2 = [2.0, 0.5, 1.1]
    u2 = [-0.8, 1.1, -1.5, log(0.7)]
    ws2 = Arrhenius.CondensedWorkspace(A2, zeros(3, 0), b2, g2, Float64[], 0.0)
    F2, J2 = Arrhenius._condensed_residual_jacobian!(ws2, u2)
    @test isapprox(J2, ForwardDiff.jacobian(uu -> reference_residual(
        A2, zeros(3, 0), g2, Float64[], b2, 0.0, uu), u2); atol=1e-10, rtol=1e-8)
    x2, abar2 = reference_pieces(A2, g2, 0.0, u2[1:3])
    Jgas = zeros(4, 4)
    for j in 1:3, i in 1:3
        moment = sum(A2[i, k] * A2[j, k] * x2[k] for k in 1:5)
        Jgas[i, j] = moment / abar2[i] - abar2[j]
    end
    for i in 1:3
        Jgas[i, 4] = 1.0
        Jgas[4, i] = abar2[i]
    end
    @test isapprox(J2, Jgas; atol=1e-12)

    # Boundary validation: zero amounts accepted, negative rejected, never clipped.
    @test Arrhenius._condensed_valid_state(ws, u1)
    @test Arrhenius._condensed_valid_state(ws, [-1.2, 0.7, -2.0, log(2.5), 0.0, 0.0])
    @test_throws ArgumentError Arrhenius._condensed_valid_state(
        ws, [-1.2, 0.7, -2.0, log(2.5), -1e-3, 0.0])
    @test_throws ArgumentError Arrhenius.CondensedWorkspace(A1, As1, [1.0, 0.0, 0.5], g1, gp1, 0.0)
    @test_throws ArgumentError Arrhenius.CondensedWorkspace(A1, As1, [1.0, -2.0, 0.5], g1, gp1, 0.0)
    @test_throws DimensionMismatch Arrhenius.CondensedWorkspace(A1, As1[1:2, :], b1, g1, gp1, 0.0)
    @test_throws DimensionMismatch Arrhenius.CondensedWorkspace(A1, As1, b1, g1[1:5], gp1, 0.0)

    # Affinity sign convention: F_p = a_pᵀλ − g_p is zero when present,
    # negative when formation is unfavorable, positive when favorable.
    λ = u1[1:3]
    a1 = As1[:, 1]
    for (shift, expected) in ((0.0, 0.0), (1.0, -1.0), (-1.0, 1.0))
        wss = Arrhenius.CondensedWorkspace(A1, As1, b1, g1, [dot(a1, λ) + shift, gp1[2]], logP1)
        @test Arrhenius._condensed_residual!(wss, u1)[5] ≈ expected
    end
end

@testset "condensed-phase standard thermodynamics" begin
    function nasa7_hs(a, T)
        h = R * (T * (a[1] + T * (a[2] / 2 + T * (a[3] / 3 + T * (a[4] / 4 + T * a[5] / 5)))) + a[6])
        s = R * (a[1] * log(T) + T * (a[2] + T * (a[3] / 2 + T * (a[4] / 3 + T * a[5] / 4))) + a[7])
        return h, s
    end
    phase = graphite_phase()
    T = 963.981435382116
    href, sref = nasa7_hs(GRAPHITE_LOW, T)
    h, s = Arrhenius._condensed_standard_hs(phase, T, one_atm)
    @test h ≈ href
    @test s ≈ sref
    href2, sref2 = nasa7_hs(GRAPHITE_HIGH, 1200.0)
    h2, s2 = Arrhenius._condensed_standard_hs(phase, 1200.0, one_atm)
    @test h2 ≈ href2
    @test s2 ≈ sref2
    # Incompressible pressure correction h = h0 + (P − Pref)·vm with vm in m^3/kmol.
    vm = phase.molar_volume
    for P in (0.5 * one_atm, 10.0 * one_atm)
        hp, sp = Arrhenius._condensed_standard_hs(phase, T, P)
        @test hp ≈ href + (P - one_atm) * vm
        @test sp ≈ sref
    end
    @test ForwardDiff.derivative(
        P -> Arrhenius._condensed_standard_hs(phase, T, P)[1], one_atm) ≈ vm
    @test ForwardDiff.derivative(
        P -> Arrhenius._condensed_gibbs(phase, T, P), one_atm) ≈ vm / (R * T)
    # Constant-cp: h = h0 + cp0·(T − T0), s = s0 + cp0·log(T/T0).
    ccp = constant_cp_phase()
    T0, h0, s0, cp0 = ccp.thermo_coefficients[1, 1:4]
    h3, s3 = Arrhenius._condensed_standard_hs(ccp, T0, one_atm)
    @test h3 ≈ h0
    @test s3 ≈ s0
    h4, s4 = Arrhenius._condensed_standard_hs(ccp, 700.0, one_atm)
    @test h4 ≈ h0 + cp0 * (700.0 - T0)
    @test s4 ≈ s0 + cp0 * log(700.0 / T0)
    @test ForwardDiff.derivative(
        P -> Arrhenius._condensed_standard_hs(ccp, 700.0, P)[1], one_atm) ≈ vm
    @test ForwardDiff.derivative(
        P -> Arrhenius._condensed_gibbs(ccp, 700.0, P), one_atm) ≈ vm / (R * 700.0)
end

@testset "TP equilibrium with an initially absent graphite phase" begin
    gas = CreateSolution(joinpath(@__DIR__, "..", "mechanism", "gri30.yaml"))
    phase = graphite_phase()
    X = set_equivalence_ratio(gas, 3.5; fuel="CH4", oxidizer="O2:1,N2:3.76")
    # Independently regenerated Cantera adiabatic graphite reference point
    # (validation harness, agreement within 3.3e-7 K); assertion only, never
    # used as a solver initial state.
    Tref, nref = 963.981435382116, 0.021389703536144382
    result = equilibrate(gas, phase; T=Tref, P=one_atm, X, mode=:TP)
    @test result.status == :condensed_present
    @test result.T == Tref
    @test result.P == one_atm
    @test result.condensed_moles ≈ nref rtol=1e-6
    @test result.condensed_moles > 0
    @test result.affinity ≈ 0 atol=1e-10
    @test abs(result.condensed_moles * result.affinity) < 1e-10
    @test result.element_error < 2e-9
    @test all(>=(0), result.gas_species_moles)
    @test length(result.gas_species_moles) == gas.n_species
    @test length(result.species_moles) == gas.n_species + 1
    @test result.species_moles[1:end-1] == result.gas_species_moles
    @test result.species_moles[end] == result.condensed_moles
    @test sum(result.gas_species_moles) ≈ result.gas_moles
    @test result.X ≈ result.gas_species_moles ./ result.gas_moles
    @test sum(result.X) ≈ 1 atol=1e-14
    @test sum(result.Y) ≈ 1 atol=1e-14
    # Complete element balance over every original row of the gas mechanism.
    before = gas.ele_matrix * X
    after = gas.ele_matrix * result.gas_species_moles
    after[findfirst(==("C"), gas.elements)] += result.condensed_moles
    @test maximum(abs.(after - before)) <= 2e-9 * maximum(abs, before)
    # Carbon moved into the condensed phase relative to the gas-only root.
    gas_only = equilibrate(gas; T=Tref, P=one_atm, X)
    @test !(isapprox(result.X, gas_only.X; atol=1e-12))
    # Amount scaling: the per-kmol solve scales consistently, temperatures and
    # compositions do not depend on initial_gas_moles.
    scaled = equilibrate(gas, phase; T=Tref, P=one_atm, X, mode=:TP, initial_gas_moles=2.5)
    @test scaled.T == result.T
    @test scaled.X ≈ result.X atol=1e-14
    @test scaled.gas_species_moles == 2.5 .* result.gas_species_moles
    @test scaled.condensed_moles == 2.5 * result.condensed_moles
    @test scaled.gas_moles == 2.5 * result.gas_moles
    @test scaled.species_moles == 2.5 .* result.species_moles
    tiny = equilibrate(gas, phase; T=Tref, P=one_atm, X, mode=:TP, initial_gas_moles=1e-300)
    @test tiny.X == result.X
    @test tiny.Y == result.Y

    # Lean mixture: graphite stays absent and the gas-only state is returned.
    Xlean = set_equivalence_ratio(gas, 0.5; fuel="CH4", oxidizer="O2:1,N2:3.76")
    lean = equilibrate(gas, phase; T=1500.0, P=one_atm, X=Xlean)
    gas_lean = equilibrate(gas; T=1500.0, P=one_atm, X=Xlean)
    @test lean.status == :gas_only
    @test lean.condensed_moles == 0
    @test lean.affinity < 0
    @test lean.T == 1500.0
    @test lean.X ≈ gas_lean.X atol=1e-13
    @test lean.gas_species_moles ≈ lean.gas_moles .* gas_lean.X atol=1e-13

    # A constant-cp phase with shifted standard enthalpy: strongly unfavorable
    # stays absent, strongly favorable forms with all invariants intact.
    unfavorable = constant_cp_phase(; h0_shift=2e8)
    away = equilibrate(gas, unfavorable; T=Tref, P=one_atm, X)
    @test away.status == :gas_only
    @test away.condensed_moles == 0
    @test away.X ≈ gas_only.X atol=1e-13
    favorable = constant_cp_phase(; s0_shift=2 * R)
    formed = equilibrate(gas, favorable; T=Tref, P=one_atm, X)
    @test formed.status == :condensed_present
    @test formed.condensed_moles > 0
    @test formed.affinity ≈ 0 atol=1e-10
    @test formed.element_error < 2e-9
    @test all(>=(0), formed.gas_species_moles)
end

@testset "HP equilibrium with an initially absent graphite phase" begin
    gas = CreateSolution(joinpath(@__DIR__, "..", "mechanism", "gri30.yaml"))
    phase = graphite_phase()
    X = set_equivalence_ratio(gas, 3.5; fuel="CH4", oxidizer="O2:1,N2:3.76")
    # Same independently regenerated Cantera reference point as the TP testset.
    Tref, nref = 963.981435382116, 0.021389703536144382
    rich = equilibrate(gas, phase; T=300.0, P=one_atm, X, mode=:HP)
    @test rich.status == :condensed_present
    @test rich.T ≈ Tref atol=1e-5
    @test rich.P == one_atm
    @test rich.condensed_moles ≈ nref rtol=1e-6
    @test rich.enthalpy_error < 1e-9
    @test rich.element_error < 2e-9
    @test abs(rich.condensed_moles * rich.affinity) < 1e-10
    @test all(>=(0), rich.gas_species_moles)
    @test all(>=(0), rich.X)
    # Complete enthalpy conservation against the 300 K, 1 atm inlet, reassembled
    # from the returned amounts and public thermodynamics only.
    hinlet = dot(X, cal_h_RT(gas, 300.0, one_atm, X)) * (R * 300.0)
    hgas = cal_h_RT(gas, rich.T, one_atm, rich.X) .* (R * rich.T)
    hphase, = Arrhenius._condensed_standard_hs(phase, rich.T, rich.P)
    hout = dot(rich.gas_species_moles, hgas) + rich.condensed_moles * hphase
    @test abs(hout - hinlet) <= 1e-9 * max(abs(hinlet), 1e6 * dot(gas.MW, X))
    # Complete element balance over every original row.
    before = gas.ele_matrix * X
    after = gas.ele_matrix * rich.gas_species_moles
    after[findfirst(==("C"), gas.elements)] += rich.condensed_moles
    @test maximum(abs.(after - before)) <= 2e-9 * maximum(abs, before)
    # Amount scaling.
    scaled = equilibrate(gas, phase; T=300.0, P=one_atm, X, mode=:HP, initial_gas_moles=4.0)
    @test scaled.T == rich.T
    @test scaled.X ≈ rich.X atol=1e-14
    @test scaled.gas_species_moles == 4.0 .* rich.gas_species_moles
    @test scaled.condensed_moles == 4.0 * rich.condensed_moles
    @test scaled.species_moles == 4.0 .* rich.species_moles

    # Below the graphite onset the adiabatic state is the gas-only one.
    Xlean = set_equivalence_ratio(gas, 1.0; fuel="CH4", oxidizer="O2:1,N2:3.76")
    lean = equilibrate(gas, phase; T=300.0, P=one_atm, X=Xlean, mode=:HP)
    gas_lean = equilibrate(gas; T=300.0, P=one_atm, X=Xlean, mode=:HP)
    @test lean.status == :gas_only
    @test lean.condensed_moles == 0
    @test lean.T ≈ gas_lean.T atol=1e-6
    @test lean.X ≈ gas_lean.X atol=1e-12
    Xmid = set_equivalence_ratio(gas, 3.0; fuel="CH4", oxidizer="O2:1,N2:3.76")
    mid = equilibrate(gas, phase; T=300.0, P=one_atm, X=Xmid, mode=:HP)
    @test mid.status == :gas_only
    @test mid.condensed_moles == 0
end

@testset "conservation-driven inactivity and input rejection" begin
    gas = CreateSolution(joinpath(@__DIR__, "..", "mechanism", "gri30.yaml"))
    phase = graphite_phase()
    # Zero initial carbon inventory: graphite formation is excluded by
    # conservation in both modes; no affinity is invented for the removed row.
    X = set_equivalence_ratio(gas, 1.0; fuel="H2", oxidizer="O2:1,N2:3.76")
    for mode in (:TP, :HP)
        result = equilibrate(gas, phase; T=300.0, P=one_atm, X, mode)
        gas_only = equilibrate(gas; T=300.0, P=one_atm, X, mode)
        @test result.status == :inactive_by_elements
        @test result.condensed_moles == 0
        @test isnan(result.affinity)
        @test result.T ≈ gas_only.T atol=1e-6
        @test result.X ≈ gas_only.X atol=1e-10
    end
    # Elements absent from the gas mechanism are rejected.
    fe_phase = StoichiometricCondensedPhase(;
        name="iron", species="Fe(s)", elements=Dict("Fe" => 1.0),
        molecular_weight=55.845, molar_volume=55.845 / 7874.0,
        thermo_model="constant-cp", thermo_coefficients=[298.15, 0.0, 27280.0, 25100.0],
        temperature_range=[200.0, 6000.0])
    Xc = mole_fractions(gas, "CH4")
    @test_throws ArgumentError equilibrate(gas, fe_phase; T=1000.0, P=one_atm, X=Xc)
    # Unsupported modes and invalid controls.
    @test_throws ArgumentError equilibrate(gas, phase; T=1000.0, X=Xc, mode=:UV)
    @test_throws ArgumentError equilibrate(gas, phase; T=1000.0, X=Xc, mode=:TV)
    @test_throws ArgumentError equilibrate(gas, phase; T=1000.0, X=Xc, mode=:SP)
    @test_throws ArgumentError equilibrate(gas, phase; T=1000.0, X=Xc, mode="SV")
    @test_throws ArgumentError equilibrate(gas, phase; T=1000.0, X=Xc, initial_gas_moles=0.0)
    @test_throws ArgumentError equilibrate(gas, phase; T=1000.0, X=Xc, initial_gas_moles=-1.0)
    @test_throws ArgumentError equilibrate(gas, phase; T=1000.0, X=Xc, initial_gas_moles=NaN)
    @test_throws ArgumentError equilibrate(gas, phase; T=-1.0, X=Xc)
    @test_throws ArgumentError equilibrate(gas, phase; T=1000.0, P=0.0, X=Xc)
    @test_throws ArgumentError equilibrate(gas, phase; T=1000.0, X=Xc, property_rtol=0.0)
    @test_throws ArgumentError equilibrate(gas, phase; T=1000.0, X=Xc, mode=:HP,
        temperature_bounds=(500.0, 100.0))
end

@testset "dependent element rows" begin
    # A single-species H2O gas has dependent H/O rows. A condensed composition
    # that cannot satisfy the removed row is rejected; a consistent one stays
    # admissible.
    gas = CreateSolution(joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml"))
    inds = [findfirst(==("H2O"), gas.species_names)]
    thermo = Arrhenius.IdealGasThermo(gas.thermo.nasa_low[inds, :],
        gas.thermo.nasa_high[inds, :], gas.thermo.Trange[inds, :], true)
    single = Arrhenius.Solution(1, gas.n_reactions, gas.MW[inds],
        gas.species_names[inds], gas.elements, gas.ele_matrix[:, inds],
        thermo, gas.trans, gas.reaction)
    hydrogen = StoichiometricCondensedPhase(;
        name="hydrogen-solid", species="H2(s)", elements=Dict("H" => 2.0),
        molecular_weight=2.016, molar_volume=2.016 / 70.0,
        thermo_model="constant-cp", thermo_coefficients=[298.15, 0.0, 1.3e5, 2.9e4],
        temperature_range=[200.0, 6000.0])
    @test_throws ArgumentError equilibrate(single, hydrogen; T=800.0, P=one_atm, X=[1.0])
    ice = StoichiometricCondensedPhase(;
        name="ice", species="H2O(s)", elements=Dict("H" => 2.0, "O" => 1.0),
        molecular_weight=18.015, molar_volume=18.015 / 917.0,
        thermo_model="constant-cp", thermo_coefficients=[298.15, 0.0, 7.0e4, 3.8e4],
        temperature_range=[200.0, 6000.0])
    result = equilibrate(single, ice; T=800.0, P=one_atm, X=[1.0])
    @test result.status == :gas_only
    @test result.condensed_moles == 0
    @test result.X == [1.0]
    @test result.element_error < 2e-9
end

end # module CondensedEquilibriumTests
