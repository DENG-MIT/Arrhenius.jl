using Arrhenius
using ForwardDiff
using LinearAlgebra
using SHA
using Test

include("plasma_eedf_loader.jl")
include("plasma_eedf_numerics.jl")
include("plasma_native.jl")
include("plasma_thermochemistry.jl")
include("plasma_energy.jl")
include("equilibrium_native.jl")
include("signed_equilibrium_native.jl")
include("condensed_equilibrium_native.jl")
include("ideal_gas_states.jl")
include("ideal_gas_mixing_native.jl")
include("species_thermo_native.jl")
include("thermo_means_native.jl")
include("reactors_native.jl")
include("reactor_networks_native.jl")
include("moving_wall_native.jl")
include("network_allocations_native.jl")
include("ic_engine_native.jl")
include("surfaces_native.jl")
include("coverage_thermo_native.jl")
include("surface_flow_native.jl")
include("pure_water_native.jl")
include("critical_properties_native.jl")
include("real_gas_native.jl")
include("real_gas_reactors_native.jl")
include("real_gas_ad_jacobian_native.jl")
include("real_gas_trial_native.jl")
include("reverse_rate_refresh_native.jl")
include("flames_native.jl")
include("conservative_flames_native.jl")
include("flame_derivatives_native.jl")
include("flame_profile_grid_native.jl")
include("counterflow_native.jl")
include("counterflow_premixed.jl")
include("catalytic_flames_native.jl")
include("dusty_gas.jl")
include("kinetics_cache.jl")
include("blowers_masel_native.jl")
include("inert_phase.jl")
include("imported_gas_native.jl")
include("ion_transport_native.jl")
include("ionized_flames_native.jl")
include("ion_conservative_native.jl")

@testset "sidecar provenance" begin
    mechanism, stream = mktemp()
    write(stream, "mechanism fixture\n")
    close(stream)
    digest = bytes2hex(SHA.sha256(read(mechanism)))
    metadata = Dict(
        "sidecar_format_utf8" => collect(codeunits("arrhenius-sidecar-v2")),
        "source_sha256_utf8" => collect(codeunits(digest)),
    )
    @test Arrhenius._validate_sidecar_metadata(metadata, mechanism) === nothing
    metadata_v3 = copy(metadata)
    metadata_v3["sidecar_format_utf8"] = collect(codeunits("arrhenius-sidecar-v3"))
    @test Arrhenius._validate_sidecar_metadata(metadata_v3, mechanism) === nothing

    stale = copy(metadata)
    stale["source_sha256_utf8"] = collect(codeunits(repeat("0", 64)))
    @test_throws ArgumentError Arrhenius._validate_sidecar_metadata(stale, mechanism)

    unknown = copy(metadata)
    unknown["sidecar_format_utf8"] = collect(codeunits("arrhenius-sidecar-v99"))
    @test_throws ArgumentError Arrhenius._validate_sidecar_metadata(unknown, mechanism)

    for source in ("species:\n  - H2\n", "species:\r\n  - H2\r\n")
        metadata["source_sha256_utf8"] = collect(codeunits(bytes2hex(SHA.sha256(source))))
        for checkout in ("species:\n  - H2\n", "species:\r\n  - H2\r\n")
            write(mechanism, checkout)
            @test Arrhenius._validate_sidecar_metadata(metadata, mechanism) === nothing
        end
        write(mechanism, "species:\n  - O2\n")
        @test_throws ArgumentError Arrhenius._validate_sidecar_metadata(metadata, mechanism)
    end
end

@testset "single-region NASA7" begin
    yaml = Dict(
        "phases" => [Dict(
            "species" => ["A", "B"],
            "elements" => ["X"],
        )],
        "reactions" => Any[],
        "species" => [
            Dict(
                "name" => "A",
                "composition" => Dict("X" => 1),
                "thermo" => Dict(
                    "temperature-ranges" => [200.0, 3000.0],
                    "data" => [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]],
                ),
            ),
            Dict(
                "name" => "B",
                "composition" => Dict("X" => 1),
                "thermo" => Dict(
                    "temperature-ranges" => [200.0, 800.0, 3000.0],
                    "data" => [
                        [8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0],
                        [15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0],
                    ],
                ),
            ),
        ],
    )
    thermo = Arrhenius.IdealGasThermo(yaml)
    @test thermo.nasa_low[1, :] == thermo.nasa_high[1, :]
    @test thermo.Trange[1, :] == [200.0, 1000.0, 3000.0]
    @test thermo.Trange[2, :] == [200.0, 800.0, 3000.0]
end

@testset "PLOG interpolation and duplicate rates" begin
    plog = Arrhenius.PlogData(
        [7],
        [2],
        [1, 3],
        [1.0e5, 1.0e6],
        [1, 3, 4],
        [2.0 0.0 0.0; 3.0 0.0 0.0; 20.0 0.0 0.0],
    )
    T = 1000.0
    @test Arrhenius._plog_rate(plog, 1, T, 1.0e4, log(T)) ≈ 5.0
    @test Arrhenius._plog_rate(plog, 1, T, sqrt(1.0e11), log(T)) ≈ 10.0
    @test Arrhenius._plog_rate(plog, 1, T, 1.0e7, log(T)) ≈ 20.0
    concentrations = [3.0, 0.25]
    @test Arrhenius._plog_rate_with_collider(
        plog, 1, T, 1.0e4, log(T), concentrations,
    ) ≈ 1.25
    legacy_plog = Arrhenius.PlogData(
        [7], [1, 3], [1.0e5, 1.0e6], [1, 3, 4],
        [2.0 0.0 0.0; 3.0 0.0 0.0; 20.0 0.0 0.0],
    )
    @test legacy_plog.collider_indices == [0]
    plog32 = Arrhenius.PlogData(
        [7], [2], [1, 3], Float32[1.0e5, 1.0e6], [1, 3, 4],
        Float32[2.0 0.0 0.0; 3.0 0.0 0.0; 20.0 0.0 0.0],
    )
    T32 = 1000.0f0
    @test eltype(plog32.Arrhenius_coeffs) == Float32
    @test Arrhenius._plog_rate(plog32, 1, T32, 1.0f4, log(T32)) ≈ 5.0f0
    log_rates = LogRateData(
        Float32[], Float32[], Float32[],
        Float32[], Float32[], Float32[],
        Float32.(log.([2.0, 3.0, 20.0])),
        zeros(Float32, 3),
        zeros(Float32, 3),
    )
    @test Arrhenius._plog_rate(
        log_rates,
        plog,
        1,
        T,
        sqrt(1.0e11),
        log(T32),
        Float32(4184.0 / R) / T32,
    ) ≈ 10.0 rtol=1.0e-5
end

@testset "jl" begin
    # Write your tests here.

    gas = CreateSolution(joinpath(@__DIR__,"..","mechanism","gri30.yaml"))
    ns = gas.n_species

    Y0 = ones(ns) ./ ns
    T0 = 1200.0
    P = one_atm
    wdot = set_states(gas, T0, P, Y0)
    @show size(wdot)
    @show wdot[1:3]

    mean_MW = 1.0 / dot(Y0, 1 ./ gas.MW)
    density = P / R / T0 * mean_MW
    X = Y2X(gas, Y0, mean_MW)
    C = Y2C(gas, Y0, density)
    h_mole = get_H(gas, T0, Y0, X)
    S0 = get_S(gas, T0, P, X)
    thermo_buffer = zeros(ns)
    @test cal_h_RT!(thermo_buffer, gas, T0, P, X) ≈ cal_h_RT(gas, T0, P, X)
    @test cal_s0_R!(thermo_buffer, gas, T0, P, X) ≈ cal_s0_R(gas, T0, P, X)
    @test cal_cp_R!(thermo_buffer, gas, T0, P, X) ≈ cal_cp_R(gas, T0, P, X)
    kinetics_workspace = KineticsWorkspace(gas.reaction)
    wdot_buffer = zeros(ns)
    @test wdot!(
        wdot_buffer, gas.reaction, T0, C, S0, h_mole, kinetics_workspace,
    ) ≈ wdot
    @test_throws ArgumentError convert_precision(gas, Float32)
    baseline_qdot = copy(wdot!(
        wdot_buffer,
        gas.reaction,
        T0,
        C,
        S0,
        h_mole,
        kinetics_workspace;
        get_qdot=true,
    ))
    log_rate_data = LogRateData(gas.reaction, Float32)
    mixed_qdot = copy(wdot!(
        wdot_buffer,
        gas.reaction,
        T0,
        C,
        S0,
        h_mole,
        kinetics_workspace;
        get_qdot=true,
        log_rate_data=log_rate_data,
    ))
    mixed_wdot = gas.reaction.vk * mixed_qdot
    @test all(isfinite, mixed_wdot)
    @test norm(mixed_wdot - wdot, Inf) / norm(wdot, Inf) < 1.0e-3
    rate_multipliers = ones(gas.n_reactions)
    rate_multipliers[1] = 2.0
    perturbed_qdot = wdot!(
        wdot_buffer,
        gas.reaction,
        T0,
        C,
        S0,
        h_mole,
        kinetics_workspace;
        get_qdot=true,
        rate_multipliers=rate_multipliers,
    )
    @test perturbed_qdot[1] ≈ 2.0 * baseline_qdot[1]
    @test_throws DimensionMismatch wdot!(
        wdot_buffer,
        gas.reaction,
        T0,
        C,
        S0,
        h_mole,
        kinetics_workspace;
        rate_multipliers=ones(gas.n_reactions - 1),
    )
    function first_rate_of_progress(multiplier)
        multipliers = ones(typeof(multiplier), gas.n_reactions)
        multipliers[1] = multiplier
        return wdot_func(
            gas.reaction,
            T0,
            C,
            S0,
            h_mole;
            get_qdot=true,
            rate_multipliers=multipliers,
        )[1]
    end
    @test ForwardDiff.derivative(first_rate_of_progress, 1.0) ≈ baseline_qdot[1]

    X_buffer = similar(Y0)
    Y2X!(X_buffer, gas, Y0)
    @test @allocated(Y2X!(X_buffer, gas, Y0)) == 0

    u0 = vcat(Y0, T0)
    function f(u)
        T = u[end]
        Y = @view(u[1:ns])
        mean_MW = 1.0 / dot(Y, 1 ./ gas.MW)
        ρ_mass = P / R / T * mean_MW
        X = Y2X(gas, Y, mean_MW)
        C = Y2C(gas, Y, ρ_mass)
        cp_mole, cp_mass = get_cp(gas, T, X, mean_MW)
        h_mole = get_H(gas, T, Y, X)
        S0 = get_S(gas, T, P, X)
        wdot = wdot_func(gas.reaction, T, C, S0, h_mole)
        Tdot = -dot(h_mole, wdot) / ρ_mass / cp_mass
        return Tdot
    end
    grad = ForwardDiff.gradient(f, u0)
    @show size(grad)
    @show grad[1:3]

end
