using Arrhenius
using ForwardDiff
using LinearAlgebra
using Test

struct _NegativeTrialReactorTag end

function _negative_trial_ad_column(reactor, gas, u, column)
    D = ForwardDiff.Dual{_NegativeTrialReactorTag,Float64,1}
    workspace = Arrhenius.ReactorWorkspace(gas, D)
    scratch = ntuple(_ -> zeros(D, length(u)), 4)
    rhs = Arrhenius.ReactorRHS(reactor, workspace, scratch...)
    dual_u = [ForwardDiff.Dual{_NegativeTrialReactorTag}(
        u[k], k == column ? 1.0 : 0.0) for k in eachindex(u)]
    dual_du = similar(dual_u)
    rhs(dual_du, dual_u, nothing, 0.0)
    return [ForwardDiff.partials(value)[1] for value in dual_du]
end

@testset "closed native ideal-gas reactors" begin
    gas = CreateSolution(joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml"))
    X = Dict("H2" => 2.0, "O2" => 1.0, "AR" => 4.0)
    initial = (; temperature=1001.0, pressure=one_atm, mole_fractions=X)
    reactor = IdealGasReactor(gas; initial...)
    state = reactor_state(reactor)
    properties = reactor_properties(reactor)
    @test properties.mass_fraction_sum ≈ 1.0
    @test properties.density ≈ one_atm * dot(properties.mole_fractions, gas.MW) / (R * 1001.0)
    @test properties.enthalpy ≈ cal_hmass_mean(gas, 1001.0, one_atm, properties.mole_fractions)
    @test properties.internal_energy ≈ cal_umass_mean(gas, 1001.0, one_atm, properties.mole_fractions)
    @test properties.cp ≈ cal_cpmass_mean(gas, 1001.0, one_atm, properties.mole_fractions)
    @test properties.cv ≈ cal_cvmass_mean(gas, 1001.0, one_atm, properties.mole_fractions)
    state[1] = 0.0
    @test reactor_state(reactor)[1] > 0
    from_mass = IdealGasReactor(gas; temperature=1001.0, mass_fractions=reactor.mass_fractions)
    @test reactor_state(from_mass) ≈ reactor_state(reactor)

    @test_throws ArgumentError IdealGasReactor(gas; temperature=1001.0)
    @test_throws ArgumentError IdealGasReactor(gas; initial..., mass_fractions=ones(gas.n_species))
    @test_throws ArgumentError IdealGasReactor(gas; initial..., constraint=:unknown)
    @test_throws ArgumentError IdealGasReactor(gas; initial..., energy=:unknown)
    @test_throws ArgumentError IdealGasReactor(gas; temperature=-1, mole_fractions=X)
    @test_throws ArgumentError IdealGasReactor(gas; temperature=1001, mole_fractions=Dict("missing" => 1))
    @test_throws ArgumentError IdealGasReactor(gas; temperature=1001, mass_fractions=zeros(gas.n_species))
    @test_throws ArgumentError IdealGasReactor(gas; temperature=1001, mass_fractions=fill(-1, gas.n_species))
    @test_throws DimensionMismatch IdealGasReactor(gas; temperature=1001, mass_fractions=[1.0])
    @test_throws DimensionMismatch IdealGasReactor(gas; initial..., rate_multipliers=[1.0])
    @test_throws ArgumentError IdealGasReactor(gas; initial..., rate_multipliers=fill(-1.0, gas.n_reactions))

    # Positive radical fractions exercise every Jacobian column and the energy
    # cancellation at a reacting state, independently of initial ignition rates.
    Y = collect(range(1.0, 2.0; length=gas.n_species))
    Y ./= sum(Y)
    for constraint in (:constant_pressure, :constant_volume), energy in (:adiabatic, :isothermal)
        r = IdealGasReactor(gas; temperature=1400.0, pressure=3one_atm,
                            mass_fractions=Y, constraint, energy)
        u = reactor_state(r)
        rhs = reactor_rhs(r)
        du = similar(u)
        rhs(du, u, nothing, 0.0)
        prop = reactor_properties(r, u)
        molar_source = set_states(gas, prop.temperature, prop.pressure, Y)
        @test du[1:end-1] ≈ molar_source .* gas.MW ./ prop.density rtol=1e-12
        @test abs(sum(du[1:end-1])) < 1e-12 * norm(du[1:end-1], 1)
        elemental_rate = gas.ele_matrix * (du[1:end-1] ./ gas.MW)
        @test norm(elemental_rate, Inf) < 1e-12 * norm(du[1:end-1], 1)
        if energy === :adiabatic
            species_energy = constraint === :constant_pressure ?
                cal_hmass(gas, prop.temperature, prop.pressure, prop.mole_fractions) :
                cal_umass(gas, prop.temperature, prop.pressure, prop.mole_fractions)
            capacity = constraint === :constant_pressure ? prop.cp : prop.cv
            species_source = dot(species_energy, du[1:end-1])
            @test abs(capacity * du[end] + species_source) < 1e-12 * abs(species_source)
        else
            @test du[end] == 0
            changed_temperature = copy(u)
            changed_temperature[end] += 100
            other = similar(u)
            rhs(other, changed_temperature, nothing, 0.0)
            @test other == du
        end
        J = zeros(length(u), length(u))
        before = copy(u)
        reactor_jacobian!(J, u, rhs)
        @test u == before
        direction = collect(range(-0.1, 0.1; length=length(u)))
        direction[end] = 100.0
        plus, minus = similar(u), similar(u)
        h = 1e-5
        rhs(plus, u + h * direction, nothing, 0.0)
        rhs(minus, u - h * direction, nothing, 0.0)
        @test J * direction ≈ (plus - minus) / (2h) rtol=1e-6
        @test_throws DimensionMismatch reactor_jacobian!(zeros(2, 2), u, rhs)

        shifted = copy(u)
        shifted[end] *= 1.1
        current = reactor_properties(r, shifted)
        if constraint === :constant_pressure
            @test current.pressure == r.pressure
        else
            @test current.density == r.density
        end
    end

    inactive = IdealGasReactor(gas; initial..., rate_multipliers=zeros(gas.n_reactions))
    rhs = reactor_rhs(inactive)
    du = similar(reactor_state(inactive))
    rhs(du, reactor_state(inactive), nothing, 0)
    @test all(iszero, du)
    J0 = zeros(length(du), length(du))
    reactor_jacobian!(J0, reactor_state(reactor), reactor_rhs(reactor))
    @test all(isfinite, J0)
    @test norm(J0, Inf) > 0

    problem = reactor_problem(reactor, (0.0, 0.001))
    other = reactor_problem(reactor, (0.0, 0.001))
    @test problem.f.workspace !== other.f.workspace
    @test problem.u0 == reactor_state(reactor)
    @test problem.tspan == (0.0, 0.001)
    @test_throws ArgumentError reactor_problem(reactor, (0.001, 0.0))
    @test_throws ArgumentError reactor_problem(reactor, (0.0, Inf))
    fill!(du, 1)
    problem.tgrad(du, problem.u0, nothing, 0)
    @test all(iszero, du)
    received = Ref(false)
    function inspect_problem(p; marker)
        received[] = marker && p.u0 == reactor_state(reactor)
        return :caller_owned_result
    end
    @test solve_reactor(reactor, (0.0, 0.001); integrator=inspect_problem, marker=true) === :caller_owned_result
    @test received[]
end

@testset "negative trial-state reactor Jacobian columns" begin
    gas = CreateSolution(joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml"))
    mixture = Dict("H2" => 2.0, "O2" => 1.0, "AR" => 4.0)
    column = findfirst(==("H"), gas.species_names)
    @test !isnothing(column)
    column = something(column)
    base_step = cbrt(eps(Float64)) * 1.0e-6

    for constraint in (:constant_pressure, :constant_volume)
        reactor = IdealGasReactor(gas; temperature=1400.0, pressure=3one_atm,
                                  mole_fractions=mixture, constraint)
        initial = reactor_state(reactor)
        @test initial[column] == 0.0
        for (regime, negative_fraction) in ((:within_one_step, -0.5 * base_step),
                                             (:between_one_and_two_steps, -1.5 * base_step),
                                             (:farther_negative, -1.0e-7))
            u = copy(initial)
            u[column] = negative_fraction
            raw_step = cbrt(eps(eltype(u))) * max(abs(u[column]), 1.0e-6)
            step = (u[column] + raw_step) - u[column]
            if regime === :within_one_step
                @test -step < u[column] < 0
            elseif regime === :between_one_and_two_steps
                @test -2 * step < u[column] < -step
            else
                @test u[column] < -2 * step
            end

            before = copy(u)
            rhs = reactor_rhs(reactor)
            value = similar(u)
            rhs(value, u, nothing, 0.0)
            J = zeros(length(u), length(u))
            reactor_jacobian!(J, u, rhs)
            oracle = _negative_trial_ad_column(reactor, gas, u, column)
            @test u == before
            @test all(isfinite, value)
            @test all(isfinite, @view J[:, column])
            @test all(isfinite, oracle)

            column_scale = max(norm(oracle, Inf), 1.0)
            tolerance = 2.0e-3 .* abs.(oracle) .+ 5.0e-7 * column_scale
            @test all(abs.(@view(J[:, column]) .- oracle) .<= tolerance)
        end
    end
end
