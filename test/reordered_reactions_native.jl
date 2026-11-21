using Arrhenius, LinearAlgebra, Test

# Fixed-step explicit RK4 -- deliberately not this package's own adaptive
# solver infrastructure, so these tests carry no solver dependency (matching
# this test suite's existing convention of exercising `reactor_rhs`/
# `reactor_jacobian!` directly rather than integrating through an external
# ODE package). The point is not realistic ignition-delay timing; it is
# exercising the reordered-kinetics opt-in across many repeated calls along
# a chemically-active trajectory, not just at a single instant.
function _rk4_trajectory!(rhs, u, dt, nsteps)
    du1, du2, du3, du4, tmp = similar(u), similar(u), similar(u), similar(u), similar(u)
    for _ in 1:nsteps
        rhs(du1, u, nothing, 0.0)
        tmp .= u .+ (dt / 2) .* du1
        rhs(du2, tmp, nothing, 0.0)
        tmp .= u .+ (dt / 2) .* du2
        rhs(du3, tmp, nothing, 0.0)
        tmp .= u .+ dt .* du3
        rhs(du4, tmp, nothing, 0.0)
        u .+= (dt / 6) .* (du1 .+ 2du2 .+ 2du3 .+ du4)
    end
    return u
end

@testset "ReorderedReactions ignition trajectories match the generic path" begin
    cases = (
        (mechanism="gri30", path=joinpath(@__DIR__, "..", "mechanism", "gri30.yaml"),
            T=1800.0, P=10one_atm, X=Dict("CH4" => 1.0, "O2" => 2.0, "N2" => 7.52)),
        (mechanism="stagni_2023 (PLOG)", path=joinpath(@__DIR__, "..", "mechanism", "stagni_2023.yaml"),
            T=1600.0, P=5one_atm, X=Dict("NH3" => 1.0, "O2" => 0.75, "N2" => 2.82)),
    )
    for case in cases
        gas = CreateSolution(case.path)
        reactor = IdealGasReactor(gas; temperature=case.T, pressure=case.P, mole_fractions=case.X)
        u_generic = reactor_state(reactor)
        u_reordered = copy(u_generic)
        rhs_generic = reactor_rhs(reactor)
        rhs_reordered = reactor_rhs(reactor; use_reordered_kinetics=true)
        dt, nsteps = 1e-10, 5000
        _rk4_trajectory!(rhs_generic, u_generic, dt, nsteps)
        _rk4_trajectory!(rhs_reordered, u_reordered, dt, nsteps)

        @test all(isfinite, u_generic) && all(isfinite, u_reordered)
        @test u_generic[end] != case.T  # sanity: chemistry is actually active over this trajectory
        @test u_reordered ≈ u_generic rtol=1e-6
        @test u_reordered[end] ≈ u_generic[end] rtol=1e-8
        @test sum(u_reordered[1:end-1]) ≈ 1.0 rtol=1e-8  # mass fractions still sum to 1
    end
end

@testset "ReorderedReactions kernel matches the generic path at flame-representative states" begin
    # Pull real T/Y states from a converged premixed flame solve (test/flames_native.jl's
    # own construction pattern) rather than hand-picked numbers -- Flames.jl itself always
    # calls wdot! against gas.reaction (plain Reactions; it isn't wired to opt into
    # ReorderedReactions), so this checks that the kernel a future flame integration would
    # engage already agrees with the one Flames.jl uses today, at states a flame solve
    # actually produces (multiple reaction regimes simultaneously active: cold inflow,
    # reaction zone, hot burned gas).
    gas = CreateSolution(joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml"))
    X = mole_fractions(gas, "H2:1.1,O2:1,AR:5")
    f = FreeFlame(gas; X, discretization=:finite_difference)
    solve!(f)
    @test f.converged

    rr = ReorderedReactions(gas.reaction)
    n = gas.n_species
    T = temperature(f)
    Y = mass_fractions(f)
    sample_points = unique((1, argmax(T), size(Y, 2)))  # cold inflow, hottest point, burned-gas end

    for j in sample_points
        Tj = T[j]
        Xj = Y[:, j] ./ gas.MW
        Xj ./= sum(Xj)
        Cj = Xj .* (one_atm / (R * Tj))
        S0 = cal_s0_R(gas, Tj, one_atm, Xj) .* R
        h = cal_h_RT(gas, Tj, one_atm, Xj) .* (R * Tj)

        wdot_generic = zeros(n)
        wdot!(wdot_generic, gas.reaction, Tj, Cj, S0, h, KineticsWorkspace(gas.reaction); pressure=one_atm)
        wdot_reordered = zeros(n)
        wdot!(wdot_reordered, rr, Tj, Cj, S0, h, KineticsWorkspace(rr); pressure=one_atm)

        @test wdot_reordered ≈ wdot_generic rtol=1e-8 atol=1e-6 * maximum(abs, wdot_generic)
    end
end
