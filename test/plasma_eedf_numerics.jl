module PlasmaEEDFNumericsTests
using Test, LinearAlgebra, Arrhenius
using Arrhenius: ElectronCollision, EffectiveCollision, ElasticCollision,
    ExcitationCollision, IonizationCollision, AttachmentCollision, EEDFWorkspace,
    _eedf_simpson, _integral_pq, _sg_coefficients, _add_collision_pq!, _eedf_norm

@testset "EEDF numerical identities" begin
    @testset "nonuniform Simpson quadrature" begin
        x = [0.0, 0.2, 0.9, 1.7, 2.0]
        @test _eedf_simpson(x .^ 2, x) ≈ 8/3 atol=2e-15
        @test _eedf_simpson([2.0, 6.0], [1.0, 3.0]) == 8.0
        @test_throws ArgumentError _eedf_simpson([1.0, 2.0], [1.0, 1.0])
    end

    @testset "piecewise-linear cross-section integral" begin
        a, b = 0.37, 1.91
        u0, u1 = 2.1e-20, 7.4e-20
        g, x0 = -0.63, 1.25
        analytic = _integral_pq(a, b, u0, u1, g, x0)
        # Independent high-order Gauss-Legendre quadrature (8 point).
        nodes = [-0.9602898564975363, -0.7966664774136267,
                 -0.5255324099163290, -0.1834346424956498,
                  0.1834346424956498,  0.5255324099163290,
                  0.7966664774136267,  0.9602898564975363]
        weights = [0.1012285362903763, 0.2223810344533745,
                   0.3137066458778873, 0.3626837833783620,
                   0.3626837833783620, 0.3137066458778873,
                   0.2223810344533745, 0.1012285362903763]
        midpoint, halfwidth = (a+b)/2, (b-a)/2
        quadrature = 0.0
        for (node, weight) in zip(nodes, weights)
            energy = midpoint + halfwidth*node
            sigma = u0 + (u1-u0)*(energy-a)/(b-a)
            quadrature += weight * energy * sigma * exp((x0-energy)*g)
        end
        quadrature *= halfwidth
        @test analytic ≈ quadrature rtol=2e-13
        @test _integral_pq(a, b, u0, u1, 0.0, x0) ≈
              halfwidth * sum(weight * (midpoint + halfwidth*node) *
              (u0 + (u1-u0)*(midpoint + halfwidth*node-a)/(b-a))
              for (node, weight) in zip(nodes, weights)) rtol=2e-14
    end

    @testset "Scharfetter-Gummel stable limits" begin
        D, h = 3.2, 0.7
        a0, a1 = _sg_coefficients(0.0, D, h)
        @test a0 ≈ D/h rtol=2e-15
        @test a1 ≈ -D/h rtol=2e-15
        a0p, a1p = _sg_coefficients(1e4, D, h)
        @test isfinite(a0p) && isfinite(a1p)
        @test a0p ≈ 1e4 rtol=2e-15
        @test abs(a1p) < 1e-100
        a0n, a1n = _sg_coefficients(-1e4, D, h)
        @test abs(a0n) < 1e-100
        @test a1n ≈ -1e4 rtol=2e-15
    end
end

function _toy_collision(kind; threshold=0.0, sigma=[2e-20, 2e-20])
    ElectronCollision("N2", kind, threshold, [0.0, 4.0], sigma, :test)
end

@testset "EEDF collision bookkeeping" begin
    edges = [0.0, 1.0, 2.0, 3.0, 4.0]
    state_for(model; EN=2e-19) = EEDFState(model; T=300.0, P=101325.0,
        mole_fractions=Dict("N2" => 1.0), molecular_weights=Dict("N2" => 28.014),
        reduced_field=EN)

    excitation = _toy_collision(ExcitationCollision; threshold=1.0)
    ionization = _toy_collision(IonizationCollision; threshold=1.0)
    attachment = _toy_collision(AttachmentCollision; threshold=0.0)
    effective = _toy_collision(EffectiveCollision)

    for (collision, incoming) in ((excitation, 1.0), (ionization, 2.0), (attachment, 0.0))
        model = EEDFModel(edges, [effective, collision], ["N2"])
        ws = EEDFWorkspace(model, state_for(model))
        cache = ws.collision_cache[2]
        @test cache.incoming_factor == incoming
        M = zeros(4, 4)
        g = zeros(4)
        _add_collision_pq!(M, cache, g, ws.centers)
        if collision.kind == AttachmentCollision
            @test all(M - Diagonal(diag(M)) .== 0)
            @test sum(M) < 0
        elseif collision.kind == IonizationCollision
            @test sum(M) > 0
        else
            @test sum(M) ≈ 0 atol=1e-27
        end
    end
end

@testset "EEDF validation and low-field solution" begin
    elastic = _toy_collision(EffectiveCollision)
    model = EEDFModel([0.0, 1.0, 2.0, 3.0], [elastic], ["N2"])
    original_edges = copy(model.energy_edges)
    original_sigma = copy(model.collisions[1].cross_section)
    state = EEDFState(model; T=300.0, P=101325.0,
        mole_fractions=Dict("N2" => 0.7, "unused" => 0.3),
        molecular_weights=Dict("N2" => 28.014), reduced_field=0.0)
    @test state.mole_fractions["N2"] == 1.0
    result = solve_eedf(model, state)
    @test result.converged
    @test result.iterations == 0
    @test all(isfinite, result.center_eedf) && all(result.center_eedf .> 0)
    @test all(isfinite, result.edge_eedf) && all(result.edge_eedf .> 0)
    @test isfinite(result.mobility) && result.mobility > 0
    @test _eedf_norm(result.center_eedf, result.centers) ≈ 1.0 atol=1e-12
    @test model.energy_edges == original_edges
    @test model.collisions[1].cross_section == original_sigma

    @test_throws ArgumentError EEDFState(model; T=300.0, P=101325.0,
        mole_fractions=Dict("N2" => 0.0), molecular_weights=Dict("N2" => 28.014),
        reduced_field=2e-19)
    @test_throws ArgumentError EEDFState(model; T=300.0, P=101325.0,
        mole_fractions=Dict{String,Float64}(), molecular_weights=Dict("N2" => 28.014),
        reduced_field=2e-19)
    bad_grid = EEDFModel([0.0, 1.0], [elastic], ["N2"])
    @test_throws ArgumentError EEDFWorkspace(bad_grid,
        EEDFState(bad_grid; T=300.0, P=101325.0,
            mole_fractions=Dict("N2" => 1.0), molecular_weights=Dict("N2" => 28.014),
            reduced_field=2e-19))
    duplicate = EEDFModel([0.0, 1.0, 2.0], [elastic, _toy_collision(ElasticCollision)], ["N2"])
    duplicate_state = EEDFState(duplicate; T=300.0, P=101325.0,
            mole_fractions=Dict("N2" => 1.0), molecular_weights=Dict("N2" => 28.014),
            reduced_field=2e-19)
    @test_throws ArgumentError EEDFWorkspace(duplicate, duplicate_state)
end

@testset "EEDF density and continuation contracts" begin
    model = EEDFModel(collect(0.0:40.0), [_toy_collision(EffectiveCollision)], ["N2"])
    common = (T=300.0, P=101325.0, mole_fractions=Dict("N2"=>1.0),
        molecular_weights=Dict("N2"=>28.014))
    low = EEDFState(model; common..., reduced_field=0.0)
    density = low.number_density
    @test density == low.P/(Arrhenius._EEDF_BOLTZMANN*low.T)
    old_positional = EEDFState(low.T, low.P, low.mole_fractions, low.molecular_weights,
        low.reduced_field, low.frequency)
    @test old_positional.number_density == density
    cold = solve_eedf(model, low)
    # The original pulse grid at300K underflows in the Maxwellian tail.
    @test any(iszero, cold.center_eedf)
    @test all(cold.center_eedf .>= 0)
    @test Arrhenius._eedf_norm(cold.center_eedf, cold.centers) ≈ 1.0 atol=1e-12
    same = solve_eedf(model, EEDFState(model; common..., reduced_field=0.0, number_density=density))
    doubled = solve_eedf(model, EEDFState(model; common..., reduced_field=0.0, number_density=2density))
    @test same.center_eedf == cold.center_eedf
    @test same.edge_eedf == cold.edge_eedf
    @test same.mobility == cold.mobility
    @test doubled.center_eedf == cold.center_eedf
    @test doubled.mobility == cold.mobility/2
    for invalid in (0.0, -1.0, Inf, NaN)
        @test_throws ArgumentError EEDFState(model; common..., reduced_field=0.0, number_density=invalid)
    end
    malformed = EEDFState(low.T, low.P, low.mole_fractions, low.molecular_weights, 0.0, 0.0, NaN)
    @test_throws ArgumentError solve_eedf(model, malformed)

    # Holding frequency/N fixed must preserve the AC operator, including its field term.
    ac = EEDFState(model; common..., reduced_field=2e-19, frequency=1e9, number_density=density)
    ac2 = EEDFState(model; common..., reduced_field=2e-19, frequency=2e9, number_density=2density)
    ws = EEDFWorkspace(model, ac); ws2 = EEDFWorkspace(model, ac2)
    Arrhenius._assemble_operator!(ws, model, ac, cold.center_eedf, 1e-300)
    Arrhenius._assemble_operator!(ws2, model, ac2, cold.center_eedf, 1e-300)
    @test ws.operator == ws2.operator

    # A deliberately different, normalized prior verifies low-field reset semantics.
    prior = deepcopy(cold)
    prior.center_eedf .= Arrhenius._maxwellian(prior.centers, 2.0)
    saved = deepcopy(prior)
    reset = solve_eedf(model, low; initial=prior)
    @test reset.center_eedf == cold.center_eedf
    @test reset.edge_eedf == cold.edge_eedf
    @test reset.mobility == cold.mobility
    @test prior.center_eedf == saved.center_eedf
    @test reset.center_eedf !== prior.center_eedf
    for change in (r->(r.center_eedf[1]=-1.0), r->(r.center_eedf[1]=NaN),
                   r->(r.center_eedf .*= 2), r->pop!(r.center_eedf),
                   r->(r.edges[end]+=1), r->(r.centers[end]+=1))
        bad = deepcopy(prior); change(bad)
        @test_throws ArgumentError solve_eedf(model, low; initial=bad)
    end
    # Finite edge samples cannot substitute for a compatible center distribution.
    hot = EEDFState(model; common..., reduced_field=2e-19)
    continued = solve_eedf(model, hot; initial=prior)
    @test continued.converged
    @test prior.center_eedf == saved.center_eedf
    @test continued.center_eedf !== prior.center_eedf
end
end # module
