using Arrhenius, SHA, Test

@testset "Cantera FreeFlame flux gradient basis" begin
    mechanism = joinpath(@__DIR__, "..", "data", "cantera_flux_gradient_mass", "h2o2.yaml")
    @test bytes2hex(sha256(read(mechanism))) ==
          "0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4"
    @test bytes2hex(sha256(read(mechanism * ".npz"))) ==
          "1cc087ea45a8451ebf575ee35d9edee90445c6d08260f0a92669ee4457f528b5"

    gas = CreateSolution(mechanism)
    flame = FreeFlame(gas; T=300.0, P=one_atm, X="H2:1.1, O2:1, AR:5.3",
        width=0.05, transport_model=:mixture_averaged)
    refinement = (; ratio=3.0, slope=0.08, curve=0.12, prune=0.0)

    @test flame.flux_gradient_basis == :mole
    solve!(flame; auto=true, refinement...)
    @test flame.converged
    molar_speed = flame_speed(flame)
    molar_grid, molar_state = copy(flame.grid), copy(flame.state)

    set_transport!(flame, :mixture_averaged; flux_gradient_basis=:mass)
    @test flame.flux_gradient_basis == :mass
    @test isequal(flame.grid, molar_grid)
    @test isequal(flame.state, molar_state)
    solve!(flame; auto=false, refinement...)
    @test flame.converged
    mass_speed = flame_speed(flame)

    # Cantera 4.0.0a2 (726522be4e2a13454d8415b7ef799d621f665cf3),
    # TestFreeFlame.test_flux_gradient_mass on the same YAML, in m/s.
    for (name, molar, mass) in (("Cantera reference", 0.6070213069296682, 0.6175342760444065),
                               ("Arrhenius", molar_speed, mass_speed))
        @testset "$name" begin
            @test isfinite(molar) && molar > 0
            @test isfinite(mass) && mass > 0
            @test molar != mass
            # pytest.approx(sl_mass, rel=0.1) uses the mass speed as its scale.
            @test abs(molar - mass) <= max(0.1 * abs(mass), 1e-12)
        end
    end
end
