using Arrhenius, Test, NPZ

@testset "Cantera fixed-temperature burner first steady velocity" begin
    reference = NPZ.npzread(joinpath(@__DIR__, "data", "cantera_h2o2_fixed_t_velocity.npz"))
    gas = CreateSolution(joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml"))
    grid = 0.03 .* [0.0, 0.25, 0.5, 0.75, 1.0]
    burner = BurnerFlame(gas; mdot=0.07, T=300.0, P=one_atm,
        X="H2:1.1,O2:1,AR:5", grid, discretization=:finite_difference)
    set_temperature_profile!(burner, [0.0, 0.25, 1.0], [300.0, 1600.0, 300.0])
    Arrhenius._fresh_analytic_guess!(burner, grid; cantera_compatible=true)
    initial_T_rel = maximum(abs.(temperature(burner) ./ reference["initial_T"] .- 1))
    initial_u_rel = maximum(abs.(velocity(burner) ./ reference["initial_velocity"] .- 1))
    @test burner.grid == reference["initial_grid"]
    @test temperature(burner) ≈ reference["initial_T"] rtol=3e-11 atol=1e-10
    @test velocity(burner) ≈ reference["initial_velocity"] rtol=5e-11 atol=1e-10
    @info "Cantera fixed-T startup relative errors" T=initial_T_rel velocity=initial_u_rel
    solve!(burner; cantera_compatible=true, refine_grid=false)
    @test burner.converged
    @test burner.grid == reference["grid"]
    @test temperature(burner) ≈ reference["T"] rtol=2e-5 atol=1e-10
    first_steady_u_rel = maximum(abs.(velocity(burner) ./ reference["velocity"] .- 1))
    @info "Cantera fixed-T first-steady velocity relative error" velocity=first_steady_u_rel
    @test first_steady_u_rel ≤ 2.2e-6
end
