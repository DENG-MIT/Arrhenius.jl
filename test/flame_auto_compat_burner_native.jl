using Arrhenius, Test

@testset "Cantera-compatible burner auto solve" begin
    mechanism = joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml")
    gas = CreateSolution(mechanism)
    X = mole_fractions(gas, "H2:1.1,O2:1,AR:5")
    grid = 0.03 .* [0.0, 0.2, 0.3, 0.6, 0.8, 1.0]
    burner = BurnerFlame(gas; mdot=0.06, T=300.0, P=one_atm, X, grid,
        discretization=:finite_difference)
    solve!(burner; cantera_compatible=true, refine_grid=false)
    @test burner.converged
    sum_y_error = maximum(abs.(sum(burner.state[2:end-1, :]; dims=1) .- 1))
    mass_flux_error = maximum(abs.(burner.state[end, :] ./ burner.mass_flux .- 1))
    @info "burner invariants" sum_y_error mass_flux_error
    @test sum_y_error < 1e-10
    @test mass_flux_error < 1e-12
end
