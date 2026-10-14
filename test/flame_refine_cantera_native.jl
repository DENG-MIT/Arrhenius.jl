@testset "Cantera flame refinement same state" begin
    gas = CreateSolution(joinpath(pkgdir(Arrhenius), "mechanism", "gri30.yaml"))
    flame = BurnerFlame(gas; mdot=0.09, T=300.0, P=one_atm,
        X="CH4:0.4, O2:1, N2:3.76", width=0.02, discretization=:finite_difference)
    data = Arrhenius.npzread(joinpath(@__DIR__, "cantera_refine_stage5.npz"))
    @test flame.grid ≈ 0.02 .* [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0] atol=1e-15
    grid = Vector{Float64}(data["grid"])
    refined_grid = Vector{Float64}(data["refined_grid"])
    temperature_profile = Vector{Float64}(data["temperature"])
    velocity_profile = Vector{Float64}(data["velocity"])
    Y = Matrix{Float64}(data["mass_fractions"])
    mean_molecular_weight = vec(gas.MW' * Y)
    mass_flux = flame.pressure .* velocity_profile ./
        (R .* temperature_profile .* mean_molecular_weight)
    state = zeros(gas.n_species + 2, length(grid))
    state[1, :] .= temperature_profile ./ 1000
    state[2:end-1, :] .= Y
    state[end, :] .= mass_flux
    flame.grid = grid
    flame.state = state
    flame.anchor = 1
    @test Arrhenius._refine_flame!(flame; ratio=3.0, slope=0.25,
        curve=0.35, prune=0.05)
    @test flame.grid == refined_grid
end
