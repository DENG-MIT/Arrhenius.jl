using Arrhenius
using Test

@testset "Arrhenius.jl" begin
    # Write your tests here.
    TODO: test not implemented yet
    gas = Arrhenius.CreateSolution("python/h2o2.yaml")
    mgas = Arrhenius.CreateMSolution(gas)

    Y = ones(gas.n_species) ./ gas.n_species
    Arrhenius.set_states(gas, mgas, 1200., Arrhenius.one_atm, Y)

    @show mgas.ρ_mass
    @show mgas.wdot
end
