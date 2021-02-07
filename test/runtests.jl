using Arrhenius
using Test

@testset "Arrhenius.jl" begin
    # Write your tests here.
    # TODO: test not implemented yet

    gas = Arrhenius.CreateSolution("python/JP10skeletal.yaml")
    mgas = Arrhenius.CreateMSolution(gas)

    Arrhenius.set_states(gas, mgas, 1200., Arrhenius.one_atm, ones(gas.n_species) ./ gas.n_species)

    @show mgas.ρ_mass
    @show mgas.wdot

end
