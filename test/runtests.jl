using Arrhenius
using Test

@testset "jl" begin
    # Write your tests here.
    # TODO: test not implemented yet

    gas = CreateSolution("python/JP10skeletal.yaml")
    mgas = CreateMSolution(gas)

    set_states(gas, mgas, 1200., one_atm, ones(gas.n_species) ./ gas.n_species)

    @show mgas.ρ_mass
    @show mgas.wdot

end
