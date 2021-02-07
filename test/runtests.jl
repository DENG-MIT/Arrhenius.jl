using Arrhenius
using Arrhenius:one_atm,CreateSolution,CreateMSolution,set_states
using Test

@testset "jl" begin
    # Write your tests here.

    gas = CreateSolution("../python/JP10skeletal.yaml")
    mgas = CreateMSolution(gas)

    set_states(gas, mgas, 1200., one_atm, ones(gas.n_species) ./ gas.n_species)

    @show mgas.ρ_mass
    @show mgas.wdot

end
