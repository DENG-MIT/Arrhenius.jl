using Arrhenius, NPZ, Test
gas = CreateSolution(ARGS[1])
reference = ARGS[2]
w = TransportWorkspace(gas)
n = gas.n_species
@testset "Cantera mixture transport" begin
    for row in eachrow(npzread(joinpath(reference, "transport.npz"))["states"])
        T, P = row[1:2]
        X = row[3:2+n]
        viscosity, conductivity = mixture_transport!(w, gas, P, T, X)
        @test viscosity ≈ row[3+n] rtol=2e-12
        @test conductivity ≈ row[4+n] rtol=2e-12
        @test w.diffusion ≈ row[5+n:end] rtol=2e-11
    end
end
@testset "Cantera equilibrium" begin
    oracle = npzread(joinpath(reference,"equilibrium.npz"))
    for (i, mode) in enumerate((:TP,:TP,:TP,:HP,:HP,:HP))
        eq = equilibrate(gas; T=oracle["Tin"][i], P=Arrhenius.one_atm,
            X=Dict("H2"=>1.1,"O2"=>1.,"AR"=>5.), mode)
        @test eq.T ≈ oracle["T"][i] atol=1e-3
        @test eq.X ≈ oracle["X"][:,i] atol=1e-7 rtol=1e-6
    end
end
