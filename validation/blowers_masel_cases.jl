using Arrhenius, Test, NPZ, LinearAlgebra
root = ARGS[1]
gas = CreateSolution(joinpath(root,"blowers-masel.yaml"))
reference = npzread(joinpath(root,"reference.npz"))
@testset "Cantera 4 Blowers–Masel rates" begin
    @test gas.reaction.blowers_masel.reaction_indices == [2,3]
    for (j,T) in enumerate(reference["T"])
        rates = reaction_rate_constants(gas;T,X="H2:1,O2:1,CH4:1,AR:1")
        @test rates.forward ≈ reference["forward"][:,j] rtol=2e-12
        @test rates.reverse ≈ reference["reverse"][:,j] rtol=2e-12
    end
    rate = BlowersMaselRate(3.87e1,2.7,6260*1000*4.184,1e9)
    @test activation_energy.(Ref(rate),reference["deltaH"]) ≈ reference["barriers"] rtol=2e-12 atol=1e-6
    # Changing thermochemistry at fixed T must change both BM rates and Kc.
    T = 1000.
    x = mole_fractions(gas,"H2:1,O2:1,CH4:1,AR:1")
    c = one_atm/(R*T).*x
    h = cal_h_RT(gas,T,one_atm,x).*(R*T)
    s = cal_s0_R(gas,T,one_atm,x).*R
    w = KineticsWorkspace(gas.reaction)
    cache = Arrhenius._KineticsTemperatureCache(gas.reaction)
    for shift in (0.,1e7,-1e7)
        changed = copy(h)
        changed[species_index(gas,"H")] += shift
        rates = wdot!(similar(c),gas.reaction,T,c,s,changed,w;temperature_cache=cache,get_rate_constants=true)
        dh = dot(gas.reaction.vk[:,2],changed)
        @test rates.forward[2] ≈ rate_constant(rate,T,dh) rtol=2e-12
    end
end
