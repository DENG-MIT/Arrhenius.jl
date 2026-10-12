using Arrhenius, NPZ, LinearAlgebra, Test
BLAS.set_num_threads(1)
mechanisms, references = ARGS[1:2]
modes = (:TP,:TV,:HP,:UV,:SP,:SV)

@testset "Cantera ideal-gas equilibrium reference" begin
    for name in ("h2o2","gri30")
        gas = CreateSolution(joinpath(mechanisms,name*".yaml"))
        data = npzread(joinpath(references,"equilibrium-"*name*".npz"))
        version = String(vec(UInt8.(data["cantera_version_utf8"])))
        println("Checking ",name," against Cantera ",version)
        for i in eachindex(data["T0"])
            mode = modes[Int(data["mode"][i])]
            @testset "$name/$i/$mode" begin
                T0, P0, X0 = data["T0"][i], data["P0"][i], data["X0"][:,i]
                result = equilibrate(gas; T=T0,P=P0,X=X0,mode)
                @test result.T ≈ data["T"][i] rtol=2e-7 atol=2e-5
                @test result.P ≈ data["P"][i] rtol=2e-7
                @test maximum(abs.(result.X-data["X"][:,i])) < 1e-7
                @test maximum(abs.(result.Y-data["Y"][:,i])) < 1e-7
                for (key, property) in (("h",cal_hmass_mean),("u",cal_umass_mean),("s",cal_smass_mean))
                    @test property(gas,result.T,result.P,result.X) ≈ data[key][i] rtol=2e-7 atol=2e-2
                end
                rho = result.P*dot(gas.MW,result.X)/(R*result.T)
                @test rho ≈ data["rho"][i] rtol=2e-7
                initial_elements = gas.ele_matrix*X0/dot(gas.MW,X0)
                final_elements = gas.ele_matrix*result.X/dot(gas.MW,result.X)
                @test final_elements ≈ initial_elements rtol=2e-9 atol=1e-14
                for e in eachindex(initial_elements)
                    if initial_elements[e] > 0
                        @test final_elements[e] ≈ initial_elements[e] rtol=2e-8 atol=0
                    else
                        @test final_elements[e] == 0
                    end
                end
                @test all(>=(0),result.X)
                if mode in (:TV,:UV,:SV)
                    @test rho ≈ P0*dot(gas.MW,X0)/(R*T0) rtol=1e-10
                end
                # TP/TV equilibria must minimize the appropriate free energy.
                property = mode == :TP ? cal_gmass_mean : cal_amass_mean
                if mode in (:TP,:TV)
                    @test property(gas,result.T,result.P,result.X) <= property(gas,T0,P0,X0)+1e-5
                end
            end
        end
    end
end

@testset "Cantera equivalence-ratio example" begin
    gas = CreateSolution(joinpath(mechanisms,"gri30.yaml"))
    data = npzread(joinpath(references,"composition-gri30.npz"))["X"]
    cases = (
        (1.,"CH4","O2:.21,N2:.79",(;)),
        (1.,"CH4","O2:.233,N2:.767",(;basis=:mass)),
        (2.5,"CH4:1,O2:.01,CO:.05,N2:.1","O2:.2,N2:.8,CO2:.05,CH4:.01",(;)),
        (2.,"H2","O2",(;diluent="H2O",fraction=(diluent=.3,))),
        (2.,"H2","O2",(;diluent="CO2:.5,H2O:.5",fraction=(fuel=.1,),basis=:mass)),
        (2.,"H2:.5,H2O:.5","O2:.21,N2:.79",(;)),
    )
    for (i,(phi,fuel,oxidizer,kwargs)) in enumerate(cases)
        X = set_equivalence_ratio(gas,phi; fuel,oxidizer,kwargs...)
        @test X ≈ data[:,i] rtol=1e-12 atol=1e-14
        if i in (4,5)
            @test equivalence_ratio(gas,X;fuel,oxidizer,include_species=["H2","O2"]) ≈ phi
        else
            basis = get(kwargs,:basis,:mole)
            @test equivalence_ratio(gas,X;fuel,oxidizer,basis) ≈ phi
            eq = equilibrate(gas;T=300.,P=one_atm,X,mode=:HP)
            @test equivalence_ratio(gas,eq.X;fuel,oxidizer,basis) ≈ phi rtol=1e-8
        end
    end
    X = set_mixture_fraction(gas,.055;fuel="CH4",oxidizer="O2:.21,N2:.79")
    @test X ≈ data[:,7] rtol=1e-12 atol=1e-14
    @test mixture_fraction(gas,X;fuel="CH4",oxidizer="O2:.21,N2:.79") ≈ .055
    @test mixture_fraction(gas,X;fuel="CH4",oxidizer="O2:.21,N2:.79",element="C") ≈ .055
end
