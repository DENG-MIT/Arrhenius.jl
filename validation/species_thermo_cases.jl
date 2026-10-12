using Arrhenius, YAML, NPZ, LinearAlgebra, Test, ForwardDiff
BLAS.set_num_threads(1)
references = ARGS[1]
base = CreateSolution(joinpath(@__DIR__,"..","mechanism","gri30.yaml"))

@testset "Cantera species thermodynamic models" begin
    for name in ("air-nasa9","mixed-thermo","constant-cp-units","single-shomate","neutral-nasa9")
        yaml = YAML.load_file(joinpath(references,name*".yaml"))
        data = npzread(joinpath(references,name*".npz"))
        thermo = IdealGasThermo(yaml)
        names = String.(yaml["phases"][1]["species"])
        gas = Arrhenius.Solution(length(names),base.n_reactions,data["MW"],names,
            String.(yaml["phases"][1]["elements"]),data["elements"],thermo,base.trans,base.reaction)
        # The borrowed kinetic/transport records are never evaluated: these
        # fixtures exercise species and mixture thermodynamics only.
        X = fill(1/gas.n_species,gas.n_species)
        version = String(vec(UInt8.(data["version_utf8"])))
        println(name,": ",gas.n_species," species / ",length(data["T"])," states; Cantera ",version)
        for i in eachindex(data["T"])
            T, P = data["T"][i], data["P"][i]
            properties = species_thermo(thermo,T;P)
            @test properties.cp_R ≈ data["cp_R"][:,i] rtol=2e-12 atol=1e-12
            @test properties.h_RT ≈ data["h_RT"][:,i] rtol=2e-12 atol=1e-12
            @test properties.s_R ≈ data["s_R"][:,i] rtol=2e-12 atol=1e-12
            @test cal_cp_R(gas,T,P,X) ≈ properties.cp_R rtol=1e-14
            @test cal_h_RT(gas,T,P,X) ≈ properties.h_RT rtol=1e-14
            @test cal_s0_R(gas,T,P,X) .- log(P/one_atm) ≈ properties.s_R rtol=1e-14
            for (allocate,inplace) in ((cal_cp_R,cal_cp_R!),(cal_h_RT,cal_h_RT!),(cal_s0_R,cal_s0_R!))
                output = zeros(gas.n_species)
                inplace(output,gas,T,P,X)
                @test output == allocate(gas,T,P,X)
            end
        end
        for T in (700.,1700.,4000.)
            # Thermodynamic identities catch independently incorrect integrals.
            dhdT = ForwardDiff.derivative(t -> cal_hmass_mean(gas,t,one_atm,X),T)
            dsdT = ForwardDiff.derivative(t -> cal_smass_mean(gas,t,one_atm,X),T)
            @test dhdT ≈ cal_cpmass_mean(gas,T,one_atm,X) rtol=1e-11 atol=1e-10
            @test T*dsdT ≈ cal_cpmass_mean(gas,T,one_atm,X) rtol=1e-11 atol=1e-10
        end
        thermo32 = Arrhenius._convert_precision(thermo,Float32)
        @test !isnothing(thermo32.extra)
        @test thermo32.extra.models == thermo.extra.models
        @test thermo32.extra.models !== thermo.extra.models
        for T in (700.f0,1700.f0,4000.f0)
            @test species_thermo(thermo32,T;P=Float32(one_atm)).cp_R ≈ cal_cp_R(gas,T,one_atm,X) rtol=2e-5 atol=1e-6
        end
        if name == "neutral-nasa9"
            equilibrium = npzread(joinpath(references,"neutral-nasa9-equilibrium.npz"))
            for i in eachindex(equilibrium["T"])
                result = equilibrate(gas;T=equilibrium["T"][i],X="N2:.79,O2:.21")
                @test maximum(abs.(result.X-equilibrium["X"][:,i])) < 1e-8
            end
        end
    end
end
