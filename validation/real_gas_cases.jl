using Arrhenius, NPZ, Test, LinearAlgebra
BLAS.set_num_threads(1)
isdefined(Arrhenius,:RedlichKwongThermo) || Base.include(Arrhenius,joinpath(@__DIR__,"..","src","RealGasThermo.jl"))
references = ARGS[1]
@testset "Cantera 4 Redlich–Kwong properties" begin
    for name in ("dodecane","binary-temperature","critical-parameters","ideal-limit","co2-roots")
        model = RedlichKwongThermo(joinpath(references,name*".yaml"))
        data = npzread(joinpath(references,name*".npz"))
        work = RedlichKwongWorkspace(model)
        @test model.MW ≈ data["MW"] rtol=1e-14
        println(name,": ",length(data["T"])," states, ",length(model.MW)," species")
        for i in eachindex(data["T"])
            T,P,X = data["T"][i],data["P"][i],data["X"][:,i]
            root = data["root"][i] == 0 ? :gas : :liquid
            result = redlich_kwong_state(model;T,P,X,root)
            @test result.rho ≈ data["density"][i] rtol=3e-10
            for (ours,theirs) in ((:h,"enthalpy_mole"),(:u,"int_energy_mole"),(:s,"entropy_mole"),
                (:cp,"cp_mole"),(:cv,"cv_mole"),(:compressibility,"isothermal_compressibility"),
                (:expansion,"thermal_expansion_coeff"),(:sound_speed,"sound_speed"))
                @test getproperty(result,ours) ≈ data[theirs][i] rtol=3e-10 atol=1e-8
            end
            @test result.u_TV ≈ data["partial_molar_int_energies_TV"][:,i] rtol=3e-10 atol=1e-6
            @test result.partial_molar_enthalpies ≈ data["partial_molar_enthalpies"][:,i] rtol=3e-10 atol=1e-6
            @test result.partial_molar_volumes ≈ data["partial_molar_volumes"][:,i] rtol=3e-10 atol=1e-11
            ref = species_thermo(model.reference,T;P)
            chemical = R*T.*(ref.h_RT-ref.s_R+log.(max.(X,1e-300))+result.lnphi)
            @test chemical ≈ data["chemical_potentials"][:,i] rtol=3e-10 atol=2e-5
            fixed = redlich_kwong_properties!(work,model,T,data["density"][i],X)
            @test fixed.P ≈ P rtol=3e-10
            @test fixed.u ≈ result.u rtol=3e-10 atol=1e-7
        end
    end
end
