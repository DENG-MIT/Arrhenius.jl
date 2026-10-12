using Arrhenius, NPZ, LinearAlgebra, Test
BLAS.set_num_threads(1)
for (name,file) in ((:RedlichKwongThermo,"RealGasThermo.jl"),(:RedlichKwongReactor,"RealGasReactors.jl"))
    isdefined(Arrhenius,name) || Base.include(Arrhenius,joinpath(@__DIR__,"..","src",file))
end
directory = ARGS[1]
@testset "Cantera 4 real-gas reaction rates" begin
    for name in ("dodecane","h2-plog")
        gas = CreateSolution(joinpath(directory,name*"_IG.yaml"))
        model = RedlichKwongThermo(joinpath(directory,name*"_RK.yaml"))
        work = RealGasKineticsWorkspace(gas,model)
        if name == "h2-plog"
            @test !isempty(gas.reaction.plog.reaction_indices)
            @test !isempty(gas.reaction.index_three_body)
            @test !isempty(gas.reaction.index_falloff)
        end
        for mode in (name == "dodecane" ? ("RK","IG") : ("RK",))
            data = npzread(joinpath(directory,"rates_"*name*"_"*mode*".npz"))
            for i in eachindex(data["T"])
                T,rho,X = data["T"][i],data["density"][i],data["X"][:,i]
                if mode == "RK"
                    state = redlich_kwong_rates!(work.wdot,gas,model,T,rho,X,work)
                    @test work.activity/(state.P/(R*T)) ≈ data["activities"][:,i] rtol=2e-11
                    @test state.cv_mass ≈ data["cv_mass"][i] rtol=2e-11
                    @test state.u_TV ≈ data["partial_molar_int_energies_TV"][:,i] rtol=2e-11 atol=1e-6
                else
                    ref = species_thermo(gas.thermo,T)
                    work.C .= data["concentrations"][:,i]
                    wdot!(work.wdot,gas.reaction,T,work.C,R*ref.s_R,R*T*ref.h_RT,work.kinetics)
                end
                @test work.C ≈ data["concentrations"][:,i] rtol=2e-11
                for (ours,key) in ((work.kinetics.kf,"forward_rates_of_progress"),
                    (work.kinetics.kr,"reverse_rates_of_progress"),
                    (work.kinetics.equilibrium_constants,"equilibrium_constants"))
                    @test all(isapprox.(ours,data[key][:,i];rtol=3e-10,atol=1e-100))
                end
                @test work.wdot ≈ data["net_production_rates"][:,i] rtol=3e-10 atol=1e-10
                @test abs(dot(work.wdot,gas.MW)) <= 1e-12*sum(abs.(work.wdot.*gas.MW))
                @test maximum(abs,gas.ele_matrix*work.wdot) <= 1e-12*sum(abs,work.wdot)
            end
            println(name," ",mode,": ",length(data["T"])," independent rate states passed")
        end
    end
end
