using Arrhenius, NPZ, LinearAlgebra, Test
module NativeTransportUnderTest
    using NPZ, LinearAlgebra, SHA
    using Arrhenius: Solution, cal_cp_R!
    include(joinpath(@__DIR__,"..","src","MulticomponentTransport.jl"))
end
BLAS.set_num_threads(1)
directory = ARGS[1]
worst = Dict("diffusion"=>0.,"thermal_diffusion"=>0.,"conductivity"=>0.,"flux"=>0.,
    "mixture_thermal_diffusion"=>0.,"mixture_diffusion"=>0.)
case_count=0
@testset "Cantera 4 multicomponent coefficients and C++ mass fluxes" begin
    for name in ("h2o2","gri30")
        data=NativeTransportUnderTest.MultiTransportData(joinpath(directory,name*"-transport.npz"))
        w=NativeTransportUnderTest.MultiTransportWorkspace(data)
        wm=NativeTransportUnderTest.MixtureThermalDiffusionWorkspace(data)
        ref=npzread(joinpath(directory,name*"-reference.npz"))
        n=length(data.molecular_weights)
        flux=zeros(n)
        for k in eachindex(ref["T"])
            T,P=ref["T"][k],ref["P"][k]
            X=vec(ref["X"][k,:])
            lambda=NativeTransportUnderTest.multicomponent_transport!(w,data,P,T,X,vec(ref["cp_R"][k,:]))
            NativeTransportUnderTest.mixture_thermal_diffusion!(wm,data,P,T,X)
            NativeTransportUnderTest.multicomponent_fluxes!(flux,w,data,P,T,X,vec(ref["grad_X"][k,:]),ref["grad_T"][k])
            for (key,actual,expected) in (
                ("diffusion",w.diffusion,ref["diffusion"][k,:,:]),
                ("thermal_diffusion",w.thermal_diffusion,ref["thermal_diffusion"][k,:]),
                ("conductivity",[lambda],[ref["conductivity"][k]]),
                ("flux",flux,ref["flux"][k,:]),
                ("mixture_thermal_diffusion",wm.thermal_diffusion,ref["mixture_thermal_diffusion"][k,:]),
                ("mixture_diffusion",wm.diffusion,ref["mixture_diffusion"][k,:]))
                # Normwise absolute floor only for physically zero Soret/flux vectors.
                floor=key=="diffusion" ? 1e-30 : key=="conductivity" ? 1e-30 : 1e-15
                error=maximum(abs.(actual.-expected))/max(maximum(abs.(expected)),floor)
                worst[key]=max(worst[key],error)
                @test error < 2e-8
            end
            @test all(isfinite,w.diffusion)
            @test all(diag(w.diffusion).==0)
            @test lambda>0
            @test abs(sum(w.thermal_diffusion)) < 1e-12
            @test abs(sum(wm.thermal_diffusion)) < 1e-12
            @test abs(sum(flux)) < 1e-12
            global case_count+=1
        end
    end
end
println("Validated states: ",case_count)
for key in sort(collect(keys(worst)))
    println(key," worst normalized error: ",worst[key])
end
