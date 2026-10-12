using Arrhenius, LinearAlgebra, NPZ, Test
if !isdefined(Arrhenius,:DustyGasTransport)
    Base.include(Arrhenius,joinpath(@__DIR__,"..","src","DustyGasTransport.jl"))
end
BLAS.set_num_threads(1)
mechanism,multisidecar,reference,output=ARGS[1:4]
gas=CreateSolution(mechanism)
data=MultiTransportData(multisidecar,gas)
w=DustyGasTransport(gas;porosity=.2,tortuosity=4.,mean_pore_radius=1.5e-7,
    mean_particle_diameter=1.5e-6,multicomponent_data=data)
ref=npzread(reference)
errors=zeros(3,length(ref["T"]))
normalized(a,b)=maximum(abs.(a.-b))/max(maximum(abs.(b)),1e-300)
@testset "Cantera4 dusty gas transport" begin
    @test split(String(ref["species_utf8"]),'\n')==gas.species_names
    for i in eachindex(ref["T"])
        por,tort,radius,diameter,permeability=ref["medium"][:,i]
        set_porous_medium!(w;porosity=por,tortuosity=tort,mean_pore_radius=radius,
            mean_particle_diameter=diameter,permeability=permeability<0 ? nothing : permeability)
        P,T,X=ref["P"][i],ref["T"][i],ref["X"][:,i]
        D=dusty_gas_diffusion!(w,P,T,X)
        lambda=dusty_gas_thermal_conductivity(w,P,T,X)
        errors[1,i]=normalized(D,ref["diffusion"][:,:,i])
        errors[2,i]=abs(lambda/ref["conductivity"][i]-1)
        @test errors[1,i]<1e-12
        @test errors[2,i]<1e-12
        flux=dusty_gas_molar_fluxes(w,T,ref["T2"][i],ref["rho1"][i],ref["rho2"][i],
            ref["Y1"][:,i],ref["Y2"][:,i],ref["delta"][i])
        errors[3,i]=normalized(flux,ref["flux"][:,i])
        @test errors[3,i]<1e-12
        @test normalized(w.resistance*flux,-w.rhs)<1e-12
    end
end
npzwrite(output,Dict("errors"=>errors))
println("Maximum relative errors (diffusion, conductivity, flux): ",vec(maximum(errors;dims=2)))

source=Dict("OH"=>1.,"H"=>2.,"O2"=>3.,"O"=>1e-8,"H2"=>1e-8,"H2O"=>1e-8,
    "H2O2"=>1e-8,"HO2"=>1e-8,"AR"=>1e-8)
X=mole_fractions(gas,source)
Y=X.*gas.MW./dot(X,gas.MW)
rho1=one_atm*dot(X,gas.MW)/(R*500)
rho2=1.2*rho1
set_porous_medium!(w;porosity=.2,tortuosity=4.,mean_pore_radius=1.5e-7,
    mean_particle_diameter=1.5e-6,permeability=nothing)
D=copy(dusty_gas_diffusion!(w,one_atm,500.,X))
lambda=dusty_gas_thermal_conductivity(w,one_atm,500.,X)
zero=dusty_gas_molar_fluxes(w,500.,500.,rho1,rho1,Y,Y,.001)
pressure=dusty_gas_molar_fluxes(w,500.,500.,rho1,rho2,Y,Y,.001)
expected=npzread(replace(reference,".npz"=>"-example.npz"))
@testset "exact source dusty_gas.py sequence" begin
    @test normalized(D,expected["diffusion"])<1e-12
    @test abs(lambda/expected["conductivity"]-1)<1e-12
    @test all(iszero,zero)
    @test normalized(pressure,expected["pressure_flux"])<1e-12
    # Native conductivity must remain independent of prior flux temperatures.
    dusty_gas_molar_fluxes(w,500.,585.,rho1,rho2,Y,Y,.001)
    dusty_gas_diffusion!(w,one_atm,500.,X)
    @test dusty_gas_thermal_conductivity(w,one_atm,500.,X)==lambda
    set_porous_medium!(w;porosity=.4,tortuosity=2.)
    @test dusty_gas_thermal_conductivity(w,10one_atm,500.,X)==lambda
end
println("Exact example conductivity=",lambda," pressure flux=",pressure)
