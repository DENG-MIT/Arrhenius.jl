# Native execution of the official input values and transport continuation.
using Arrhenius,LinearAlgebra,NPZ
BLAS.set_num_threads(1)
parameters,output,case=ARGS[1:3]
free=case=="free";fixed=case=="fixed"
name=fixed ? "gri30" : "h2o2"
mechanism=joinpath(parameters,name*".yaml");mkpath(output)
gas=CreateSolution(mechanism);data=MultiTransportData(mechanism*".multicomponent.npz",gas)
levels=parse(Int,get(ENV,"CONSERVATIVE_LEVELS",free ? "3" : fixed ? "2" : "5"))
width=parse(Float64,get(ENV,"CONSERVATIVE_FREE_WIDTH",".03"))
X=free ? "H2:1.1,O2:1,AR:5" : fixed ? "CH4:.65,O2:1,N2:3.76" : "H2:1.5,O2:1,AR:7"
f=free ? FreeFlame(gas;T=300.,P=one_atm,X,width,discretization=:conservative) :
    BurnerFlame(gas;T=fixed ? 373.7 : 373.,P=fixed ? one_atm : .05one_atm,X,
        mdot=fixed ? .04 : .06,width=fixed ? .01 : .5,discretization=:conservative)
if fixed
    profile=npzread(joinpath(parameters,"fixed-profile.npz"))
    set_temperature_profile!(f,vec(profile["positions"]),vec(profile["temperatures"]);relative=false)
end
atomic=gas.ele_matrix'\gas.MW;E=Diagonal(atomic)*gas.ele_matrix*Diagonal(1 ./gas.MW)
modes=free ? ["mass","mass-soret","multi","multi-soret"] : ["mole","multi"]
for mode in modes
    multi=startswith(mode,"multi")
    set_transport!(f,multi ? :multicomponent : :mixture_averaged;data,soret=endswith(mode,"soret"),
        flux_gradient_basis=free ? :mass : :mole)
    slope=free ? .06 : fixed ? (multi ? .1 : .3) : .05
    curve=free ? .12 : fixed ? (multi ? .2 : 1.) : .1
    println("Official sequence ",case," ",mode," species=",gas.n_species);flush(stdout)
    solve!(f;slope,curve,max_points=3000,max_time_steps=600,loglevel=1)
    refined=deepcopy(f)
    for level in 0:levels
        if level>0
            old=refined.state;N=length(refined.grid)
            refined.grid=sort(vcat(refined.grid,.5 .* (refined.grid[1:end-1].+refined.grid[2:end])))
            refined.state=zeros(size(old,1),2N-1)
            refined.state[:,1:2:end].=old;refined.state[:,2:2:end].=.5 .* (old[:,1:end-1].+old[:,2:end])
            refined.anchor=2refined.anchor-1
            if fixed
                refined.imposed_temperature=[Arrhenius._interpolate_profile(refined.profile_positions,refined.profile_temperatures,z) for z in refined.grid]
                refined.state[1,:].=refined.imposed_temperature./1000
            end
            solve!(refined;refine_grid=false,max_time_steps=600)
        end
        w=Arrhenius.FlameWorkspace(refined)
        Arrhenius._flame_newton!(refined,w;tolerance=1e-9,maxiters=50)
        r=flame_residual!(similar(refined.state),refined,refined.state,w);c=w.conservative
        norm(r,Inf)<1e-8 || error("residual exceeds production tolerance")
        mdot=refined.state[end,1];Y=mass_fractions(refined);T=temperature(refined)
        println((;case,mode,level,points=length(refined.grid),speed=flame_speed(refined),Tmax=maximum(T),
            element_outlet_error=maximum(abs,E*(Y[:,end]-refined.inlet_Y)),
            element_flux_error=maximum(abs,E*c.species_flux.-mdot.*(E*refined.inlet_Y))/mdot,residual=norm(r,Inf)));flush(stdout)
        npzwrite(joinpath(output,"$case-$mode-$level.npz"),Dict("grid"=>refined.grid,"T"=>T,"Y"=>Y,"inlet_Y"=>refined.inlet_Y,
            "state"=>refined.state,"velocity"=>velocity(refined),"P"=>[refined.pressure],"element_matrix"=>E,"residual"=>r,
            "species_flux"=>c.species_flux,"enthalpy_flux"=>c.enthalpy_flux,"diffusive_flux"=>w.flux,"conductivity"=>w.conductivity))
    end
end
