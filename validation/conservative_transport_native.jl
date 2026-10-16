# Direct native initializations followed by controlled equal-domain refinement.
using Arrhenius,LinearAlgebra,NPZ
BLAS.set_num_threads(1)
root,output=ARGS[1:2];mkpath(output)
cases=length(ARGS)>2 ? ARGS[3:end] : ["free"]
levels=parse(Int,get(ENV,"CONSERVATIVE_LEVELS","3"))
burner_mdot=parse(Float64,get(ENV,"CONSERVATIVE_BURNER_MDOT","0.06"))
modes=split(get(ENV,"CONSERVATIVE_MODES","mole,mass,mole-soret,mass-soret,multi,multi-soret"),",")
gas=CreateSolution(joinpath(root,"mechanism","h2o2.yaml"))
data=MultiTransportData(joinpath(root,"mechanism","h2o2.yaml.multicomponent.npz"),gas)
atomic=gas.ele_matrix'\gas.MW
E=Diagonal(atomic)*gas.ele_matrix*Diagonal(1 ./gas.MW)
for case in cases, mode in modes
    free=case=="free"
    kwargs=(;discretization=:conservative,transport_model=startswith(mode,"multi") ? :multicomponent : :mixture_averaged,
        multicomponent_data=data,soret=endswith(mode,"soret"),flux_gradient_basis=startswith(mode,"mass") ? :mass : :mole)
    f=free ? FreeFlame(gas;T=300.,P=one_atm,X="H2:1.1,O2:1,AR:5",width=.06,kwargs...) :
        BurnerFlame(gas;T=373.,P=.05one_atm,mdot=burner_mdot,width=.5,X="H2:1.5,O2:1,AR:7",kwargs...)
    case=="fixed" && set_temperature_profile!(f,[0.,.005,.01,.02,.05,.1,1.],[373.,650.,1000.,1350.,1650.,1750.,1750.])
    for initial_level in 1:parse(Int,get(ENV,"CONSERVATIVE_INITIAL_BISECTIONS","0"))
        old=f.state;N=length(f.grid)
        f.grid=sort(vcat(f.grid,.5 .* (f.grid[1:end-1].+f.grid[2:end])))
        f.state=zeros(size(old,1),2N-1)
        f.state[:,1:2:end].=old;f.state[:,2:2:end].=.5 .* (old[:,1:end-1].+old[:,2:end])
        f.anchor=2f.anchor-1
        if case=="fixed"
            f.imposed_temperature=[Arrhenius._interpolate_profile(f.profile_positions,f.profile_temperatures,z) for z in f.grid]
            f.state[1,:].=f.imposed_temperature./1000
        end
    end
    println("Direct ",case," ",mode);flush(stdout)
    try
        solve!(f;slope=free ? .06 : .05,curve=free ? .12 : .1,max_points=3000,max_time_steps=600,
            loglevel=parse(Int,get(ENV,"CONSERVATIVE_LOGLEVEL","0")))
    catch
        npzwrite(joinpath(output,"failed-$case-$mode.npz"),Dict("grid"=>f.grid,"state"=>f.state,"anchor"=>[f.anchor],
            "fixed_temperature"=>[f.fixed_temperature]))
        rethrow()
    end
    for level in 0:levels
        if level>0
            old=f.state;N=length(f.grid)
            f.grid=sort(vcat(f.grid,.5 .* (f.grid[1:end-1].+f.grid[2:end])))
            f.state=zeros(size(old,1),2N-1)
            f.state[:,1:2:end].=old;f.state[:,2:2:end].=.5 .* (old[:,1:end-1].+old[:,2:end])
            f.anchor=2f.anchor-1
            if case=="fixed"
                f.imposed_temperature=[Arrhenius._interpolate_profile(f.profile_positions,f.profile_temperatures,z) for z in f.grid]
                f.state[1,:].=f.imposed_temperature./1000
            end
            try
                solve!(f;refine_grid=false,max_time_steps=600)
            catch
                r=flame_residual!(similar(f.state),f)
                println((;case,mode,level,failure_residual=norm(r,Inf),location=argmax(abs.(r)),minY=minimum(f.state[2:end-1,:])))
                npzwrite(joinpath(output,"failed-$case-$mode-$level.npz"),Dict("grid"=>f.grid,"state"=>f.state,
                    "anchor"=>[f.anchor],"fixed_temperature"=>[f.fixed_temperature],"residual"=>r))
                rethrow()
            end
        end
        w=Arrhenius.FlameWorkspace(f)
        polished=Arrhenius._flame_newton!(f,w;tolerance=1e-9,maxiters=50)
        r=flame_residual!(similar(f.state),f,f.state,w)
        norm(r,Inf)<1e-8 || error("steady residual failed after polish")
        c=w.conservative;mdot=f.state[end,1]
        element_flux_error=maximum(abs,E*c.species_flux .- mdot.*(E*f.inlet_Y))/mdot
        Y=mass_fractions(f);T=temperature(f)
        println((;case,mode,level,points=length(f.grid),speed=flame_speed(f),Tmax=maximum(T),
            element_drift=maximum(abs,E*(Y[:,end]-f.inlet_Y)),element_flux_error,
            enthalpy_flux_range=maximum(c.enthalpy_flux)-minimum(c.enthalpy_flux),residual=norm(r,Inf),polished));flush(stdout)
        npzwrite(joinpath(output,"$case-$mode-$level.npz"),Dict("grid"=>f.grid,"T"=>T,"Y"=>Y,
            "state"=>f.state,"inlet_Y"=>f.inlet_Y,"velocity"=>velocity(f),"element_matrix"=>E,"P"=>[f.pressure],
            "residual"=>r,"species_flux"=>c.species_flux,"enthalpy_flux"=>c.enthalpy_flux,
            "diffusive_flux"=>w.flux,"conductivity"=>w.conductivity))
    end
end
