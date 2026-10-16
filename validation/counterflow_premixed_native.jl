using Arrhenius, LinearAlgebra, NPZ
if !isdefined(Arrhenius,:CounterflowPremixedFlame)
    Base.include(Arrhenius,joinpath(@__DIR__,"..","src","PremixedCounterflowFlames.jl"))
end
const CF=Arrhenius
BLAS.set_num_threads(1)
root,kind,output=ARGS[1:3]
gas=CreateSolution(joinpath(root,"mechanism",kind=="twin" ? "gri30.yaml" : "h2o2.yaml"))
if kind=="premixed"
    f=CF.CounterflowPremixedFlame(gas;reactants="H2:1.6,O2:1,AR:7",
        mdot_reactants=.12,mdot_products=.06,T_reactants=373.,P=.05*one_atm,width=.2)
elseif kind=="twin"
    composition="CH4:.75,O2:2,N2:7.52"
    X=mole_fractions(gas,composition)
    mdot=2*one_atm*dot(X,gas.MW)/(R*300.)
    f=CF.CounterflowTwinPremixedFlame(gas;reactants=composition,mdot,width=.025)
elseif kind=="wall"
    f=CF.ImpingingJet(gas;reactants="H2:1.8,O2:1,AR:7",mdot=.06,T_inlet=373.,
        T_surface=500.,P=.05*one_atm,width=.2)
else
    error("kind must be premixed, twin or wall")
end
CF.solve!(f;ratio=kind=="twin" ? 2 : 3,slope=kind=="twin" ? .3 : .1,
    curve=kind=="twin" ? .3 : .2,prune=kind=="twin" ? .05 : kind=="wall" ? .06 : .02,
    grid_min=kind=="wall" ? 1e-4 : 1e-10,loglevel=1)
r=similar(f.state)
CF.counterflow_residual!(r,f)
npzwrite(output,Dict("grid"=>f.grid,"T"=>CF.temperature(f),"Y"=>CF.mass_fractions(f),
    "velocity"=>CF.velocity(f),"spread_rate"=>CF.spread_rate(f),"Lambda"=>CF.pressure_curvature(f),
    "state"=>f.state,"residual"=>r))
println(kind," points=",length(f.grid)," Tmax=",maximum(CF.temperature(f)),
    " Lambda=",CF.pressure_curvature(f)[1]," residual=",norm(r,Inf))
if kind=="twin"
    diagnostics=CF.twin_flame_diagnostics(f)
    npzwrite(replace(output,".npz"=>"-diagnostics.npz"),Dict(
        "consumption_speed"=>diagnostics.consumption_speed,
        "characteristic_strain_rate"=>diagnostics.characteristic_strain_rate,
        "strain_rate_point"=>diagnostics.strain_rate_point,
        "strain_rate_profile"=>diagnostics.strain_rate_profile,
        "density"=>CF.density(f),"heat_release_rate"=>CF.heat_release_rate(f)))
    println("Twin consumption speed=",diagnostics.consumption_speed,
        " characteristic strain rate=",diagnostics.characteristic_strain_rate)
end
