using Arrhenius, LinearAlgebra, NPZ
if !isdefined(Arrhenius,:CounterflowPremixedFlame)
    Base.include(Arrhenius,joinpath(@__DIR__,"..","src","PremixedCounterflowFlames.jl"))
end
const CF=Arrhenius
BLAS.set_num_threads(1)
root,output=ARGS[1:2]
gas=CreateSolution(joinpath(root,"mechanism","h2o2.yaml"))
f=CF.ImpingingJet(gas;reactants="H2:1.8,O2:1,AR:7",mdot=.06,T_inlet=373.,
    T_surface=500.,P=.05*one_atm,width=.2)
for mdot in (.06,.07,.08,.09,.10,.11,.12)
    CF.set_mass_flux!(f;reactants=mdot)
    CF.solve!(f;ratio=3,slope=.1,curve=.2,prune=.06,grid_min=1e-4,loglevel=1)
    r=similar(f.state); CF.counterflow_residual!(r,f)
    npzwrite(joinpath(output,"counterflow-wall-$(round(Int,100mdot))-native.npz"),
        Dict("grid"=>f.grid,"T"=>CF.temperature(f),"Y"=>CF.mass_fractions(f),
        "velocity"=>CF.velocity(f),"spread_rate"=>CF.spread_rate(f),
        "Lambda"=>CF.pressure_curvature(f),"residual"=>r,"mdot"=>mdot))
    println("Wall sequence mdot=",mdot," points=",length(f.grid)," Tmax=",maximum(CF.temperature(f)),
        " Lambda=",CF.pressure_curvature(f)[1]," residual=",norm(r,Inf))
end
