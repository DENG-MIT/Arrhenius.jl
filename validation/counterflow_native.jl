using Arrhenius, NPZ, LinearAlgebra
if !isdefined(Arrhenius,:CounterflowDiffusionFlame)
    Base.include(Arrhenius,joinpath(@__DIR__,"..","src","CounterflowFlames.jl"))
end
BLAS.set_num_threads(1)
root=ARGS[1]
kind=length(ARGS)>1 ? ARGS[2] : "h2"
mechanism=joinpath(root,"mechanism",kind=="h2" ? "h2o2.yaml" : "gri30.yaml")
gas=CreateSolution(mechanism)
fuel=kind=="h2" ? "H2:1,AR:1" : kind=="methane" ? "CH4:1" : "C2H6:1"
oxidizer=kind=="h2" ? "O2:.2,AR:.8" : "O2:.21,N2:.78,AR:.01"
f=Arrhenius.CounterflowDiffusionFlame(gas;fuel,oxidizer,mdot_fuel=.24,mdot_oxidizer=.72,
    T_fuel=300.,T_oxidizer=300.,P=101325.,width=.02)
println("Native initial peak T: ",maximum(temperature(f)))
solve!(f;loglevel=2,slope=.2,curve=.3,max_time_steps=800)
r=similar(f.state)
Arrhenius.counterflow_residual!(r,f)
println("Native converged: ",f.converged," points=",length(f.grid)," Tmax=",maximum(temperature(f)),
    " Lambda=",Arrhenius.pressure_curvature(f)[1]," residual=",norm(r,Inf))
npzwrite(joinpath(@__DIR__,"counterflow-"*kind*"-native.npz"),Dict(
    "grid"=>f.grid,"T"=>temperature(f),"Y"=>mass_fractions(f),"velocity"=>velocity(f),
    "spread_rate"=>Arrhenius.spread_rate(f),"Lambda"=>Arrhenius.pressure_curvature(f),
    "state"=>f.state,"residual"=>r))
