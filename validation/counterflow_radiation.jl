using Arrhenius, LinearAlgebra, NPZ

const CF=Arrhenius
BLAS.set_num_threads(1)
root,output=ARGS[1:2]
emissivities=length(ARGS)>2 ? Tuple(parse.(Float64,split(ARGS[3],','))) : (0.,0.)
gas=CreateSolution(joinpath(root,"mechanism","gri30.yaml"))
f=CF.CounterflowDiffusionFlame(gas;fuel="C2H6:1",oxidizer="O2:.21,N2:.78,AR:.01",
    mdot_fuel=.24,mdot_oxidizer=.72,width=.02,boundary_emissivities=emissivities)
CF.solve!(f;ratio=4,slope=.2,curve=.3,loglevel=1)
nonradiating_T=CF.temperature(f)
nonradiating_Lambda=CF.pressure_curvature(f)
grid=copy(f.grid)
CF.set_radiation!(f,true)
CF.solve!(f;refine_grid=false,loglevel=1)
@assert f.grid==grid
r=similar(f.state)
CF.counterflow_residual!(r,f)
sources=CF.radiation_source(f)
npzwrite(output,Dict("grid"=>f.grid,"T"=>CF.temperature(f),"Y"=>CF.mass_fractions(f),
    "velocity"=>CF.velocity(f),"spread_rate"=>CF.spread_rate(f),"Lambda"=>CF.pressure_curvature(f),
    "state"=>f.state,"residual"=>r,"radiative_heat_loss"=>sources.heat_loss,
    "planck_absorption"=>sources.planck_absorption,"nonradiating_T"=>nonradiating_T,
    "nonradiating_Lambda"=>nonradiating_Lambda,"boundary_emissivities"=>collect(emissivities)))
println("Radiating ethane: points=",length(f.grid)," nonradiating Tmax=",maximum(nonradiating_T),
    " radiating Tmax=",maximum(CF.temperature(f))," maxloss=",maximum(sources.heat_loss),
    " residual=",norm(r,Inf))

# Optional isolated property check on reference states. These states are never
# supplied to the native nonlinear solve above.
if length(ARGS)>3
    using Test
    ref=npzread(ARGS[4])
    probe=CF.CounterflowDiffusionFlame(gas;fuel="C2H6:1",oxidizer="O2:.21,N2:.78,AR:.01",
        mdot_fuel=.24,mdot_oxidizer=.72,grid=ref["grid"],radiation=true,
        boundary_emissivities=Tuple(ref["boundary_emissivities"]))
    probe.state[1,:] .= ref["T"]./1000
    probe.state[2:gas.n_species+1,:] .= ref["Y"]
    q=CF.radiative_heat_loss(probe)
    normalized=maximum(abs.(q.-ref["radiative_heat_loss"]))/maximum(abs.(ref["radiative_heat_loss"]))
    @test normalized<1e-12
    println("Radiation source isolated normalized error: ",normalized)
end
