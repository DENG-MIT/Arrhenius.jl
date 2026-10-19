using Arrhenius

mechanism=joinpath(@__DIR__,"..","..","mechanism","h2o2.yaml")
gas=CreateSolution(mechanism)
flame=CounterflowDiffusionFlame(gas;fuel="H2:1",oxidizer="O2:1",
    mdot_fuel=.5,mdot_oxidizer=3.,T_oxidizer=500.,width=.018,
    grid=collect(range(0,.018;length=12)))
solve!(flame;refine_grid=false,max_time_steps=100)
curve=pseudo_arclength!(flame;parameter=:inlet_flux_scale,steps=20,ds=.05,parameter_scale=1.,state_scale=1000.)
for p in curve
    println((parameter=p.parameter,Tmax_K=p.max_temperature,
        strain_rate_per_s=p.max_strain_rate,residual_inf=p.residual_inf,
        turning_point=p.turning_point))
end
