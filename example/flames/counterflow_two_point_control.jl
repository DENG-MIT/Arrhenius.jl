using Arrhenius

mechanism=joinpath(dirname(pathof(Arrhenius)),"..","mechanism","h2o2.yaml")
gas=CreateSolution(mechanism)
grid=collect(range(0.,.018;length=48))
flame=CounterflowDiffusionFlame(gas;fuel="H2:1",oxidizer="O2:1",
    mdot_fuel=.5,mdot_oxidizer=3.,T_fuel=300.,T_oxidizer=500.,
    P=one_atm,width=.018,grid)
solve!(flame;refine_grid=false)
curve=two_point_control!(flame;max_steps=100)
peak=argmax(point.max_strain_rate for point in curve)
point=curve[peak]
println("peak strain = $(point.max_strain_rate) 1/s")
println("Tmax at peak = $(point.max_temperature) K")
