using Arrhenius, NPZ, Printf
using LinearAlgebra

np=parse(Int,get(ENV,"CONTINUATION_GRID_POINTS","266"))
steps=parse(Int,get(ENV,"CONTINUATION_STEPS","140"))
parameter_scale=parse(Float64,get(ENV,"CONTINUATION_PARAMETER_SCALE","1"))
state_scale=parse(Float64,get(ENV,"CONTINUATION_STATE_SCALE","1000"))
mechanism=joinpath(dirname(pathof(Arrhenius)),"..","mechanism","h2o2.yaml")
fixture=joinpath(@__DIR__,"..","data","cantera_h2o2_two_point_$(np)pt.npz")
ref=NPZ.npzread(fixture)
gas=CreateSolution(mechanism)
width=.018
flame=CounterflowDiffusionFlame(gas;fuel="H2:1",oxidizer="O2:1",
    mdot_fuel=.5,mdot_oxidizer=3.,T_fuel=300.,T_oxidizer=500.,P=1e5,width,
    grid=collect(range(0,width;length=np)))
solve!(flame;refine_grid=false,max_time_steps=200,max_points=np)
curve=pseudo_arclength!(flame;parameter=:inlet_flux_scale,steps,ds=.05,
    min_ds=.05/16,max_ds=.1,tolerance=1e-8,parameter_scale,state_scale)
peak=argmax([p.max_strain_rate for p in curve])
p=curve[peak]
@printf("grid=%d points=%d parameter_start=%.9g parameter_end=%.9g native peak strain=%.9g Tmax=%.9g residual=%.3g turning_flag=%s\n",
    np,length(curve),curve[1].parameter,curve[end].parameter,p.max_strain_rate,p.max_temperature,p.residual_inf,string(p.turning_point))
@printf("Cantera fixed-grid two-point peak strain=%.9g Tmax=%.9g turn_index=%d fuel_mdot=%.9g oxidizer_mdot=%.9g\n",
    only(ref["turning_strain_rate_per_s"]),only(ref["turning_Tmax_K"]),only(ref["turning_index"]),only(ref["turning_fuel_mdot"]),only(ref["turning_oxidizer_mdot"]))
@printf("relative peak strain difference=%.6g relative Tmax difference=%.6g\n",
    abs(p.max_strain_rate-only(ref["turning_strain_rate_per_s"]))/only(ref["turning_strain_rate_per_s"]),
    abs(p.max_temperature-only(ref["turning_Tmax_K"]))/only(ref["turning_Tmax_K"]))
@printf("max residual over native curve=%.3g\n",maximum(x.residual_inf for x in curve))
let sumYmax=0., element_source_relmax=0., mass_flux_bc_max=0.
for x in curve
    Y=x.state[2:gas.n_species+1,:]
    sumYmax=max(sumYmax,maximum(abs.(sum(Y;dims=1).-1)))
    mass_flux_bc_max=max(mass_flux_bc_max,abs(x.state[end,1]-.5*x.parameter),abs(x.state[end,end]+3*x.parameter))
    for j in axes(Y,2)
        source=Arrhenius.set_states(gas,1000*x.state[1,j],1e5,@view(Y[:,j]))
        for e in axes(gas.ele_matrix,1)
            atoms=@view gas.ele_matrix[e,:]
            element_source_relmax=max(element_source_relmax,abs(dot(atoms,source))/max(dot(abs.(atoms),abs.(source)),floatmin(Float64)))
        end
    end
end
@printf("invariants sumY_max=%.3g element_source_relmax=%.3g mass_flux_boundary_error=%.3g\n",sumYmax,element_source_relmax,mass_flux_bc_max)
end
