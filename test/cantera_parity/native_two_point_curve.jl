using Arrhenius, NPZ, Printf, LinearAlgebra
np=parse(Int,get(ENV,"CONTINUATION_GRID_POINTS","12"))
maxsteps=parse(Int,get(ENV,"CONTINUATION_STEPS","250"))
out=get(ENV,"TWO_POINT_OUT",joinpath(@__DIR__,"..","data","native_two_point_$(np)pt.npz"))
fixture=get(ENV,"TWO_POINT_FIXTURE",joinpath(@__DIR__,"..","data","cantera_h2o2_two_point_$(np)pt.npz"))
ref=NPZ.npzread(fixture)
mechanism_name=get(ENV,"CONTINUATION_MECHANISM","h2o2.yaml")
mechanism=isabspath(mechanism_name) ? mechanism_name : joinpath(dirname(pathof(Arrhenius)),"..","mechanism",mechanism_name)
fuel=get(ENV,"CONTINUATION_FUEL","H2:1")
oxidizer=get(ENV,"CONTINUATION_OXIDIZER","O2:1")
T_fuel=parse(Float64,get(ENV,"CONTINUATION_T_FUEL","300"))
T_oxidizer=parse(Float64,get(ENV,"CONTINUATION_T_OXIDIZER","500"))
gas=CreateSolution(mechanism)
P=Float64(only(ref["pressure_Pa"]))
grid=vec(Float64.(ref["initial_grid_m"]))
mdotf=Float64(only(ref["initial_fuel_mdot"]))
mdoto=Float64(only(ref["initial_oxidizer_mdot"]))
flame=CounterflowDiffusionFlame(gas;fuel,oxidizer,mdot_fuel=mdotf,
    mdot_oxidizer=mdoto,T_fuel,T_oxidizer,P,width=grid[end]-grid[1],grid)
Y=Matrix{Float64}(ref["initial_Y"])
T=vec(Float64.(ref["initial_T_K"]))
MWmix=1 ./ (Y*(1 ./ gas.MW))
rho=P .* MWmix ./ (Arrhenius.R .* T)
flame.state[1,:] .= T ./ 1000
flame.state[2:gas.n_species+1,:] .= permutedims(Y)
flame.state[gas.n_species+2,:] .= vec(Float64.(ref["initial_spreadRate"])) ./ Arrhenius._counterflow_Vscale
flame.state[gas.n_species+3,:] .= vec(Float64.(ref["initial_Lambda"])) ./ Arrhenius._counterflow_Lscale
flame.state[end,:] .= rho .* vec(Float64.(ref["initial_velocity_m_s"]))
flame.converged=true
r=similar(flame.state); counterflow_residual!(r,flame,flame.state,CounterflowWorkspace(flame))
@printf("initial shared state N=%d Tmax=%.9g residual_inf=%.6g\n",np,maximum(T),norm(r,Inf))
t0=time()
curve=two_point_control!(flame;max_steps=maxsteps,
    temperature_decrement_schedule=vec(Float64.(ref["curve_temperature_decrement_K"])),
    loglevel=parse(Int,get(ENV,"CONTINUATION_LOGLEVEL","0")))
total=time()-t0
isempty(curve) && error("two-point continuation returned no converged points")
peak=argmax([p.max_strain_rate for p in curve]); p=curve[peak]
@printf("two-point native N=%d points=%d peak_index=%d peak_strain=%.12g Tmax_at_peak=%.12g mdot_f=%.12g mdot_o=%.12g total_s=%.6g sum_solve_s=%.6g\n",
    np,length(curve),peak,p.max_strain_rate,p.max_temperature,p.fuel_mass_flux,p.oxidizer_mass_flux,total,sum(x.solve_seconds for x in curve))
refT=vec(Float64.(ref["Tmax_curve_K"])); refA=vec(Float64.(ref["strain_rate_curve_per_s"]))
@printf("Cantera N=%d points=%d peak_index=%d peak_strain=%.12g Tmax_at_peak=%.12g\n",
    np,length(refT),argmax(refA),maximum(refA),refT[argmax(refA)])
@printf("relative peak strain error=%.9g relative peak Tmax error=%.9g\n",
    abs(p.max_strain_rate-maximum(refA))/maximum(refA),
    abs(p.max_temperature-refT[argmax(refA)])/refT[argmax(refA)])
nmatch=min(length(curve),length(refT))
relT=maximum(abs(curve[i].max_temperature-refT[i])/max(abs(refT[i]),1.) for i in 1:nmatch)
relA=maximum(abs(curve[i].max_strain_rate-refA[i])/max(abs(refA[i]),1.) for i in 1:nmatch)
@printf("pointwise matched=%d native=%d cantera=%d max_relative_Tmax=%.9g max_relative_strain=%.9g\n",
    nmatch,length(curve),length(refT),relT,relA)
let sumYmax=0., element_source_relmax=0.
for q in curve
    Yq=q.state[2:gas.n_species+1,:]
    sumYmax=max(sumYmax,maximum(abs.(sum(Yq;dims=1).-1)))
    for j in axes(Yq,2)
        source=Arrhenius.set_states(gas,1000q.state[1,j],P,@view(Yq[:,j]))
        for e in axes(gas.ele_matrix,1)
            atoms=@view gas.ele_matrix[e,:]
            element_source_relmax=max(element_source_relmax,
                abs(dot(atoms,source))/max(dot(abs.(atoms),abs.(source)),floatmin(Float64)))
        end
    end
end
println("invariants ", (sumY_max=sumYmax, element_source_relmax=element_source_relmax))
end
NPZ.npzwrite(out,Dict("Tmax_curve_K"=>[q.max_temperature for q in curve],
    "strain_rate_curve_per_s"=>[q.max_strain_rate for q in curve],
    "fuel_mdot"=>[q.fuel_mass_flux for q in curve],"oxidizer_mdot"=>[q.oxidizer_mass_flux for q in curve],
    "control_temperature_K"=>[q.control_temperature for q in curve],
    "temperature_decrement_K"=>[q.temperature_decrement for q in curve],
    "solve_seconds"=>[q.solve_seconds for q in curve],"total_seconds"=>[total],
    "turning_index"=>[peak],"turning_T"=>p.state[1,:].*1000,
    "turning_Y"=>permutedims(p.state[2:gas.n_species+1,:]),
    "turning_velocity"=>p.state[end,:]./(P .* (1 ./ (permutedims(p.state[2:gas.n_species+1,:])*(1 ./ gas.MW))) ./ (Arrhenius.R .* (p.state[1,:].*1000))),"turning_spreadRate"=>p.state[gas.n_species+2,:].*Arrhenius._counterflow_Vscale,
    "turning_Lambda"=>p.state[gas.n_species+3,:].*Arrhenius._counterflow_Lscale,
    "grid"=>copy(flame.grid)))
