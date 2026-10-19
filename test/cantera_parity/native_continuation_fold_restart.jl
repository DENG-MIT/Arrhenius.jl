using Arrhenius, NPZ, Printf, LinearAlgebra

np=parse(Int,get(ENV,"CONTINUATION_GRID_POINTS","266"))
fixture=joinpath(@__DIR__,"..","data","cantera_h2o2_two_point_$(np)pt.npz")
ref=NPZ.npzread(fixture)
mechanism=joinpath(dirname(pathof(Arrhenius)),"..","mechanism","h2o2.yaml")
gas=CreateSolution(mechanism)
z=vec(Float64.(ref["grid"]))
f=CounterflowDiffusionFlame(gas;fuel="H2:1",oxidizer="O2:1",
    mdot_fuel=only(ref["turning_fuel_mdot"]),mdot_oxidizer=only(ref["turning_oxidizer_mdot"]),
    T_fuel=300.,T_oxidizer=500.,P=only(ref["pressure_Pa"]),width=.018,grid=z)
n=gas.n_species; T=vec(Float64.(ref["T"])); Y=Matrix{Float64}(ref["Y"])
V=vec(Float64.(ref["spreadRate"])); Lam=vec(Float64.(ref["Lambda"])); velocity=vec(Float64.(ref["velocity"]))
meanMW=1 ./ (Y*(1 ./ gas.MW))
rho=only(ref["pressure_Pa"]) .* meanMW ./ (Arrhenius.R .* T)
u=zeros(n+4,length(z)); u[1,:]=T./1000; u[2:n+1,:]=permutedims(Y)
u[n+2,:]=V./Arrhenius._counterflow_Vscale; u[n+3,:]=Lam./Arrhenius._counterflow_Lscale
u[end,:]=rho.*velocity
f.state=u; f.converged=false
r=similar(u); w=CounterflowWorkspace(f); counterflow_residual!(r,f,u,w)
@printf("restart_seed N=%d mdot=(%.12g,%.12g) Tmax=%.9g initial_residual_inf=%.6g\\n",
    length(z),f.fuel_mass_flux,f.oxidizer_mass_flux,maximum(T),norm(r,Inf))
solve_error=nothing
start=time()
try
    solve!(f;refine_grid=false,max_time_steps=100,max_points=length(z))
catch e
    solve_error=sprint(showerror,e)
end
elapsed=time()-start
counterflow_residual!(r,f,f.state,w)
@printf("native_seed_solve converged=%s elapsed_s=%.6g residual_inf=%.6g Tmax=%.9g points=%d error=%s\\n",
    string(f.converged),elapsed,norm(r,Inf),maximum(temperature(f)),length(f.grid),something(solve_error,"none"))
println("stats=",solver_stats(f))
