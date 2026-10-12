using Arrhenius

# Cantera's diffusion_flame.py sequence: ethane opposed by air, radiation off/on.
# Preprocess mechanism/gri30.yaml once with mechanism/export_sidecar.py.
gas=CreateSolution(joinpath(@__DIR__,"..","..","mechanism","gri30.yaml"))
flame=CounterflowDiffusionFlame(gas;
    fuel="C2H6:1",oxidizer="O2:.21,N2:.78,AR:.01",
    mdot_fuel=.24,mdot_oxidizer=.72,T_fuel=300.,T_oxidizer=300.,
    P=one_atm,width=.02)
solve!(flame;ratio=4,slope=.2,curve=.3,loglevel=1)
nonradiating_T=temperature(flame)
set_radiation!(flame,true;boundary_emissivities=(0.,0.))
solve!(flame;refine_grid=false,loglevel=1)
println("Peak temperature without radiation: ",maximum(nonradiating_T)," K")
println("Peak temperature with radiation: ",maximum(temperature(flame))," K")
println("Pressure curvature: ",pressure_curvature(flame)[1]," Pa/m²")

output=length(ARGS)>0 ? ARGS[1] : "counterflow_ethane.csv"
T=temperature(flame); U=velocity(flame); V=spread_rate(flame); Y=mass_fractions(flame)
qrad=radiative_heat_loss(flame)
open(output,"w") do io
    println(io,join(vcat(["z_m","T_nonradiating_K","T_radiating_K","u_m_s","V_per_s","radiative_loss_W_m3"],"Y_".*gas.species_names),','))
    for j in eachindex(flame.grid)
        println(io,join(vcat([flame.grid[j],nonradiating_T[j],T[j],U[j],V[j],qrad[j]],Y[:,j]),','))
    end
end
