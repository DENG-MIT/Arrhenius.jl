using Arrhenius, LinearAlgebra, NPZ
if !isdefined(Arrhenius,:CatalyticImpingingJet)
    Base.include(Arrhenius,joinpath(@__DIR__,"..","src","CatalyticFlames.jl"))
end
include(joinpath(@__DIR__,"..","example","reactors","surface_solver.jl"))
const CF=Arrhenius
BLAS.set_num_threads(1)
root,output=ARGS[1:2]
mkpath(output)
m=SurfaceMechanism(joinpath(root,"mechanism","pt_h2.surface.npz"))
gas=CreateSolution(m.gas_file)
f=CF.CatalyticImpingingJet(gas,m;reactants="H2:.05,O2:.21,N2:.78,AR:.01",
    mdot=.06,T_inlet=300.,T_surface=900.,width=.1,coverages=Dict("PT(S)"=>.5,"O(S)"=>.5))
CF.initialize_catalytic_coverages!(f;integrator=native_surface_bdf)
npzwrite(joinpath(output,"initial-coverages.npz"),Dict("coverages"=>f.coverages))
println("Initial stationary coverages ",f.coverages)
function save_case(name)
    d=CF.catalytic_wall_diagnostics(f)
    npzwrite(joinpath(output,name*".npz"),Dict("grid"=>f.grid,"state"=>f.state,
        "T"=>CF.temperature(f),"Y"=>CF.mass_fractions(f),"velocity"=>CF.velocity(f),
        "spread_rate"=>CF.spread_rate(f),"Lambda"=>CF.pressure_curvature(f),
        "coverages"=>f.coverages,"coverage_rates"=>d.coverage_rates,
        "surface_production_rates"=>d.surface_production_rates,
        "diffusive_mass_flux"=>d.diffusive_mass_flux,"wall_species_residual"=>d.wall_species_residual,
        "elemental_production"=>d.elemental_production,"total_mass_production"=>d.total_mass_production,
        "surface_heat_release"=>d.surface_heat_release,"flow_residual"=>d.flow_residual))
    println(name," N=",length(f.grid)," Tmax=",maximum(CF.temperature(f))," coverage residual=",
        norm(d.coverage_rates,Inf)," wall residual=",norm(d.wall_species_residual,Inf)," flow residual=",d.flow_residual)
end
CF.set_catalytic_reactions!(f;gas_multiplier=0.,surface_multiplier=0.,coverage_enabled=false)
CF.solve!(f;loglevel=1)
save_case("inert")
get(ENV,"CATALYTIC_STOP","")=="inert" && exit()
for exponent in -5:0
    scale=10.0^exponent
    CF.set_catalytic_reactions!(f;gas_multiplier=scale,surface_multiplier=scale,coverage_enabled=true)
    CF.solve!(f;loglevel=1)
    save_case("h2-"*string(exponent))
end
get(ENV,"CATALYTIC_STOP","")=="h2" && exit()
CF.set_catalytic_inlet!(f,"CH4:.095,O2:.21,N2:.78,AR:.01")
CF.solve!(f;ratio=100.,slope=.15,curve=.2,prune=0.,loglevel=1)
save_case("methane")
