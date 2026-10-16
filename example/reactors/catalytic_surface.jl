# Prepare once: python mechanism/export_surface.py ptcombust.yaml Pt_surf --output pt.surface.npz
# Run: julia --project=YOUR_SCIML_ENV example/reactors/catalytic_surface.jl pt.surface.npz
using Arrhenius
include("surface_solver.jl")

function catalytic_surface(path)
    mechanism = SurfaceMechanism(path)
    surface = IdealSurface(mechanism; temperature=900.0, pressure=one_atm,
        mole_fractions=Dict("H2"=>0.05,"O2"=>0.21,"N2"=>0.78,"AR"=>0.01),
        coverages=Dict("PT(S)"=>0.5,"O(S)"=>0.5))
    result = steady_coverages(surface; integrator=native_surface_bdf,
        reltol=1e-10, abstol=1e-18, maxiters=1_000_000)
    rates = surface_rates(surface,result.coverages)
    println("species,coverage")
    for (name, theta) in zip(mechanism.species_names[1:mechanism.n_surface],result.coverages)
        println(name,",",theta)
    end
    println("Maximum stationary coverage rate: ",result.residual," s^-1")
    return (;surface,result,rates)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("supply the prepared Pt surface parameter archive")
    catalytic_surface(ARGS[1])
end
