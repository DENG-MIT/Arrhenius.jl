# Prepare once: python mechanism/export_surface.py diamond.yaml diamond_100 --output diamond.surface.npz
# Run: julia --project=YOUR_SCIML_ENV example/reactors/diamond_cvd.jl diamond.surface.npz
using Arrhenius
include("surface_solver.jl")

function diamond_cvd(path)
    m = SurfaceMechanism(path)
    X = copy(m.initial_mole_fractions)
    kH = findfirst(==("H"),m.species_names[m.n_surface+1:m.n_surface+m.n_gas])
    theta = copy(m.initial_coverages)
    results = NamedTuple[]
    println("atomic_hydrogen_mole_fraction,growth_micrometres_per_hour")
    for _ in 1:20
        surface = IdealSurface(m;temperature=1200.0,pressure=20one_atm/760,
                               mole_fractions=X,coverages=theta)
        stationary = steady_coverages(surface;integrator=native_surface_bdf,
            reltol=1e-10,abstol=1e-18,maxiters=1_000_000)
        theta = stationary.coverages
        rates = surface_rates(surface,theta)
        growth = only(rates.bulk)*only(m.bulk_molar_volumes)*1e6*3600
        push!(results,(hydrogen=X[kH],growth=growth,coverages=copy(theta)))
        println(X[kH],",",growth)
        X[kH] /= 1.4
        X ./= sum(X)
    end
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("supply the prepared diamond surface parameter archive")
    diamond_cvd(ARGS[1])
end
