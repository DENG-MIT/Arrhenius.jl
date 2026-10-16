# Prepare ptcombust.yaml / Pt_surf with mechanism/export_surface.py first.
using Arrhenius
include("surface_solver.jl")

function catalytic_reactor(path;flow=false)
    surface = SurfaceMechanism(path)
    gas = CreateSolution(surface.gas_file)
    initial = IdealGasReactor(gas;temperature=900.0,pressure=one_atm,
        mole_fractions=Dict("H2"=>0.05,"O2"=>0.21,"N2"=>0.78,"AR"=>0.01),
        constraint=:constant_volume,energy=:isothermal)
    vessel = WellStirredReactor(initial;volume=1e-6)
    if flow
        feed = MassFlowController(:inlet,:reactor;mdot=initial.density*vessel.volume/0.001)
        outlet = PressureController(:reactor,:exhaust;primary=feed,K=1e-7)
        network = ReactorNetwork((inlet=Reservoir(initial),reactor=vessel,exhaust=Reservoir(initial));
                                 flows=(feed,outlet))
    else
        network = ReactorNetwork((reactor=vessel,))
    end
    system = CatalyticNetwork(network;surfaces=(ReactorSurface(:reactor,surface;area=0.001,
                              coverages=Dict("PT(S)"=>0.5,"O(S)"=>0.5)),))
    return system
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 1 || error("supply the prepared Pt surface parameter archive")
    system = catalytic_reactor(ARGS[1];flow="--flow" in ARGS)
    solution = solve_catalytic(system,(0.0,0.02);integrator=native_surface_bdf,
        reltol=1e-10,abstol=1e-18,saveat=0.0001,maxiters=1_000_000)
    report = catalytic_diagnostics(system,solution.u[end],solution.t[end])
    println("Combined gas/surface mass: ",report.mass," kg")
    println("Combined elemental inventory: ",Dict(zip(report.element_names,report.element_inventory)))
    println("Thermostat power: ",report.thermostat_power," W")
end
