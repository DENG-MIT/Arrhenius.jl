"""
Lightweight `SciMLBase`/`OrdinaryDiffEqBDF`/`OrdinaryDiffEqSDIRK` integrator
adapters: `native_network_bdf`, `native_surface_bdf`, `native_surface_flow_bdf`,
`native_engine_sdirk`, `shocktube_sdirk`, `shocktube_integrator`,
`shocktube_qndf_integrator`, `shocktube_ignition`, `shocktube_calculations`,
and `non_ideal_shock_tube`.

Unlike `ArrheniusStructuredReactorExt`, these never touch QNDF's internal
cache interfaces, so their weakdeps are loosely bounded rather than
exact-pinned.
"""
module ArrheniusOrdinaryDiffEqExt

using Arrhenius
using SciMLBase
using OrdinaryDiffEqBDF
using OrdinaryDiffEqSDIRK

include("network_solver.jl")
include("surface_solver.jl")
include("surface_flow_solver.jl")
include("real_gas_ode_solver.jl")
include("real_gas_qndf_solver.jl")
include("ic_engine_solver.jl")
include("non_ideal_shock_tube.jl")

end
