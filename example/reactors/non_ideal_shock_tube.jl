# Usage: julia --project=SOLVER_ENV example/reactors/non_ideal_shock_tube.jl MECHANISM_DIR
#        [--smoke] [--finite-difference] [--qndf]
# The caller-owned solver environment needs SciMLBase, OrdinaryDiffEqSDIRK and
# ForwardDiff to activate ArrheniusOrdinaryDiffEqExt, which provides
# non_ideal_shock_tube, shocktube_calculations and shocktube_ignition. The
# optional --qndf route additionally needs OrdinaryDiffEqBDF loaded.
# MECHANISM_DIR contains dodecane_IG.yaml, its native kinetic sidecar, and
# dodecane_RK.yaml. validation/real_gas_reactor_cases.py prepares the paired input
# files and independent Cantera references. No Cantera call occurs in this file.
using Arrhenius, LinearAlgebra, Printf
using SciMLBase, OrdinaryDiffEqSDIRK
BLAS.set_num_threads(1)

if abspath(PROGRAM_FILE) == @__FILE__
    !isempty(ARGS) && all(arg -> arg in ("--smoke","--sweep","--finite-difference","--qndf"),ARGS[2:end]) ||
        error("supply MECHANISM_DIR [--smoke] [--finite-difference] [--qndf]")
    "--qndf" in ARGS && (@eval using OrdinaryDiffEqBDF)
    results = non_ideal_shock_tube(ARGS[1];sweep=!("--smoke" in ARGS),
        jacobian="--finite-difference" in ARGS ? :finite_difference : :ad,
        solver="--qndf" in ARGS ? :qndf : :sdirk)
    for result in results
        @printf("%s T=%4d K: ignition delay %.9g s\n",result.phase,result.temperature,result.delay)
    end
    @printf("1000 K ideal-gas error: %.4g%%\n",100*(results[2].delay/results[1].delay-1))
end
