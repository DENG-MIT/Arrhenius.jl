# Native equivalent of Cantera's piston.py: two reacting gases and release at 0.1 s.
include("network_solver.jl")
include("moving_wall_setup.jl")
if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("provide a prepared standard-mechanism directory")
    run_moving_wall_example("piston",ARGS[1],get(ARGS,2,"piston.csv");integrator=native_network_bdf)
end
