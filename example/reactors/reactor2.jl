# Native equivalent of Cantera's reactor2.py: pressure-driven piston and heat loss.
# Prepare air/gri30/h2o2 YAML+NPZ once with validation/moving_wall_cases.py.
# julia --project=<SciML environment> example/reactors/reactor2.jl <mechanisms> [output.csv]
include("network_solver.jl")
include("moving_wall_setup.jl")
if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("provide a prepared standard-mechanism directory")
    run_moving_wall_example("reactor2",ARGS[1],get(ARGS,2,"reactor2.csv");integrator=native_network_bdf)
end
