# Native equivalent of Cantera's custom2.py, including integrated piston inertia.
using Arrhenius
using SciMLBase, OrdinaryDiffEqBDF
if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("provide a prepared standard-mechanism directory")
    run_moving_wall_example("custom2",ARGS[1],get(ARGS,2,"custom2.csv");integrator=native_network_bdf)
end
