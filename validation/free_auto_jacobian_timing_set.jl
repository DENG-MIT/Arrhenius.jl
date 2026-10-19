length(ARGS) > 0 && length(ARGS) % 3 == 0 || error(
    "supply one or more triples: saved_state.npz mechanism.yaml pressure_atm")
cases = copy(ARGS)
for start in 1:3:length(cases)
    empty!(ARGS)
    append!(ARGS, cases[start:start+2])
    include(joinpath(@__DIR__, "..", "validation", "free_auto_jacobian_timing.jl"))
end
