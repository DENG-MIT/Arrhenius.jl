# Native steady burning states at decreasing residence times, using GRI-Mech 3.0.
using Arrhenius
include("network_cases.jl")
include("network_solver.jl")

function main()
  gas = CreateSolution(isempty(ARGS) ? joinpath(@__DIR__,"..","..","mechanism","gri30.yaml") : ARGS[1])
  residence = Ref(0.1)
  network = combustor_network(gas;residence_time=residence)
  state = network_state(network)
  for iteration in 1:100
    result = solve_network_steady(network; integrator=native_network_bdf, initial_state=state,
        interval=max(5residence[],0.01),max_time=max(100residence[],1.0),
        steady_tolerance=1e-7,reltol=1e-8,abstol=1e-15,
        save_everystep=false,save_start=false)
    state = result.state
    chamber = network_diagnostics(network,state).nodes.combustor
    X = Y2X(gas,chamber.mass_fractions)
    heat_release = -sum(cal_h(gas,chamber.temperature,chamber.pressure,X) .*
        set_states(gas,chamber.temperature,chamber.pressure,chamber.mass_fractions))
    println((residence_time_s=residence[],temperature_K=chamber.temperature,
             pressure_Pa=chamber.pressure,heat_release_rate_W_m3=heat_release,residual=result.residual))
    chamber.temperature < 500 && break
    residence[] *= 0.9
  end
end
main()
