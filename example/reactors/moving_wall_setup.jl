using Arrhenius

"Build the published moving-wall examples using prepared standard Cantera mechanisms."
function moving_wall_example(name, mechanism_directory)
    gas(file) = CreateSolution(joinpath(mechanism_directory,file*".yaml"))
    if name == "reactor2"
        argon = WellStirredReactor(gas("air"); temperature=1000.0, pressure=20one_atm,
                                  mole_fractions=Dict("AR"=>1.0))
        reacting = WellStirredReactor(gas("gri30"); temperature=500.0, pressure=0.2one_atm,
            mole_fractions=Dict("CH4"=>1.1,"O2"=>2.0,"N2"=>7.52))
        environment = Reservoir(gas("air");temperature=300.0,pressure=one_atm,
            mole_fractions=Dict("O2"=>0.21,"N2"=>0.78,"AR"=>0.01))
        base = ReactorNetwork((argon=argon,reacting=reacting,environment=environment);
            walls=(HeatTransferWall(:reacting,:environment;U=500.0),))
        model = MovingWallNetwork(base;walls=(MovingWall(:reacting,:argon;K=0.5e-4,U=100.0),))
        return (; model, times=collect(range(0,0.12;length=301)), stops=Float64[])
    elseif name == "piston"
        left = WellStirredReactor(gas("h2o2");temperature=900.0,pressure=one_atm,volume=0.5,
            mole_fractions=Dict("H2"=>2.0,"O2"=>1.0,"AR"=>20.0))
        right = WellStirredReactor(gas("gri30");temperature=900.0,pressure=one_atm,volume=0.1,
            mole_fractions=Dict("CO"=>2.0,"H2O"=>0.01,"O2"=>5.0))
        velocity = (states,t) -> t < 0.1 ? 0.0 : 1e-4*(states.left.pressure-states.right.pressure)
        model = MovingWallNetwork(ReactorNetwork((left=left,right=right));
            walls=(MovingWall(:left,:right;velocity),))
        return (; model, times=collect(range(0,0.2;length=201)), stops=[0.1])
    elseif name == "custom2"
        phase = gas("h2o2")
        mixture = Dict("H2"=>1.0,"O2"=>1.0,"N2"=>3.76)
        reactor = WellStirredReactor(phase;temperature=920.0,pressure=one_atm,mass_fractions=mixture)
        environment = Reservoir(phase;temperature=920.0,pressure=one_atm,mass_fractions=mixture)
        # Source acceleration is 0.01*(P_reactor-P_environment); A/m = 0.01 m²/kg.
        model = MovingWallNetwork(ReactorNetwork((reactor=reactor,environment=environment));
            walls=(InertialWall(:reactor,:environment;mass=100.0),))
        return (; model, times=collect(range(0,0.5;length=101)), stops=Float64[])
    end
    throw(ArgumentError("unknown moving-wall example $name"))
end

function run_moving_wall_example(name, mechanism_directory, output; integrator,
                                reltol=name=="reactor2" ? 1e-12 : 1e-10,
                                abstol=name=="reactor2" ? 1e-21 : 1e-18)
    case = moving_wall_example(name,mechanism_directory)
    # Species masses and joule-valued quadrature states have different units.
    # A trace-species tolerance on a newly activated work ledger can force a
    # restarted multistep solver below floating-point time resolution.
    absolute_tolerances = fill(abstol,length(moving_wall_state(case.model)))
    ledger = case.model.ledger_offset
    absolute_tolerances[ledger] = max(abstol,1e-15)
    absolute_tolerances[ledger+1:ledger+2] .= max(abstol,1e-5)
    if name == "piston"
        # Restart with the released wall law at 0.1 s. A multistep solver's
        # history must not straddle the jump, including the final left-side RHS.
        held = MovingWallNetwork(case.model.network;walls=(MovingWall(:left,:right),))
        released = MovingWallNetwork(case.model.network;walls=(MovingWall(:left,:right;K=1e-4),))
        first_part = solve_moving_wall(held,(0.0,0.1);integrator,reltol,abstol=absolute_tolerances,
            saveat=filter(t->t<=0.1,case.times),dt=1e-10)
        second_part = solve_moving_wall(released,(0.1,last(case.times));integrator,reltol,abstol=absolute_tolerances,
            initial_state=first_part.u[end],saveat=filter(t->t>=0.1,case.times),dt=1e-10)
        solution = (;t=vcat(first_part.t,second_part.t[2:end]),
                     u=vcat(first_part.u,second_part.u[2:end]))
    else
        solution = solve_moving_wall(case.model,(first(case.times),last(case.times));integrator,
            saveat=case.times,reltol,abstol=absolute_tolerances,dt=1e-10)
    end
    rhs = moving_wall_rhs(case.model)
    vessel_names = [key for (key,node) in pairs(case.model.network.nodes) if node isa WellStirredReactor]
    open(output,"w") do io
        header = ["time_s"]
        for key in vessel_names
            append!(header,["$(key)_T_K","$(key)_P_Pa","$(key)_V_m3"])
        end
        append!(header,["gas_internal_energy_J","wall_kinetic_energy_J","energy_input_J","pressure_work_output_J"])
        println(io,join(header,','))
        for (t,u) in zip(solution.t,solution.u)
            d = moving_wall_diagnostics(rhs,u,t)
            row = Float64[t]
            for key in vessel_names
                state = getproperty(d.nodes,key)
                append!(row,[state.temperature,state.pressure,state.volume])
            end
            append!(row,[d.total_internal_energy,d.wall_kinetic_energy,d.energy_input,d.pressure_work_output])
            println(io,join(row,','))
        end
    end
    return (; model=case.model,solution)
end
