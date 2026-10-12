abstract type ReactorMovingWall end

"""
    MovingWall(left, right; area=1, K=0, velocity=0, U=0, heat_flux=0)

A massless piston between named network nodes. Positive velocity expands the
left vessel and compresses the right. Its velocity is
`K*(P_left-P_right) + velocity(states,t)`, in m/s; prescribed velocity may also
be a number or `f(t)`. Heat flowing left to right is
`area*(U*(T_left-T_right) + heat_flux(t))`, in W.

Pressure work uses each vessel's own pressure. Consequently a pressure-driven
massless piston removes mechanical energy from the gases at
`area*velocity*(P_left-P_right)`; this work is recorded by the network's ledger.
The piston has neither inertia nor a stored thermal energy.
"""
struct MovingWall{F,Q} <: ReactorMovingWall
    left::Symbol
    right::Symbol
    area::Float64
    K::Float64
    velocity::F
    U::Float64
    heat_flux::Q
end
function MovingWall(left::Symbol, right::Symbol; area=1.0, K=0.0, velocity=0.0,
                    U=0.0, heat_flux=0.0)
    isfinite(area) && area > 0 && all(x -> isfinite(x) && x >= 0, (K, U)) ||
        throw(ArgumentError("positive wall area and nonnegative finite K/U required"))
    return MovingWall(left, right, Float64(area), Float64(K), velocity, Float64(U), heat_flux)
end

"""
    InertialWall(left, right; area=1, mass, initial_velocity=0, U=0, heat_flux=0)

A piston with one integrated velocity and equation `mass*dv/dt = area*(P_left-P_right)`.
Wall mass is in kg and positive velocity expands the left vessel. Gas pressure
work transfers energy to the piston's kinetic energy. Heat transfer follows
`MovingWall`. No friction, wall heat capacity, or contact stops are modeled.
"""
struct InertialWall{Q} <: ReactorMovingWall
    left::Symbol
    right::Symbol
    area::Float64
    mass::Float64
    initial_velocity::Float64
    U::Float64
    heat_flux::Q
end
function InertialWall(left::Symbol, right::Symbol; area=1.0, mass,
                      initial_velocity=0.0, U=0.0, heat_flux=0.0)
    all(x -> isfinite(x) && x > 0, (area, mass)) && isfinite(initial_velocity) &&
        isfinite(U) && U >= 0 || throw(ArgumentError("invalid inertial-wall parameters"))
    return InertialWall(left, right, Float64(area), Float64(mass),
                        Float64(initial_velocity), Float64(U), heat_flux)
end

"""
    MovingWallNetwork(network::ReactorNetwork; walls)

Add moving walls to an ideal-gas network, retaining its chemistry, flow devices,
and fixed heat-transfer walls. Reservoirs remain fixed boundaries. State order:
the base network's species masses/temperatures, vessel volumes, inertial-wall
velocities, then cumulative external mass input, energy input (including
isothermal thermostats), and gas pressure work output. Use `moving_wall_state`
and `moving_wall_diagnostics` instead of interpreting indices directly.

All runtime calculations are Julia. Callers supply a stiff ODE integrator to
`solve_moving_wall`. At a prescribed discontinuity, restart with each segment's
continuous wall law; `example/reactors/piston.jl` demonstrates this for a release.
"""
struct MovingWallNetwork{N,W}
    network::N
    walls::W
    wall_endpoints::Vector{Tuple{Int,Int}}
    volume_indices::Vector{Int}
    velocity_indices::Vector{Int}
    ledger_offset::Int
    initial_state::Vector{Float64}
end
function MovingWallNetwork(network::ReactorNetwork; walls)
    wall_tuple = Tuple(walls)
    all(w -> w isa ReactorMovingWall, wall_tuple) || throw(ArgumentError("moving walls required"))
    names = keys(network.nodes)
    endpoints = Tuple{Int,Int}[]
    for wall in wall_tuple
        left, right = _network_endpoint(names, wall.left), _network_endpoint(names, wall.right)
        left != right || throw(ArgumentError("a wall needs two distinct nodes"))
        network.offsets[left] > 0 || network.offsets[right] > 0 ||
            throw(ArgumentError("a wall must connect to a vessel"))
        push!(endpoints, (left, right))
    end
    initial = network_state(network)
    volume_indices = zeros(Int, length(network.nodes))
    for (i, node) in enumerate(network.nodes)
        network.offsets[i] == 0 && continue
        push!(initial, node.volume)
        volume_indices[i] = length(initial)
    end
    velocity_indices = zeros(Int, length(wall_tuple))
    for (i, wall) in enumerate(wall_tuple)
        wall isa InertialWall || continue
        push!(initial, wall.initial_velocity)
        velocity_indices[i] = length(initial)
    end
    ledger_offset = length(initial) + 1
    append!(initial, zeros(3))
    return MovingWallNetwork(network, wall_tuple, endpoints, volume_indices,
                             velocity_indices, ledger_offset, initial)
end

"Independent initial species-mass, temperature, volume, velocity, and ledger state."
moving_wall_state(model::MovingWallNetwork) = copy(model.initial_state)

struct MovingWallRHS{M,B}
    model::M
    network_rhs::B
    volumes::Vector{Float64}
    volume_rates::Vector{Float64}
    wall_velocities::Vector{Float64}
    wall_heat_rates::Vector{Float64}
    motion_heat::Vector{Float64}
    jac_state::Vector{Float64}
    base::Vector{Float64}
    plus::Vector{Float64}
    minus::Vector{Float64}
end
function moving_wall_rhs(model::MovingWallNetwork)
    n, nn, nw = length(model.initial_state), length(model.network.nodes), length(model.walls)
    return MovingWallRHS(model, network_rhs(model.network), fill(Inf, nn), zeros(nn),
        zeros(nw), zeros(nw), zeros(nn), moving_wall_state(model), zeros(n), zeros(n), zeros(n))
end

function (rhs::MovingWallRHS)(du, u, p, t)
    model, base = rhs.model, rhs.network_rhs
    network = model.network
    length(du) == length(u) == length(model.initial_state) ||
        throw(DimensionMismatch("invalid moving-wall state length"))
    fill!(du, 0)
    fill!(rhs.volume_rates, 0)
    fill!(rhs.motion_heat, 0)
    for (i, index) in enumerate(model.volume_indices)
        index > 0 && (rhs.volumes[i] = u[index])
    end
    nbase = length(network.initial_state)
    base(@view(du[1:nbase]), @view(u[1:nbase]), p, t; volumes=rhs.volumes)
    _network_foreach(model.walls) do wall, i
        left, right = model.wall_endpoints[i]
        l, r = base.state_vector[left], base.state_vector[right]
        if wall isa InertialWall
            velocity_index = model.velocity_indices[i]
            velocity = u[velocity_index]
            du[velocity_index] = wall.area * (l.pressure - r.pressure) / wall.mass
        else
            velocity = wall.K * (l.pressure - r.pressure) + _network_signal(wall.velocity, base.states, t)
        end
        heat = wall.area * (wall.U * (l.temperature - r.temperature) + _network_time_signal(wall.heat_flux, t))
        isfinite(velocity) && isfinite(heat) || throw(DomainError((velocity, heat), "finite wall rates required"))
        rhs.wall_velocities[i], rhs.wall_heat_rates[i] = velocity, heat
        if model.volume_indices[left] > 0
            rhs.volume_rates[left] += wall.area * velocity
            rhs.motion_heat[left] -= heat
        end
        if model.volume_indices[right] > 0
            rhs.volume_rates[right] -= wall.area * velocity
            rhs.motion_heat[right] += heat
        end
    end
    index = model.ledger_offset
    _network_foreach(values(network.nodes)) do node, i
        offset = network.offsets[i]
        offset == 0 && return nothing
        state = base.state_vector[i]
        volume_rate = rhs.volume_rates[i]
        work = state.pressure * volume_rate
        motion_energy = rhs.motion_heat[i] - work
        du[model.volume_indices[i]] = volume_rate
        if node.initial.energy === :isothermal
            base.thermostat_power[i] -= motion_energy
        else
            du[offset+node.initial.gas.n_species] += motion_energy / (state.mass * state.cv)
        end
        du[index+1] += base.energy_flux[i] + rhs.motion_heat[i] + base.thermostat_power[i]
        du[index+2] += work
    end
    mass_input = 0.0
    for (i, (up, down)) in enumerate(network.flow_endpoints)
        network.offsets[up] == 0 && network.offsets[down] > 0 && (mass_input += base.mass_flow_rates[i])
        network.offsets[up] > 0 && network.offsets[down] == 0 && (mass_input -= base.mass_flow_rates[i])
    end
    du[index] = mass_input
    return nothing
end

"Reject nonpositive volume, temperature, and mass, or materially negative species mass."
function moving_wall_isoutofdomain(model::MovingWallNetwork, u; tolerance=1e-13)
    length(u) == length(model.initial_state) && all(isfinite, u) || return true
    network_isoutofdomain(model.network, @view(u[1:length(model.network.initial_state)]); tolerance) && return true
    return any(index -> index > 0 && u[index] <= 0, model.volume_indices)
end

"Dense finite-difference Jacobian with species, temperature, volume, and velocity couplings."
function moving_wall_jacobian!(J, u, rhs::MovingWallRHS, t=0.0)
    n = length(u)
    size(J) == (n,n) || throw(DimensionMismatch("invalid moving-wall Jacobian dimensions"))
    fill!(J, 0)
    copyto!(rhs.jac_state, u)
    rhs(rhs.base, u, nothing, t)
    relative_step = cbrt(eps(Float64))
    for j in 1:rhs.model.ledger_offset-1
        scale = j in rhs.model.velocity_indices ? 1.0 : 1e-6
        for i in eachindex(rhs.model.network.offsets)
            offset = rhs.model.network.offsets[i]
            offset == 0 && continue
            offset <= j < offset+length(rhs.network_rhs.state_vector[i].mass_fractions) &&
                (scale = rhs.network_rhs.state_vector[i].mass * 1e-6)
        end
        step = relative_step * max(abs(u[j]), scale)
        rhs.jac_state[j] = u[j] + step
        step = rhs.jac_state[j] - u[j]
        rhs(rhs.plus, rhs.jac_state, nothing, t)
        if u[j] >= step || j in rhs.model.velocity_indices
            rhs.jac_state[j] = u[j] - step
            rhs(rhs.minus, rhs.jac_state, nothing, t)
            @views @. J[:,j] = (rhs.plus - rhs.minus) / (2step)
        else
            rhs.jac_state[j] = u[j] + 2step
            rhs(rhs.minus, rhs.jac_state, nothing, t)
            @views @. J[:,j] = (-3rhs.base + 4rhs.plus - rhs.minus) / (2step)
        end
        rhs.jac_state[j] = u[j]
    end
    return nothing
end

"Prepare native callbacks for a caller-owned stiff ODE integrator."
function moving_wall_problem(model::MovingWallNetwork, tspan; initial_state=moving_wall_state(model))
    length(tspan) == 2 && all(isfinite, tspan) && tspan[2] > tspan[1] ||
        throw(ArgumentError("finite increasing tspan required"))
    moving_wall_isoutofdomain(model, initial_state) && throw(ArgumentError("invalid initial moving-wall state"))
    rhs = moving_wall_rhs(model)
    return (f=rhs, jac=(J,u,p,t) -> moving_wall_jacobian!(J,u,rhs,t),
        tgrad=(du,u,p,t) -> _network_tgrad!(du,u,rhs,t), u0=Float64.(initial_state),
        tspan=(float(tspan[1]),float(tspan[2])), p=nothing,
        isoutofdomain=(u,p,t) -> moving_wall_isoutofdomain(model,u))
end
function solve_moving_wall(model::MovingWallNetwork, tspan; integrator,
                          initial_state=moving_wall_state(model), kwargs...)
    return integrator(moving_wall_problem(model,tspan;initial_state); kwargs...)
end

"""
    moving_wall_diagnostics(model_or_rhs, u, t=0)

Report vessel states, wall velocity/heat/pressure work, and integrated balances.
`mass_balance = total_mass - external_mass_input` and
`energy_balance = total_internal_energy + pressure_work_output - energy_input`
are invariants. Heat and work are signed in J, masses in kg, and elemental
inventories in kmol. Internal gas heat transfers cancel from energy input.
`wall_kinetic_energy` accounts for inertial walls separately; massless piston
work need not vanish, even for an otherwise closed network.
"""
function moving_wall_diagnostics(rhs::MovingWallRHS, u, t=0.0)
    derivative = similar(u)
    rhs(derivative,u,nothing,t)
    model, base = rhs.model, rhs.network_rhs
    mass, energy, total_volume, energy_rate = 0.0, 0.0, 0.0, 0.0
    elements = Dict{String,Float64}()
    for (i, node) in enumerate(model.network.nodes)
        offset = model.network.offsets[i]
        offset == 0 && continue
        state, gas = base.state_vector[i], node.initial.gas
        mass += state.mass
        energy += state.mass * state.internal_energy
        total_volume += state.volume
        dm = @view derivative[offset:offset+gas.n_species-1]
        energy_rate += state.mass * state.cv * derivative[offset+gas.n_species]
        for k in 1:gas.n_species
            energy_rate += (base.workspaces[i].h_mole[k]-R*state.temperature) * dm[k]/gas.MW[k]
        end
        amounts = gas.ele_matrix * ((state.mass .* state.mass_fractions) ./ gas.MW)
        for (name, amount) in zip(gas.elements, amounts)
            elements[name] = get(elements, name, 0.0) + amount
        end
    end
    kinetic_energy = sum((wall isa InertialWall ? wall.mass*rhs.wall_velocities[i]^2/2 : 0.0)
                         for (i,wall) in enumerate(model.walls); init=0.0)
    nodes = map(base.states) do state
        (; mass=state.mass, temperature=state.temperature, pressure=state.pressure,
           volume=state.volume, density=state.density, mass_fractions=copy(state.mass_fractions),
           internal_energy=state.internal_energy, enthalpy=state.enthalpy)
    end
    index = model.ledger_offset
    return (; nodes, wall_velocities=copy(rhs.wall_velocities), wall_heat_rates=copy(rhs.wall_heat_rates),
        volume_rates=copy(rhs.volume_rates), mass_flow_rates=copy(base.mass_flow_rates),
        thermostat_heat_rates=copy(base.thermostat_power), total_mass=mass, total_volume,
        total_internal_energy=energy, wall_kinetic_energy=kinetic_energy, element_inventories=elements,
        external_mass_input=u[index], energy_input=u[index+1], pressure_work_output=u[index+2],
        external_mass_rate=derivative[index], energy_input_rate=derivative[index+1],
        pressure_work_rate=derivative[index+2], total_internal_energy_rate=energy_rate,
        mass_balance=mass-u[index], energy_balance=energy+u[index+2]-u[index+1])
end
moving_wall_diagnostics(model::MovingWallNetwork, u=moving_wall_state(model), t=0.0) =
    moving_wall_diagnostics(moving_wall_rhs(model),u,t)

export MovingWall, InertialWall, MovingWallNetwork, moving_wall_state, moving_wall_rhs
export moving_wall_jacobian!, moving_wall_isoutofdomain, moving_wall_problem, solve_moving_wall
export moving_wall_diagnostics
