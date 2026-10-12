"A constant-temperature, constant-pressure ideal-gas boundary for a reactor network."
struct Reservoir{R<:IdealGasReactor}
    initial::R
end
Reservoir(gas::Solution; kwargs...) = Reservoir(IdealGasReactor(gas; kwargs...))

"""
    WellStirredReactor(initial; volume=1.0, chemistry=true)
    WellStirredReactor(gas; volume=1.0, chemistry=true, kwargs...)

A homogeneous fixed-volume vessel, with volume in m³. Its state is the species
masses in kg followed by temperature in K. Total mass and mass fractions are
derived from species masses. The initial `IdealGasReactor` must use
`:constant_volume`; the gas constructor selects that constraint automatically.
Energy may be adiabatic or isothermal. Set `chemistry=false` for pure mixing.
"""
struct WellStirredReactor{R<:IdealGasReactor}
    initial::R
    volume::Float64
    chemistry::Bool
end
function WellStirredReactor(initial::IdealGasReactor; volume=1.0, chemistry=true)
    initial.constraint === :constant_volume ||
        throw(ArgumentError("network vessels require constant-volume initial states"))
    isfinite(volume) && volume > 0 || throw(ArgumentError("positive finite volume required"))
    return WellStirredReactor(initial, Float64(volume), Bool(chemistry))
end
function WellStirredReactor(gas::Solution; volume=1.0, chemistry=true, kwargs...)
    initial = IdealGasReactor(gas; constraint=:constant_volume, kwargs...)
    return WellStirredReactor(initial; volume, chemistry)
end

abstract type NetworkFlowDevice end

"""
    MassFlowController(upstream, downstream; mdot)

One-way flow in kg/s, independent of pressure. Node identifiers are symbols.
`mdot` is a number, `f(t)`, or `f(states, t)`, where named node states expose
`mass`, `temperature`, `pressure`, `density`, and `mass_fractions`. For example,
`mdot=(states,t)->states.combustor.mass/0.1` sets a 0.1 s residence time.
Negative requested rates are clipped to zero. Treat callback states as read-only.
"""
struct MassFlowController{F} <: NetworkFlowDevice
    upstream::Symbol
    downstream::Symbol
    mdot::F
end
MassFlowController(upstream::Symbol, downstream::Symbol; mdot) =
    MassFlowController(upstream, downstream, mdot)

"A one-way valve: mdot = max(0, K * time_function(t) * pressure_function(ΔP))."
struct Valve{F,G} <: NetworkFlowDevice
    upstream::Symbol
    downstream::Symbol
    K::Float64
    pressure_function::F
    time_function::G
end
function Valve(upstream::Symbol, downstream::Symbol; K, pressure_function=identity,
               time_function=t -> 1.0)
    isfinite(K) && K >= 0 || throw(ArgumentError("valve K must be finite and nonnegative"))
    return Valve(upstream, downstream, Float64(K), pressure_function, time_function)
end

"""
    PressureController(upstream, downstream; primary, K, pressure_function=identity)

One-way flow equal to the primary device's flow plus `K*pressure_function(ΔP)`.
`primary` is an earlier device in the network's flow tuple, or its one-based
index. `K` has units kg/s/Pa for the default pressure function.
"""
struct PressureController{D,F} <: NetworkFlowDevice
    upstream::Symbol
    downstream::Symbol
    primary::D
    K::Float64
    pressure_function::F
end
function PressureController(upstream::Symbol, downstream::Symbol; primary, K,
                            pressure_function=identity)
    isfinite(K) && K >= 0 || throw(ArgumentError("controller K must be finite and nonnegative"))
    return PressureController(upstream, downstream, primary, Float64(K), pressure_function)
end

"""
    HeatTransferWall(left, right; area=1.0, U=0.0, heat_flux=0.0)

A fixed wall transferring `area*(U*(T_left-T_right)+heat_flux(t))` watts from
left to right. Area is in m², U in W/m²/K, and prescribed heat flux in W/m².
`heat_flux` may be a constant or a time function. No wall motion is modeled.
"""
struct HeatTransferWall{F}
    left::Symbol
    right::Symbol
    area::Float64
    U::Float64
    heat_flux::F
end
function HeatTransferWall(left::Symbol, right::Symbol; area=1.0, U=0.0, heat_flux=0.0)
    all(x -> isfinite(x) && x >= 0, (area, U)) ||
        throw(ArgumentError("wall area and U must be finite and nonnegative"))
    return HeatTransferWall(left, right, Float64(area), Float64(U), heat_flux)
end

"""
    ReactorNetwork(nodes::NamedTuple; flows=(), walls=())

Connect named `WellStirredReactor`s and `Reservoir`s with one-way flow devices
and fixed heat-transfer walls. Species are mapped by name. Every species that
can enter a vessel must exist in its mechanism; incompatible molecular weights
are rejected. Reservoir outlets need not contain every exhaust species.

All calculations use ideal-gas thermodynamics. This network supports fixed
volumes, homogeneous chemistry, mixing, pressure regulation, and heat transfer.
Surface chemistry, moving walls, and non-ideal phases are excluded.
"""
struct ReactorNetwork{N,F,W}
    nodes::N
    flows::F
    walls::W
    offsets::Vector{Int}
    flow_endpoints::Vector{Tuple{Int,Int}}
    wall_endpoints::Vector{Tuple{Int,Int}}
    species_maps::Vector{Vector{Int}}
    primary_indices::Vector{Int}
    initial_state::Vector{Float64}
end

function _network_endpoint(names, name)
    index = findfirst(==(name), names)
    isnothing(index) && throw(ArgumentError("unknown network node: $name"))
    return index
end

function ReactorNetwork(nodes::NamedTuple; flows=(), walls=())
    all(node -> node isa Union{Reservoir,WellStirredReactor}, values(nodes)) ||
        throw(ArgumentError("nodes must be Reservoirs or WellStirredReactors"))
    any(node -> node isa WellStirredReactor, values(nodes)) ||
        throw(ArgumentError("a network must contain a reactor"))
    flows, walls = Tuple(flows), Tuple(walls)
    all(device -> device isa NetworkFlowDevice, flows) || throw(ArgumentError("invalid flow device"))
    all(wall -> wall isa HeatTransferWall, walls) || throw(ArgumentError("invalid heat-transfer wall"))
    offsets = zeros(Int, length(nodes))
    u0 = Float64[]
    for (i, node) in enumerate(nodes)
        if node isa WellStirredReactor
            offsets[i] = length(u0) + 1
            mass = node.initial.density * node.volume
            append!(u0, mass .* node.initial.mass_fractions)
            push!(u0, node.initial.temperature)
        end
    end
    flow_endpoints, wall_endpoints = Tuple{Int,Int}[], Tuple{Int,Int}[]
    maps, primary = Vector{Int}[], zeros(Int, length(flows))
    for (i, device) in enumerate(flows)
        up = _network_endpoint(keys(nodes), device.upstream)
        down = _network_endpoint(keys(nodes), device.downstream)
        up == down && throw(ArgumentError("a flow device must connect distinct nodes"))
        push!(flow_endpoints, (up, down))
        source, target = nodes[up].initial, nodes[down].initial
        mapping = zeros(Int, source.gas.n_species)
        if nodes[down] isa WellStirredReactor
            for k in eachindex(mapping)
                j = findfirst(==(source.gas.species_names[k]), target.gas.species_names)
                if isnothing(j)
                    (nodes[up] isa WellStirredReactor || source.mass_fractions[k] > 0) &&
                        throw(ArgumentError("downstream mechanism lacks $(source.gas.species_names[k])"))
                else
                    isapprox(source.gas.MW[k], target.gas.MW[j]; rtol=1e-10) ||
                        throw(ArgumentError("incompatible molecular weight for $(source.gas.species_names[k])"))
                    for element in union(source.gas.elements, target.gas.elements)
                        a = findfirst(==(element), source.gas.elements)
                        b = findfirst(==(element), target.gas.elements)
                        source_atoms = isnothing(a) ? 0 : source.gas.ele_matrix[a,k]
                        target_atoms = isnothing(b) ? 0 : target.gas.ele_matrix[b,j]
                        source_atoms == target_atoms ||
                            throw(ArgumentError("incompatible elemental composition for $(source.gas.species_names[k])"))
                    end
                    mapping[k] = j
                end
            end
        end
        push!(maps, mapping)
        if device isa PressureController
            candidates = device.primary isa Integer ? [Int(device.primary)] :
                         findall(flow -> flow === device.primary, flows)
            length(candidates) == 1 && 1 <= only(candidates) < i ||
                throw(ArgumentError("pressure controller primary must identify one earlier flow device"))
            primary[i] = only(candidates)
        end
    end
    for wall in walls
        left = _network_endpoint(keys(nodes), wall.left)
        right = _network_endpoint(keys(nodes), wall.right)
        left == right && throw(ArgumentError("a wall must connect distinct nodes"))
        push!(wall_endpoints, (left, right))
    end
    return ReactorNetwork(nodes, flows, walls, offsets, flow_endpoints, wall_endpoints,
                          maps, primary, u0)
end

"Return independent species-mass/temperature initial states in named-node order."
network_state(network::ReactorNetwork) = copy(network.initial_state)

mutable struct NetworkNodeState
    mass::Float64
    temperature::Float64
    pressure::Float64
    density::Float64
    volume::Float64
    mass_fractions::Vector{Float64}
    cp::Float64
    cv::Float64
    enthalpy::Float64
    internal_energy::Float64
end

struct NetworkRHS{N,S}
    network::N
    states::S
    state_vector::Vector{NetworkNodeState}
    workspaces::Vector{ReactorWorkspace{Float64}}
    mass_flow_rates::Vector{Float64}
    wall_heat_rates::Vector{Float64}
    energy_flux::Vector{Float64}
    thermostat_power::Vector{Float64}
    jac_state::Vector{Float64}
    base::Vector{Float64}
    plus::Vector{Float64}
    minus::Vector{Float64}
end

"Prepare an in-place native network RHS. Each concurrent solve needs its own instance."
function network_rhs(network::ReactorNetwork)
    states = map(network.nodes) do node
        initial = node.initial
        volume = node isa WellStirredReactor ? node.volume : Inf
        NetworkNodeState(initial.density * volume, initial.temperature, initial.pressure,
            initial.density, volume, Float64.(initial.mass_fractions), 0.0, 0.0, 0.0, 0.0)
    end
    workspaces = [ReactorWorkspace(node.initial.gas, Float64) for node in network.nodes]
    return NetworkRHS(network, states, collect(values(states)), workspaces, zeros(length(network.flows)),
        zeros(length(network.walls)), zeros(length(network.nodes)), zeros(length(network.nodes)),
        network_state(network), zeros(length(network.initial_state)),
        zeros(length(network.initial_state)), zeros(length(network.initial_state)))
end

_network_signal(value::Real, states, t) = value
_network_signal(f, states, t) = applicable(f, states, t) ? f(states, t) : f(t)
_network_time_signal(value::Real, t) = value
_network_time_signal(f, t) = f(t)

# Keep each heterogeneous node/device concrete. Runtime tuple indexing boxes
# the complete immutable mechanism-bearing node tuple on every iteration.
@inline _network_foreach(f, ::Tuple{}, index::Int=1) = nothing
@inline function _network_foreach(f::F, items::Tuple, index::Int=1) where F
    f(first(items), index)
    _network_foreach(f, Base.tail(items), index + 1)
    return nothing
end

function _update_network_states!(rhs, u; volumes=nothing)
    length(u) == length(rhs.network.initial_state) || throw(DimensionMismatch("invalid network state length"))
    volumes === nothing || length(volumes) == length(rhs.network.nodes) ||
        throw(DimensionMismatch("one volume per network node required"))
    _network_foreach(values(rhs.network.nodes)) do node, i
        state, workspace = rhs.state_vector[i], rhs.workspaces[i]
        gas = node.initial.gas
        offset, ns = rhs.network.offsets[i], gas.n_species
        if offset > 0
            mass = sum(@view u[offset:offset+ns-1])
            temperature = node.initial.energy === :isothermal ? node.initial.temperature : u[offset+ns]
            isfinite(mass) && mass > 0 && isfinite(temperature) && temperature > 0 ||
                throw(DomainError((mass, temperature), "positive finite mass and temperature required"))
            volume = volumes === nothing ? node.volume : volumes[i]
            isfinite(volume) && volume > 0 || throw(DomainError(volume, "positive finite vessel volume required"))
            state.volume = volume
            state.mass, state.temperature, state.density = mass, temperature, mass / volume
            @inbounds for k in 1:ns
                state.mass_fractions[k] = u[offset+k-1] / mass
            end
        end
        inverse_mw = 0.0
        @inbounds for k in 1:ns
            inverse_mw += state.mass_fractions[k] / gas.MW[k]
        end
        state.pressure = state.density * R * state.temperature * inverse_mw
        @inbounds for k in 1:ns
            workspace.X[k] = state.mass_fractions[k] / (gas.MW[k] * inverse_mw)
            workspace.C[k] = max(state.mass_fractions[k], 0.0) * state.density / gas.MW[k]
        end
        cal_cp_R!(workspace.cp_R, gas, state.temperature, state.pressure, workspace.X)
        cal_h_RT!(workspace.h_mole, gas, state.temperature, state.pressure, workspace.X)
        cal_s0_R!(workspace.entropy, gas, state.temperature, state.pressure, workspace.X)
        cp, enthalpy = 0.0, 0.0
        @inbounds for k in 1:ns
            workspace.h_mole[k] *= R * state.temperature
            workspace.entropy[k] *= R
            cp += state.mass_fractions[k] * workspace.cp_R[k] * R / gas.MW[k]
            enthalpy += state.mass_fractions[k] * workspace.h_mole[k] / gas.MW[k]
        end
        state.cp, state.cv = cp, cp - R * inverse_mw
        state.enthalpy, state.internal_energy = enthalpy, enthalpy - R * state.temperature * inverse_mw
        state.cv > 0 || throw(DomainError(state.cv, "positive heat capacity required"))
    end
    return nothing
end

function (rhs::NetworkRHS)(du, u, p, t; volumes=nothing)
    network = rhs.network
    length(du) == length(u) || throw(DimensionMismatch("network derivative length mismatch"))
    fill!(du, 0)
    fill!(rhs.energy_flux, 0)
    fill!(rhs.thermostat_power, 0)
    _update_network_states!(rhs, u; volumes)
    _network_foreach(values(network.nodes)) do node, i
        offset = network.offsets[i]
        offset == 0 && return nothing
        gas, workspace, state = node.initial.gas, rhs.workspaces[i], rhs.state_vector[i]
        if node.chemistry
            wdot!(workspace.wdot, gas.reaction, state.temperature, workspace.C,
                workspace.entropy, workspace.h_mole, workspace.kinetics;
                rate_multipliers=node.initial.rate_multipliers)
            @inbounds for k in 1:gas.n_species
                du[offset+k-1] = workspace.wdot[k] * gas.MW[k] * state.volume
            end
        end
    end
    _network_foreach(network.flows) do device, i
        up, down = network.flow_endpoints[i]
        source, target = rhs.state_vector[up], rhs.state_vector[down]
        delta_p = source.pressure - target.pressure
        rate = if device isa MassFlowController
            _network_signal(device.mdot, rhs.states, t)
        elseif device isa Valve
            device.K * device.time_function(t) * device.pressure_function(delta_p)
        else
            rhs.mass_flow_rates[network.primary_indices[i]] + device.K * device.pressure_function(delta_p)
        end
        isfinite(rate) || throw(DomainError(rate, "finite mass flow rate required"))
        rate = max(rate, 0.0)
        rhs.mass_flow_rates[i] = rate
        upstream_offset, downstream_offset = network.offsets[up], network.offsets[down]
        if upstream_offset > 0
            @inbounds for k in eachindex(source.mass_fractions)
                du[upstream_offset+k-1] -= rate * source.mass_fractions[k]
            end
            rhs.energy_flux[up] -= rate * source.enthalpy
        end
        if downstream_offset > 0
            @inbounds for k in eachindex(source.mass_fractions)
                j = network.species_maps[i][k]
                j > 0 && (du[downstream_offset+j-1] += rate * source.mass_fractions[k])
            end
            rhs.energy_flux[down] += rate * source.enthalpy
        end
    end
    _network_foreach(network.walls) do wall, i
        left, right = network.wall_endpoints[i]
        power = wall.area * (wall.U * (rhs.state_vector[left].temperature - rhs.state_vector[right].temperature) +
                            _network_time_signal(wall.heat_flux, t))
        isfinite(power) || throw(DomainError(power, "finite wall heat rate required"))
        rhs.wall_heat_rates[i] = power
        network.offsets[left] > 0 && (rhs.energy_flux[left] -= power)
        network.offsets[right] > 0 && (rhs.energy_flux[right] += power)
    end
    _network_foreach(values(network.nodes)) do node, i
        offset = network.offsets[i]
        offset == 0 && return nothing
        gas, state, workspace = node.initial.gas, rhs.state_vector[i], rhs.workspaces[i]
        composition_energy = 0.0
        @inbounds for k in 1:gas.n_species
            composition_energy += (workspace.h_mole[k] - R * state.temperature) / gas.MW[k] * du[offset+k-1]
        end
        if node.initial.energy === :isothermal
            du[offset+gas.n_species] = 0
            rhs.thermostat_power[i] = composition_energy - rhs.energy_flux[i]
        else
            du[offset+gas.n_species] = (rhs.energy_flux[i] - composition_energy) / (state.mass * state.cv)
        end
    end
    return nothing
end

"Dense second-order finite-difference Jacobian, including all device/state couplings."
function network_jacobian!(J, u, rhs::NetworkRHS, t=0.0)
    n = length(u)
    size(J) == (n, n) || throw(DimensionMismatch("invalid network Jacobian dimensions"))
    copyto!(rhs.jac_state, u)
    rhs(rhs.base, u, nothing, t)
    relative_step = cbrt(eps(Float64))
    _network_foreach(values(rhs.network.nodes)) do node, node_index
        offset = rhs.network.offsets[node_index]
        offset == 0 && return nothing
        mass_scale = rhs.state_vector[node_index].mass * 1e-6
        for j in offset:offset+node.initial.gas.n_species
            step = relative_step * max(abs(u[j]), j == offset+node.initial.gas.n_species ? 1.0 : mass_scale)
            rhs.jac_state[j] = u[j] + step
            step = rhs.jac_state[j] - u[j]
            rhs(rhs.plus, rhs.jac_state, nothing, t)
            if u[j] >= step
                rhs.jac_state[j] = u[j] - step
                rhs(rhs.minus, rhs.jac_state, nothing, t)
                @inbounds for i in 1:n
                    J[i, j] = (rhs.plus[i] - rhs.minus[i]) / (2step)
                end
            else
                rhs.jac_state[j] = u[j] + 2step
                rhs(rhs.minus, rhs.jac_state, nothing, t)
                @inbounds for i in 1:n
                    J[i, j] = (-3rhs.base[i] + 4rhs.plus[i] - rhs.minus[i]) / (2step)
                end
            end
            rhs.jac_state[j] = u[j]
        end
    end
    return nothing
end

function _network_tgrad!(du, u, rhs, t)
    step = cbrt(eps(Float64)) * max(abs(t), 1e-6)
    rhs(rhs.base, u, nothing, t)
    rhs(rhs.plus, u, nothing, t + step)
    rhs(rhs.minus, u, nothing, t + 2step)
    @. du = (-3rhs.base + 4rhs.plus - rhs.minus) / (2step)
    return nothing
end

"Check positive temperature/mass and reject species masses below a relative tolerance."
function network_isoutofdomain(network::ReactorNetwork, u; tolerance=1e-13)
    length(u) == length(network.initial_state) && all(isfinite, u) || return true
    for (i, node) in enumerate(network.nodes)
        offset = network.offsets[i]
        offset == 0 && continue
        masses = @view u[offset:offset+node.initial.gas.n_species-1]
        mass = sum(masses)
        mass > 0 && u[offset+node.initial.gas.n_species] > 0 || return true
        minimum(masses) < -tolerance * mass && return true
    end
    return false
end

"Prepare native callbacks and state for a caller-owned stiff ODE solver."
function network_problem(network::ReactorNetwork, tspan; initial_state=network_state(network))
    length(tspan) == 2 && all(isfinite, tspan) && tspan[2] > tspan[1] ||
        throw(ArgumentError("finite increasing tspan required"))
    network_isoutofdomain(network, initial_state) && throw(ArgumentError("invalid initial network state"))
    rhs = network_rhs(network)
    return (f=rhs, jac=(J,u,p,t) -> network_jacobian!(J,u,rhs,t),
            tgrad=(du,u,p,t) -> _network_tgrad!(du,u,rhs,t), u0=Float64.(initial_state),
            tspan=(float(tspan[1]), float(tspan[2])), p=nothing,
            isoutofdomain=(u,p,t) -> network_isoutofdomain(network,u))
end

"Integrate with integrator(problem; kwargs...), using the same contract as solve_reactor."
function solve_network(network::ReactorNetwork, tspan; integrator,
                       initial_state=network_state(network), kwargs...)
    return integrator(network_problem(network, tspan; initial_state); kwargs...)
end

"""
    network_diagnostics(network_or_rhs, u, t=0.0)

Return independent named node properties, device flow/heat rates, total and
external mass/energy rates, and net elemental rates. Energy accounting includes
reservoir enthalpy flow, reservoir wall heat, and imposed isothermal bath heat.
Rates are kg/s, W, and kmol of element/s. Reuse a prepared RHS to avoid repeatedly
allocating workspaces when evaluating a trajectory.
"""
function network_diagnostics(rhs::NetworkRHS, u, t=0.0)
    derivative = similar(u)
    rhs(derivative, u, nothing, t)
    network = rhs.network
    total_mass_rate, total_energy_rate, external_mass_rate = 0.0, 0.0, 0.0
    element_rates = Dict{String,Float64}()
    for (i, node) in enumerate(network.nodes)
        offset = network.offsets[i]
        offset == 0 && continue
        gas, state, workspace = node.initial.gas, rhs.state_vector[i], rhs.workspaces[i]
        dm = @view derivative[offset:offset+gas.n_species-1]
        total_mass_rate += sum(dm)
        total_energy_rate += state.mass * state.cv * derivative[offset+gas.n_species]
        @inbounds for k in 1:gas.n_species
            total_energy_rate += (workspace.h_mole[k] - R*state.temperature) * dm[k] / gas.MW[k]
        end
        elements = gas.ele_matrix * (dm ./ gas.MW)
        for (name, value) in zip(gas.elements, elements)
            element_rates[name] = get(element_rates, name, 0.0) + value
        end
    end
    for (i, (up, down)) in enumerate(network.flow_endpoints)
        if network.offsets[up] == 0 && network.offsets[down] > 0
            external_mass_rate += rhs.mass_flow_rates[i]
        elseif network.offsets[up] > 0 && network.offsets[down] == 0
            external_mass_rate -= rhs.mass_flow_rates[i]
        end
    end
    properties = map(rhs.states) do state
        (; mass=state.mass, temperature=state.temperature, pressure=state.pressure,
           density=state.density, volume=state.volume, mass_fractions=copy(state.mass_fractions),
           enthalpy=state.enthalpy, internal_energy=state.internal_energy,
           total_internal_energy=state.mass * state.internal_energy)
    end
    return (; nodes=properties, mass_flow_rates=copy(rhs.mass_flow_rates),
        wall_heat_rates=copy(rhs.wall_heat_rates), thermostat_heat_rates=copy(rhs.thermostat_power),
        total_mass_rate, external_mass_rate, total_energy_rate,
        external_energy_rate=sum(rhs.energy_flux) + sum(rhs.thermostat_power), element_rates)
end
network_diagnostics(network::ReactorNetwork, u=network_state(network), t=0.0) =
    network_diagnostics(network_rhs(network), u, t)

"""
    solve_network_steady(network; integrator, interval=1, max_time=100,
                         steady_tolerance=1e-8, initial_state=network_state(network), kwargs...)

Relax an autonomous network in finite integration intervals until species-mass
rates divided by vessel mass and temperature rates divided by temperature are
below `steady_tolerance` (s⁻¹). The caller's solution must expose its states as
`solution.u`. Returns `(state, time, residual)`; throws if convergence is not
reached by `max_time`. Time-varying boundary conditions require transient solves.
"""
function solve_network_steady(network::ReactorNetwork; integrator, interval=1.0,
        max_time=100.0, steady_tolerance=1e-8, initial_state=network_state(network), kwargs...)
    all(x -> isfinite(x) && x > 0, (interval, max_time, steady_tolerance)) ||
        throw(ArgumentError("positive finite steady integration settings required"))
    state, time, residual = Float64.(initial_state), 0.0, Inf
    rhs, derivative = network_rhs(network), similar(state)
    while time < max_time
        stop = min(time + interval, max_time)
        solution = solve_network(network, (time, stop); integrator, initial_state=state, kwargs...)
        state, time = copy(solution.u[end]), stop
        rhs(derivative, state, nothing, time)
        residual = 0.0
        for (i, node) in enumerate(network.nodes)
            offset = network.offsets[i]
            offset == 0 && continue
            ns = node.initial.gas.n_species
            residual = max(residual, maximum(abs, @view derivative[offset:offset+ns-1]) / rhs.state_vector[i].mass,
                           abs(derivative[offset+ns]) / rhs.state_vector[i].temperature)
        end
        residual <= steady_tolerance && return (; state, time, residual)
    end
    error("network did not reach steady state by $max_time s; residual=$residual s^-1")
end

export Reservoir, WellStirredReactor, ReactorNetwork
export MassFlowController, Valve, PressureController, HeatTransferWall
export network_state, network_rhs, network_jacobian!, network_problem, network_isoutofdomain
export network_diagnostics, solve_network, solve_network_steady
