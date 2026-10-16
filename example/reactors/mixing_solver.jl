# Full reacting stationary mixer with damped Newton and pseudo-transient steps.
using Arrhenius, LinearAlgebra

include(joinpath(@__DIR__, "network_cases.jl"))

"One accepted Newton/pseudo-transient iterate (iteration 0 records the initial state)."
struct MixingIteration
    iteration::Int
    step::Float64
    shift_per_s::Float64
    residual::Float64
    physical_residual::Float64
    backtracks::Int
end

function _mixing_fixed_scales(network, u)
    scale = similar(u)
    for (i, node) in enumerate(network.nodes)
        first = network.offsets[i]
        first == 0 && continue
        ns = node.initial.gas.n_species
        scale[first:first+ns-1] .= sum(@view u[first:first+ns-1])
        scale[first+ns] = u[first+ns]
    end
    scale
end

"Fixed-scale residual: max |f_i| / scale_i with scales frozen at the initial state."
function _mixing_residual!(f, rhs, u, scale)
    rhs(f, u, nothing, 0.)
    value = 0.
    @inbounds for i in eachindex(f)
        value = max(value, abs(f[i]) / scale[i])
    end
    value
end

"Physical residual: per-vessel max(|dm_k/dt|/m, |dT/dt|/T) in s^-1 at the current state."
function _mixing_physical_norm(f, rhs)
    residual = 0.
    for (i, node) in enumerate(rhs.network.nodes)
        first = rhs.network.offsets[i]
        first == 0 && continue
        ns = node.initial.gas.n_species
        state = rhs.state_vector[i]
        residual = max(residual, maximum(abs, @view f[first:first+ns-1]) / state.mass,
                       abs(f[first+ns]) / state.temperature)
    end
    residual
end

function _mixing_finish(network, pass, iteration, message, u, residual, physical,
                        history, iterates, save_history)
    states = save_history ? reduce(hcat, iterates) : nothing
    return (; network, state=copy(u), converged=pass, residual, physical_residual=physical,
            iterations=iteration, history=save_history ? history : nothing, states, message)
end

function _solve_mixing_stationary(network::ReactorNetwork; rhs=network_rhs(network),
        residual_tolerance=1e-9, max_iterations=50, max_backtracks=30, save_history=false)
    isfinite(residual_tolerance) && residual_tolerance > 0 ||
        throw(ArgumentError("positive finite tolerance required"))
    max_iterations isa Integer && max_iterations > 0 ||
        throw(ArgumentError("positive iteration count required"))
    max_backtracks isa Integer && max_backtracks >= 0 ||
        throw(ArgumentError("nonnegative backtrack count required"))
    # This example helper is called only with the source mixing_network.
    keys(network.nodes) == (:inlet_air, :inlet_fuel, :mixer, :outlet) ||
        throw(ArgumentError("source mixing network required"))
    network.nodes.mixer.chemistry || throw(ArgumentError("source chemistry must remain enabled"))
    network.nodes.mixer.initial.energy === :adiabatic ||
        throw(ArgumentError("source energy equation required"))
    rhs.network === network || throw(ArgumentError("workspace belongs to another network"))
    u = network_state(network)
    network_isoutofdomain(network, u) && throw(ArgumentError("invalid source state"))
    scale = _mixing_fixed_scales(network, u)
    f, trial, trial_f = similar(u), similar(u), similar(u)
    J = zeros(length(u), length(u))
    system = similar(J)
    scaled_f = similar(u)
    history = MixingIteration[]
    iterates = Vector{Float64}[]
    residual = _mixing_residual!(f, rhs, u, scale)
    physical = _mixing_physical_norm(f, rhs)
    if save_history
        push!(history, MixingIteration(0, 0., 0., residual, physical, 0))
        push!(iterates, copy(u))
    end
    for iteration in 1:max_iterations
        physical <= residual_tolerance && return _mixing_finish(network, true, iteration - 1,
            "converged full network residual", u, residual, physical, history, iterates, save_history)
        network_jacobian!(J, u, rhs, 0.)
        for j in axes(J, 2), i in axes(J, 1)
            J[i, j] *= scale[j] / scale[i]
        end
        @. scaled_f = f / scale
        base_shift = max(residual, 1.)
        accepted = false
        # Shifted solves are continuation steps toward the same f(u)=0;
        # they do not change the final physical residual or initial state.
        for shift in (0., (base_shift * 10.0^k for k in 0:8)...)
            @. system = -J
            for i in axes(system, 1); system[i, i] += shift; end
            direction = try
                lu!(system) \ scaled_f
            catch err
                err isa SingularException || rethrow()
                continue
            end
            all(isfinite, direction) || continue
            model_direction = J * direction
            step = 1.
            for backtrack in 0:max_backtracks
                @. trial = u + step * scale * direction
                if !network_isoutofdomain(network, trial)
                    candidate = _mixing_residual!(trial_f, rhs, trial, scale)
                    predicted_decrease = residual - maximum(abs, scaled_f .+ step .* model_direction)
                    if predicted_decrease > 0 && isfinite(candidate) &&
                       candidate <= residual - 1e-4 * predicted_decrease
                        copyto!(u, trial)
                        copyto!(f, trial_f)
                        residual = candidate
                        physical = _mixing_physical_norm(f, rhs)
                        if save_history
                            push!(iterates, copy(u))
                            push!(history, MixingIteration(iteration, step, shift,
                                                           residual, physical, backtrack))
                        end
                        accepted = true
                        break
                    end
                end
                step /= 2
            end
            accepted && break
        end
        accepted || return _mixing_finish(network, false, iteration - 1,
            "Newton and pseudo-transient line searches exhausted",
            u, residual, physical, history, iterates, save_history)
    end
    return _mixing_finish(network, physical <= residual_tolerance, max_iterations,
        "iteration limit", u, residual, physical, history, iterates, save_history)
end

"""
    solve_mixing_network(gas, air; residual_tolerance=1e-9, max_iterations=50,
                         max_backtracks=30, save_history=false)

Compute the stationary air/methane mixer from Cantera's `reactors/mix1` example.
Provide separate prepared fuel/mixer and air solution models. Each call builds
a fresh reacting, adiabatic network at the original 300 K air state.

The returned tuple contains `network`, `state`, `converged`, `iterations` and
`message`. `physical_residual` is the maximum species mass rate divided by
vessel mass, or temperature rate divided by temperature, in s^-1; convergence
requires it to be at most `residual_tolerance`. `residual` uses fixed initial
scales for the line search. An exhausted iteration or line search returns
`converged=false` and the last accepted state.

Set `save_history=true` for typed iteration records in `history` and a matrix
of corresponding `states`; otherwise both fields are `nothing`.
"""
function solve_mixing_network(gas, air; residual_tolerance=1e-9, max_iterations=50,
                              max_backtracks=30, save_history=false)
    save_history isa Bool || throw(ArgumentError("save_history must be a Bool"))
    gas === air && throw(ArgumentError("distinct gas and air models are required"))
    network = mixing_network(gas, air)
    return _solve_mixing_stationary(network; residual_tolerance, max_iterations,
                                    max_backtracks, save_history)
end
