using SciMLBase, OrdinaryDiffEqBDF, LinearAlgebra, SparseArrays

if !isdefined(@__MODULE__, :GuardedStructuredJacobian)
    include(joinpath(@__DIR__, "guarded_jacobian.jl"))
end
if !isdefined(@__MODULE__, :KLUCoreStructuredSolver)
    include(joinpath(@__DIR__, "klu_structured_solver.jl"))
end
using .KLUCoreStructuredSolver: new_structured_solver, factor_shift!, solve_factored!, solver_report

# LinearSolveFunction defaults to matrix-free. This owned callable implements
# a direct solve and requires QNDF to construct both J and concrete W.
struct ConcreteStructuredSolveFunction{F}
    f::F
end
(call::ConcreteStructuredSolveFunction)(args...; kwargs...) = call.f(args...; kwargs...)
LS.needs_concrete_A(::LS.LinearSolveFunction{<:ConcreteStructuredSolveFunction}) = true

mutable struct QNDFStructuredAdapter{G,S}
    guard::G
    solver::S
    latest_J::Matrix{Float64}
    jacobian_generations::Int
    refreshed_generations::Int
    refreshed_generation::Int
    current_mode::Symbol
    current_generation::Int
    factored_generation::Int
    factored_sigma::Float64
    factored_mode::Symbol
    pending_newW::Bool
    pending_W::Union{Nothing,Matrix{Float64}}
    pending_sigma::Float64
    pending_generation::Int
    pending_mode::Symbol
    init_hook_calls::Int
    real_hook_calls::Int
    new_w_calls::Int
    factorizations::Int
    solves::Int
    rhs_reuses::Int
    fd_factorizations::Int
    fd_solves::Int
    relation_checks::Int
    new_w_sigmas::Vector{Float64}
    new_w_generations::Vector{Int}
    new_w_modes::Vector{Symbol}
end

function _initial_structured_data(guard)
    jac = guard.analytic
    ns = jac.reactor.gas.n_species
    Dict{String,Any}(
        "J" => zeros(ns+1,ns+1),
        "core_colptr" => copy(jac.mass_core.colptr),
        "core_rowval" => copy(jac.mass_core.rowval),
        "core_nzval" => copy(jac.mass_core.nzval),
        "q" => copy(jac.q),
        "rank_vector" => copy(jac.rank_vector),
        "temperature_column" => copy(jac.temperature_column),
        "energy_row" => copy(jac.energy_row),
    )
end

function qndf_structured_adapter(guard::GuardedStructuredJacobian)
    n = guard.analytic.reactor.gas.n_species + 1
    solver = new_structured_solver(_initial_structured_data(guard))
    QNDFStructuredAdapter(guard,solver,zeros(n,n),0,0,0,:none,0,0,NaN,:none,
        false,nothing,NaN,0,:none,0,0,0,0,0,0,0,0,0,
        Float64[],Int[],Symbol[])
end

function _refresh_analytic!(adapter::QNDFStructuredAdapter, J)
    jac, solver = adapter.guard.analytic, adapter.solver
    jac.mass_core.colptr == solver.colptr0 || error("analytic sparse colptr changed")
    jac.mass_core.rowval == solver.rowval0 || error("analytic sparse rowval changed")
    solver.shifted.colptr == solver.colptr0 || error("solver sparse colptr changed")
    solver.shifted.rowval == solver.rowval0 || error("solver sparse rowval changed")
    copyto!(solver.base_nzval,jac.mass_core.nzval)
    copyto!(solver.nzval0,jac.mass_core.nzval)
    copyto!(solver.q,jac.q); copyto!(solver.q0,jac.q)
    copyto!(solver.u,jac.rank_vector); copyto!(solver.u0,jac.rank_vector)
    copyto!(solver.t,@view(jac.temperature_column[1:solver.ns]))
    copyto!(solver.t0,@view(jac.temperature_column[1:solver.ns]))
    copyto!(solver.e,jac.energy_row); copyto!(solver.e0,jac.energy_row)
    solver.a = jac.temperature_column[end]
    copyto!(solver.dense.original,J)
    adapter.refreshed_generations += 1
    adapter.refreshed_generation = adapter.current_generation
    return nothing
end

function (adapter::QNDFStructuredAdapter)(J,state,p,t)
    before = adapter.guard.generation
    adapter.guard(J,state,p,t)
    adapter.guard.generation == before + 1 || error("guard generation did not advance once")
    copyto!(adapter.latest_J,J)
    adapter.jacobian_generations += 1
    adapter.current_generation = adapter.guard.generation
    adapter.current_mode = adapter.guard.last_used_fallback ? :fd : :structured
    copyto!(adapter.solver.dense.original,J)
    adapter.current_mode === :structured && _refresh_analytic!(adapter,J)
    return nothing
end

function _exact_W(J,W,sigma)
    size(J) == size(W) || return false
    @inbounds for column in axes(J,2), row in axes(J,1)
        expected = row == column ? muladd(-1.0,sigma,J[row,column]) : J[row,column]
        isequal(W[row,column],expected) || return false
    end
    return true
end

function adapter_precs(adapter::QNDFStructuredAdapter,integrator_ref)
    function precs(W,du,u,p,t,newW,Plprev,Prprev,solverdata)
        if integrator_ref[] === nothing
            newW === nothing || error("unexpected QNDF initialization marker")
            adapter.init_hook_calls += 1
            return nothing,nothing
        end
        newW isa Bool || error("real QNDF hook lacks Boolean newW marker")
        adapter.real_hook_calls += 1
        nl = integrator_ref[].cache.nlsolver.cache
        nl.linsolve.A === W || error("QNDF linear matrix alias changed")
        sigma = inv(nl.W_γdt)
        isfinite(sigma) && sigma > 0 || error("QNDF shift must be finite and positive")
        adapter.current_generation > 0 || error("W prepared before guarded Jacobian")
        if newW
            W isa Matrix{Float64} || error("adapter requires a concrete Float64 W")
            _exact_W(adapter.latest_J,W,sigma) || error("QNDF W differs from J-sigma*I")
            adapter.relation_checks += 1
            adapter.new_w_calls += 1
            adapter.pending_newW = true
            adapter.pending_W = W
            adapter.pending_sigma = sigma
            adapter.pending_generation = adapter.current_generation
            adapter.pending_mode = adapter.current_mode
            push!(adapter.new_w_sigmas,sigma)
            push!(adapter.new_w_generations,adapter.current_generation)
            push!(adapter.new_w_modes,adapter.current_mode)
        else
            adapter.pending_newW = false
            adapter.pending_W === W || error("RHS reuse changed the QNDF W object")
            adapter.factored_generation > 0 || error("RHS reuse preceded factorization")
            adapter.factored_generation == adapter.current_generation || error("RHS reuse changed Jacobian generation")
            isequal(sigma,adapter.factored_sigma) || error("RHS reuse changed the actual QNDF shift")
            adapter.factored_mode === adapter.current_mode || error("RHS reuse changed guarded Jacobian mode")
        end
        return nothing,nothing
    end
    return precs
end

function adapter_linsolve(adapter::QNDFStructuredAdapter,integrator_ref)
    function solve_linear(A,b,x,p,isfresh,Pl,Pr,cacheval;kwargs...)
        integrator_ref[] === nothing && error("linear solve before integrator reference wiring")
        nl = integrator_ref[].cache.nlsolver.cache
        nl.linsolve.A === A || error("QNDF solve matrix alias changed")
        if isfresh
            adapter.pending_newW || error("fresh linear cache without newW hook")
            adapter.pending_W === A || error("fresh solve differs from hooked W")
            adapter.pending_generation == adapter.current_generation ||
                error("stale Jacobian generation at factorization")
            if adapter.pending_mode === :structured
                adapter.refreshed_generation == adapter.pending_generation ||
                    error("structured generation was not refreshed")
                factor_shift!(adapter.solver,adapter.pending_sigma)
            elseif adapter.pending_mode === :fd
                factor_shift!(adapter.solver.dense,adapter.pending_sigma)
                adapter.fd_factorizations += 1
            else
                error("unknown guarded Jacobian mode")
            end
            adapter.factorizations += 1
            adapter.factored_generation = adapter.pending_generation
            adapter.factored_sigma = adapter.pending_sigma
            adapter.factored_mode = adapter.pending_mode
        else
            !adapter.pending_newW || error("newW hook was not followed by a fresh solve")
            adapter.factored_generation == adapter.current_generation ||
                error("RHS solve would reuse a stale Jacobian generation")
            adapter.rhs_reuses += 1
        end
        if adapter.factored_mode === :structured
            solve_factored!(x,adapter.solver,b)
        elseif adapter.factored_mode === :fd
            solve_factored!(x,adapter.solver.dense,b)
            adapter.fd_solves += 1
        else
            error("linear solve has no factored mode")
        end
        adapter.solves += 1
        nl.linsolve.isfresh = false
        adapter.pending_newW = false
        return x
    end
    return LS.LinearSolveFunction(ConcreteStructuredSolveFunction(solve_linear))
end

function _lifecycle_valid(adapter::QNDFStructuredAdapter)
    adapter.real_hook_calls == adapter.solves &&
        adapter.new_w_calls == adapter.factorizations &&
        adapter.solves == adapter.factorizations + adapter.rhs_reuses &&
        adapter.relation_checks == adapter.new_w_calls &&
        adapter.jacobian_generations == adapter.guard.generation &&
        adapter.fd_factorizations <= adapter.factorizations &&
        adapter.fd_solves <= adapter.solves && !adapter.pending_newW
end

function adapter_report(adapter::QNDFStructuredAdapter)
    lifecycle = _lifecycle_valid(adapter)
    Dict{String,Any}(
        "init_hook_calls"=>adapter.init_hook_calls,
        "real_hook_calls"=>adapter.real_hook_calls,
        "new_w_calls"=>adapter.new_w_calls,
        "factorizations"=>adapter.factorizations,
        "solves"=>adapter.solves,
        "rhs_reuses"=>adapter.rhs_reuses,
        "jacobian_generations"=>adapter.jacobian_generations,
        "structured_refreshes"=>adapter.refreshed_generations,
        "last_structured_refresh_generation"=>adapter.refreshed_generation,
        "fd_factorizations"=>adapter.fd_factorizations,
        "fd_solves"=>adapter.fd_solves,
        "relation_checks"=>adapter.relation_checks,
        "new_w_sigmas"=>copy(adapter.new_w_sigmas),
        "new_w_generations"=>copy(adapter.new_w_generations),
        "new_w_modes"=>string.(adapter.new_w_modes),
        "last_factored_generation"=>adapter.factored_generation,
        "last_factored_sigma"=>adapter.factored_sigma,
        "last_factored_mode"=>string(adapter.factored_mode),
        "guard_analytic_calls"=>adapter.guard.analytic_calls,
        "guard_fallback_calls"=>adapter.guard.fallback_calls,
        "structured_solver"=>solver_report(adapter.solver),
        "lifecycle_passed"=>lifecycle,
    )
end
