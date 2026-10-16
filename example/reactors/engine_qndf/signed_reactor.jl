include("trial_rates.jl")

# Signed trial continuation of the source engine. Only chemistry is changed;
# the shared network computes the same signed thermo, flows, piston and ledgers.
struct SignedEngineRHS{F,M,B}
    original::F
    model::M
    network_rhs::B
    volume_rates::Vector{Float64}
    old_omega::Vector{Float64}
    base::Vector{Float64}
    plus::Vector{Float64}
    minus::Vector{Float64}
end
function (rhs::SignedEngineRHS)(du,u,p,t)
    rhs.original(du,u,p,t)
    gas=rhs.model.network.nodes.cylinder.initial.gas
    state=rhs.network_rhs.state_vector[1]
    work=rhs.network_rhs.workspaces[1]
    n=gas.n_species;T=state.temperature;V=state.volume
    copyto!(rhs.old_omega,work.wdot)
    for k in 1:n
        work.C[k]=u[k]/(V*gas.MW[k])
    end
    trial_wdot!(work.wdot,gas.reaction,T,work.C,work.entropy,work.h_mole,work.kinetics;
        rate_multipliers=rhs.model.network.nodes.cylinder.initial.rate_multipliers)
    energy_correction=0.0
    for k in 1:n
        change=V*(work.wdot[k]-rhs.old_omega[k])
        du[k]+=gas.MW[k]*change
        energy_correction+=(work.h_mole[k]-R*T)*change
    end
    du[n+1]-=energy_correction/(state.mass*state.cv)
    nothing
end
function signed_engine_problem(problem)
    original=problem.f
    check_trial_reactions(original.model.network.nodes.cylinder.initial.gas.reaction)
    n=original.model.network.nodes.cylinder.initial.gas.n_species
    rhs=SignedEngineRHS(original,original.model,original.network_rhs,original.volume_rates,
        zeros(n),zero(problem.u0),zero(problem.u0),zero(problem.u0))
    tgrad=(du,u,p,t)->Arrhenius._network_tgrad!(du,u,rhs,t)
    candidate=merge(problem,(;f=rhs,tgrad))
    jac=engine_ad_jacobian(candidate;rate_evaluator=trial_wdot!,clip_trials=false)
    merge(candidate,(;jac))
end
