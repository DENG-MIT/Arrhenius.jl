using Arrhenius, ForwardDiff, LinearAlgebra
include(joinpath(@__DIR__,"real_gas_trial_states.jl"))

"""
    shocktube_ad_jacobian(reactor; chunk=8, trial_policy=ClippedShockTubeTrials())

Prepare a caller-owned `jac!(J,u,p,t)` for a Float64 constant-volume adiabatic
ideal-gas or Redlich-Kwong reactor. The solver environment needs ForwardDiff.

Composition derivatives reuse temperature-only rate factors while differentiating
EOS pressure, fugacity activities, colliders and falloff. The temperature column
uses full automatic differentiation. At zero species the clipped policy uses
the limit from positive composition; the signed policy retains exact zeros.
Negative trial species follow the RHS's clipped
concentration branch by default. An explicit `SignedIntegerShockTubeTrials` policy
uses its source-matched signed continuation. Neither policy modifies the ODE state.
At a thermo-polynomial split, the derivative follows the selected coefficient
region. Blowers-Masel rates are outside this helper's validated scope.

The returned callback owns mutable work arrays and must not be shared by
simultaneously running solves. This helper uses internal Arrhenius workspaces;
include it from the same checkout as the Arrhenius package being used.
"""
function shocktube_ad_jacobian(reactor;chunk=8,trial_policy::ShockTubeTrialPolicy=ClippedShockTubeTrials())
    chunk isa Integer && chunk>0 || throw(ArgumentError("positive integer chunk size required"))
    _validate_trial_policy(trial_policy,reactor)
    return _shocktube_ad_jacobian(reactor,Val(min(chunk,reactor.gas.n_species)),trial_policy)
end
function _shocktube_ad_jacobian(reactor,::Val{N},trial_policy::P) where {N,P<:ShockTubeTrialPolicy}
    reactor.energy === :adiabatic || throw(ArgumentError("adiabatic reactor required"))
    reactor isa IdealGasReactor && reactor.constraint !== :constant_volume && throw(ArgumentError("constant-volume reactor required"))
    gas=reactor.gas
    isempty(gas.reaction.blowers_masel.reaction_indices) || throw(ArgumentError("this Jacobian does not support Blowers-Masel rates"))
    n=gas.n_species
    state=reactor_state(reactor)
    eltype(state)===Float64 || throw(ArgumentError("Float64 reactor state required"))
    x=state[1:n]
    out=zero(state)
    cfg=ForwardDiff.JacobianConfig(nothing,out,x,ForwardDiff.Chunk{N}())
    D=eltype(typeof(cfg))
    float_rhs=_trial_rhs(trial_policy,reactor)
    floats=float_rhs.workspace
    cache=Arrhenius._KineticsTemperatureCache(gas.reaction)
    work=reactor isa RedlichKwongReactor ? RealGasKineticsWorkspace(gas,reactor.model,D) : Arrhenius.ReactorWorkspace(gas,D)
    current_temperature=Ref(state[end])
    if reactor isa RedlichKwongReactor
        composition_rhs = function (du,Y)
            T=current_temperature[]
            total=zero(eltype(Y))
            @inbounds for k in 1:n
                work.X[k]=_trial_concentration(trial_policy,Y[k])/gas.MW[k]
                total+=work.X[k]
            end
            work.X ./= total
            properties=_trial_caloric!(trial_policy,work.thermo,reactor.model,T,reactor.density,work.X)
            @inbounds for k in 1:n
                work.C[k]=work.X[k]/properties.v
                work.activity[k]=work.X[k]*properties.P/(R*T)*exp(properties.lnphi[k])
            end
            _trial_wdot!(trial_policy,work.wdot,gas.reaction,T,work.C,floats.s0,floats.h0,work.kinetics;
                temperature_cache=cache,pressure=properties.P,activity_concentrations=work.activity,
                rate_multipliers=reactor.rate_multipliers)
            source=zero(eltype(Y))
            @inbounds for k in 1:n
                du[k]=work.wdot[k]*gas.MW[k]/reactor.density
                source+=work.wdot[k]*properties.u_TV[k]
            end
            du[end]=-source/(reactor.density*properties.cv_mass)
            nothing
        end
    else
        composition_rhs = function (du,Y)
            T=current_temperature[]
            @inbounds for k in 1:n
                work.C[k]=_trial_concentration(trial_policy,Y[k])*reactor.density/gas.MW[k]
            end
            _trial_wdot!(trial_policy,work.wdot,gas.reaction,T,work.C,floats.entropy,floats.h_mole,work.kinetics;
                temperature_cache=cache,rate_multipliers=reactor.rate_multipliers)
            capacity,source=zero(eltype(Y)),zero(eltype(Y))
            @inbounds for k in 1:n
                du[k]=work.wdot[k]*gas.MW[k]/reactor.density
                capacity+=Y[k]*(floats.cp_R[k]-1)/gas.MW[k]
                source+=work.wdot[k]*(floats.h_mole[k]-R*T)
            end
            du[end]=-source/(reactor.density*R*capacity)
            nothing
        end
    end
    temp=[state[end]]
    tcfg=ForwardDiff.JacobianConfig(nothing,out,temp,ForwardDiff.Chunk{1}())
    DT=eltype(typeof(tcfg))
    dual_state=DT.(state)
    temp_work=reactor isa RedlichKwongReactor ? RealGasKineticsWorkspace(gas,reactor.model,DT) : Arrhenius.ReactorWorkspace(gas,DT)
    temp_rhs=_trial_rhs(trial_policy,reactor,temp_work,dual_state)
    temperature_rhs = function (du,temperature)
        copyto!(dual_state,state)
        dual_state[end]=temperature[1]
        temp_rhs(du,dual_state,nothing,0.)
    end
    # Freeze the branch-selected callable in a fresh binding so it is stored
    # with its concrete closure type rather than a captured Core.Box.
    return let composition_rhs=composition_rhs
      function (J,u,p,t)
        length(u)==n+1 && size(J)==(n+1,n+1) || throw(DimensionMismatch("Jacobian and state dimensions must match the reactor"))
        copyto!(state,u)
        @inbounds for k in 1:n
            # Only the clipped branch needs the positive-side convention.
            # Signed integer-order trials retain exact absent-element zeros.
            trial_policy isa ClippedShockTubeTrials && state[k]==0 && (state[k]=1e-100)
        end
        copyto!(x,view(state,1:n))
        T=state[end]
        current_temperature[]=T
        temp[1]=T
        float_rhs(out,state,nothing,0.)
        h=reactor isa RedlichKwongReactor ? floats.h0 : floats.h_mole
        s=reactor isa RedlichKwongReactor ? floats.s0 : floats.entropy
        # Only T-dependent Arrhenius and reference-Kc factors enter this cache.
        # Composition, EOS pressure, colliders, falloff and mass action are
        # reevaluated using dual values in composition_rhs.
        wdot!(floats.wdot,gas.reaction,T,floats.C,s,h,floats.kinetics;
            temperature_cache=cache,rate_multipliers=reactor.rate_multipliers)
        ForwardDiff.jacobian!(view(J,:,1:n),composition_rhs,out,x,cfg)
        ForwardDiff.jacobian!(view(J,:,n+1:n+1),temperature_rhs,out,temp,tcfg)
        nothing
      end
    end
end
