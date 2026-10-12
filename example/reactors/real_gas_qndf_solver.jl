# Include real_gas_ode_solver.jl (or non_ideal_shock_tube.jl) first.
# Optional dependency of the caller's solver environment; not an Arrhenius dependency.
using OrdinaryDiffEqBDF

struct _ScaledShockTubeRHS{F,V}
    physical::F
    scale::V
    state::V
    derivative::V
end
function (rhs::_ScaledShockTubeRHS)(dz,z,p,t)
    @inbounds for i in eachindex(z)
        rhs.state[i]=z[i]*rhs.scale[i]
    end
    rhs.physical(rhs.derivative,rhs.state,p,t)
    @inbounds for i in eachindex(z)
        dz[i]=rhs.derivative[i]/rhs.scale[i]
    end
    nothing
end

struct _ScaledShockTubeJacobian{F,V,M}
    physical::F
    scale::V
    state::V
    matrix::M
end
function (jac::_ScaledShockTubeJacobian)(J,z,p,t)
    @inbounds for i in eachindex(z)
        jac.state[i]=z[i]*jac.scale[i]
    end
    jac.physical(jac.matrix,jac.state,p,t)
    @inbounds for j in eachindex(z), i in eachindex(z)
        J[i,j]=jac.matrix[i,j]*jac.scale[j]/jac.scale[i]
    end
    nothing
end

struct _ScaledShockTubeDomain{V}
    scale::V
end
function (domain::_ScaledShockTubeDomain)(z,p,t)
    isfinite(z[end]) && z[end]>0 || return true
    @inbounds for i in 1:length(z)-1
        isfinite(z[i]) && z[i]*domain.scale[i]>=-1e-13 || return true
    end
    false
end

"""
    shocktube_qndf_integrator(reactor; trial_policy, end_time=0.005)

Initialize the optional signed/scaled QNDF route and return `(integrator, scale)`.
`trial_policy` must be a `SignedIntegerShockTubeTrials` prepared for this gas.
The transformed state is `[C; T/1000 K]`; returned physical states are `[Y; T]`.
Concentration atol is 1e-19 kmol/m^3, temperature atol is 1e-6 K, rtol is 1e-9,
and the nonlinear convergence coefficient is 0.01. These choices do not change
physical domain gates or project accepted states. No final time stop is imposed.
"""
function shocktube_qndf_integrator(reactor;trial_policy::SignedIntegerShockTubeTrials,end_time=.005)
    isfinite(end_time) && end_time>0 || throw(ArgumentError("positive finite end time required"))
    _validate_trial_policy(trial_policy,reactor)
    u=reactor_state(reactor)
    scale=vcat(reactor.gas.MW./reactor.density,1000.)
    rhs=_ScaledShockTubeRHS(_trial_rhs(trial_policy,reactor),scale,copy(u),zero(u))
    jac=_ScaledShockTubeJacobian(shocktube_ad_jacobian(reactor;trial_policy),scale,copy(u),zeros(length(u),length(u)))
    f=ODEFunction(rhs;jac,tgrad=(out,z,p,t)->fill!(out,zero(eltype(out))))
    ode=ODEProblem(f,u./scale,(0.,Inf),nothing)
    absolute=fill(1e-19,length(u));absolute[end]=1e-6/scale[end]
    integrator=init(ode,QNDF(nlsolve=OrdinaryDiffEqBDF.NLNewton(κ=1//100));reltol=1e-9,abstol=absolute,
        isoutofdomain=_ScaledShockTubeDomain(scale),save_everystep=false,save_start=false,
        save_end=false,dense=false,maxiters=1_000_000)
    return integrator,scale
end

function _shocktube_initialize(::Val{:qndf},reactor,end_time,jacobian,trial_policy)
    jacobian===:ad || throw(ArgumentError("the signed/scaled QNDF adapter requires jacobian=:ad"))
    trial_policy isa SignedIntegerShockTubeTrials ||
        throw(ArgumentError("QNDF requires a mechanism-bound SignedIntegerShockTubeTrials policy"))
    shocktube_qndf_integrator(reactor;trial_policy,end_time)
end
