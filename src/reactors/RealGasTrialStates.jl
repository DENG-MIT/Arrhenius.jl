# Redlich–Kwong expressions adapted from Cantera RedlichKwongMFTP.cpp,
# revision 726522be4e2a13454d8415b7ef799d621f665cf3.
# https://github.com/Cantera/cantera/blob/726522be4e2a13454d8415b7ef799d621f665cf3/src/thermo/RedlichKwongMFTP.cpp
#
# Copyright (c) 2001-2009, California Institute of Technology
# All rights reserved.
# Copyright (c) 2009 Sandia Corporation. Under the terms of Contract
# AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain
# rights in this software.
# Copyright (c) 2011-2026, Cantera Developers. All rights reserved.
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# - Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# - Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# - Neither the name of the California Institute of Technology, Sandia
#   Corporation nor the names of other contributors may be used to endorse or
#   promote products derived from this software without specific prior written
#   permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
# Signed C1/C2/C3 products follow Cantera StoichManager.h at the same revision.
# https://github.com/Cantera/cantera/blob/726522be4e2a13454d8415b7ef799d621f665cf3/include/cantera/kinetics/StoichManager.h
# Local numerical continuation for the caller-owned shock-tube solver.
# Signed values are solver trials; physical state constructors remain unchanged.
using Arrhenius, LinearAlgebra, SparseArrays

abstract type ShockTubeTrialPolicy end

"Use the package reactor RHS, which clips negative trial concentrations."
struct ClippedShockTubeTrials <: ShockTubeTrialPolicy end

# Preparation compares the numerical mechanism, including nested rate and thermo
# data. Species counts or names alone cannot establish the mass-action policy.
function _trial_data_equal(a,b)
    typeof(a) === typeof(b) || return false
    if a isa AbstractArray || a isa Number || a isa AbstractString || a isa Symbol || isnothing(a)
        return isequal(a,b)
    end
    return all(name -> _trial_data_equal(getfield(a,name),getfield(b,name)),fieldnames(typeof(a)))
end

"""
    SignedIntegerShockTubeTrials(gas, mechanism_path)

Prepare the bounded signed trial policy for elementary, third-body, Lindemann,
Troe and PLOG reactions with integer stoichiometric orders of at most three.
Explicit custom orders are rejected, even when equal to the default orders:
Cantera selects a different negative-trial branch from that YAML metadata.

The supplied YAML and hash-verified native sidecar must describe `gas` exactly.
The policy is bound to that gas and checks its unchanged numerical data when
callbacks are prepared. Do not mutate a mechanism during an active solve.
This is a local solver continuation, not a thermodynamic signed-state API.
"""
struct SignedIntegerShockTubeTrials{G} <: ShockTubeTrialPolicy
    gas::G
    verified::G
    mechanism_sha256::String
    function SignedIntegerShockTubeTrials(gas,mechanism_path::AbstractString)
        input=Arrhenius.YAML.load_file(mechanism_path)
        reactions=get(input,"reactions",nothing)
        reactions isa AbstractVector || throw(ArgumentError("an inline reaction list is required"))
        for (i,reaction) in enumerate(reactions)
            reaction isa AbstractDict || throw(ArgumentError("inline reaction dictionaries required"))
            haskey(reaction,"orders") && throw(ArgumentError("reaction $i has explicit custom orders"))
            kind=get(reaction,"type","elementary")
            kind in ("elementary","three-body","falloff","pressure-dependent-Arrhenius") ||
                throw(ArgumentError("unsupported signed-trial reaction type: $kind"))
            (haskey(reaction,"SRI") || haskey(reaction,"Tsang")) &&
                throw(ArgumentError("only Lindemann and Troe falloff trial rates are supported"))
        end
        sidecar=Arrhenius.npzread(mechanism_path*".npz")
        haskey(sidecar,"source_sha256_utf8") ||
            throw(ArgumentError("signed trials require a YAML-hash-bound native sidecar"))
        verified=CreateSolution(mechanism_path)
        _trial_data_equal(gas,verified) ||
            throw(ArgumentError("supplied mechanism does not match the reactor gas numerical data"))
        _check_integer_trial_reactions(gas.reaction)
        new{typeof(gas)}(gas,verified,bytes2hex(Arrhenius.SHA.sha256(read(mechanism_path))))
    end
end

_validate_trial_policy(::ClippedShockTubeTrials,reactor)=nothing
function _validate_trial_policy(policy::SignedIntegerShockTubeTrials,reactor)
    policy.gas === reactor.gas || throw(ArgumentError("trial policy belongs to a different gas"))
    _trial_data_equal(policy.gas,policy.verified) ||
        throw(ArgumentError("gas numerical data changed after trial-policy preparation"))
    if reactor isa RedlichKwongReactor
        _trial_data_equal(reactor.model.reference,reactor.gas.thermo) ||
            throw(ArgumentError("RK and kinetic reference thermodynamics must match"))
    end
    nothing
end

_trial_concentration(::ClippedShockTubeTrials,c)=max(c,zero(c))
_trial_concentration(::SignedIntegerShockTubeTrials,c)=c
_trial_caloric!(::ClippedShockTubeTrials,args...)=redlich_kwong_properties!(args...)
_trial_caloric!(::SignedIntegerShockTubeTrials,args...)=_signed_trial_caloric!(args...)
_trial_wdot!(::ClippedShockTubeTrials,args...;kwargs...)=wdot!(args...;kwargs...)
_trial_wdot!(::SignedIntegerShockTubeTrials,args...;kwargs...)=_integer_trial_wdot!(args...;kwargs...)
_trial_rhs(::ClippedShockTubeTrials,reactor)=reactor_rhs(reactor)
_trial_rhs(::SignedIntegerShockTubeTrials,reactor)=_signed_reactor_rhs(reactor)
_trial_rhs(::SignedIntegerShockTubeTrials,reactor,work,state)=_SignedTrialRHS(reactor,work)
_trial_rhs(::ClippedShockTubeTrials,reactor,work,state)=
    Arrhenius.ReactorRHS(reactor,work,copy(state),zero(state),zero(state),zero(state))

function _check_integer_trial_reactions(reaction)
    isempty(reaction.blowers_masel.reaction_indices) || throw(ArgumentError("BM trial rates are unsupported"))
    reaction.reactant_orders==reaction.reactant_stoich_coeffs || throw(ArgumentError("custom orders unsupported"))
    all(isinteger,nonzeros(reaction.reactant_orders)) || throw(ArgumentError("integer reactant orders required"))
    all(isinteger,nonzeros(reaction.product_stoich_coeffs)) || throw(ArgumentError("integer product orders required"))
    all(vec(sum(reaction.reactant_orders;dims=1)).<=3) || throw(ArgumentError("reactant orders above three unsupported"))
    all(vec(sum(reaction.product_stoich_coeffs;dims=1))[reaction.is_reversible].<=3) ||
        throw(ArgumentError("active reverse orders above three unsupported"))
    nothing
end

# Read live CSC entries on every call; retain the generic ordered loop below.
# Only the common (1), (2), and (1,1) columns receive explicit scalar paths.
@inline function _integer_trial_product(factor,activity,rows,values,indices)
    first_index,last_index=first(indices),last(indices)
    @inbounds if first_index==last_index
        order=values[first_index]
        if order==one(order)
            concentration=activity[rows[first_index]]
            value=convert(promote_type(typeof(concentration),typeof(order)),concentration)
            return factor*value
        elseif order==oftype(order,2)
            concentration=activity[rows[first_index]]
            value=convert(promote_type(typeof(concentration),typeof(order)),concentration)
            result=factor*(value*value)
            return concentration<0 ? zero(result) : result
        end
    elseif last_index==first_index+1
        order1,order2=values[first_index],values[last_index]
        if order1==one(order1) && order2==one(order2)
            c1,c2=activity[rows[first_index]],activity[rows[last_index]]
            v1=convert(promote_type(typeof(c1),typeof(order1)),c1)
            v2=convert(promote_type(typeof(c2),typeof(order2)),c2)
            result=(factor*v1)*v2
            return c1<0 && c2<0 ? zero(result) : result
        end
    end
    negative_factors=0.
    @inbounds for j in indices
        concentration=activity[rows[j]]
        factor*=Arrhenius._concentration_power(concentration,values[j])
        negative_factors+=concentration<0 ? values[j] : 0.
    end
    negative_factors>=2 && (factor=zero(factor))
    factor
end

function _integer_trial_wdot!(output,reaction,T,C,S0,h,work;activity_concentrations=nothing,reverse_plan=nothing,kwargs...)
    # Reuse native Arrhenius / equilibrium / collider / falloff / PLOG factors.
    if isnothing(reverse_plan)
        wdot!(output,reaction,T,C,S0,h,work;activity_concentrations,get_rate_constants=true,kwargs...)
    else
        Arrhenius._rate_factors!(reaction,T,C,S0,h,work,reverse_plan;kwargs...)
    end
    activity=isnothing(activity_concentrations) ? C : activity_concentrations
    reactants,products=reaction.reactant_orders,reaction.product_stoich_coeffs
    ri,rv=rowvals(reactants),nonzeros(reactants)
    pi,pv=rowvals(products),nonzeros(products)
    @inbounds for i in 1:reaction.n_reactions
        forward,reverse=work.kf[i],work.kr[i]
        forward=_integer_trial_product(forward,activity,ri,rv,nzrange(reactants,i))
        if reaction.is_reversible[i]
            reverse=_integer_trial_product(reverse,activity,pi,pv,nzrange(products,i))
        end
        work.kf[i],work.kr[i]=forward,reverse
        work.rates_of_progress[i]=forward-reverse
    end
    mul!(output,reaction.vk,work.rates_of_progress)
    output
end

# Caloric and EOS expressions have a smooth local continuation around X=0.
# Mixing entropy has no real-valued signed extension and is not evaluated here.
function _signed_trial_caloric!(w,m,T,rho,X)
    n=length(m.MW)
    length(X)==n && length(w.Ak)==n || throw(DimensionMismatch("trial RK dimensions"))
    isfinite(T) && T>0 && isfinite(rho) && rho>0 || throw(DomainError((T,rho)))
    all(isfinite,X) && abs(sum(X)-1)<=1e-10 || throw(DomainError(X))
    a,aT,b=Arrhenius._rk_mixing!(w,m,T,X)
    MW=dot(X,m.MW);v=MW/rho
    v>b || throw(DomainError(v))
    P=Arrhenius._rk_pressure(T,v,a,b)
    dpdv=Arrhenius._rk_dpdv(T,v,a,b)
    P>0 && dpdv<0 || throw(DomainError((T,rho)))
    sqt,RT=sqrt(T),R*T
    @inbounds for i in 1:n
        if isnothing(m.reference.extra)
            coeff=Arrhenius._nasa7_coefficients(m.reference,i,T)
            w.cp0[i]=coeff[i,1]+T*(coeff[i,2]+T*(coeff[i,3]+T*(coeff[i,4]+T*coeff[i,5])))
            w.h0[i]=coeff[i,1]+T*(coeff[i,2]/2+T*(coeff[i,3]/3+T*(coeff[i,4]/4+T*coeff[i,5]/5)))+coeff[i,6]/T
            w.s0[i]=coeff[i,1]*log(T)+T*(coeff[i,2]+T*(coeff[i,3]/2+T*(coeff[i,4]/3+T*coeff[i,5]/4)))+coeff[i,7]
        else
            w.cp0[i],w.h0[i],w.s0[i]=Arrhenius._extended_thermo(m.reference,i,T)
        end
    end
    L=log1p(b/v);Lb=b==0 ? inv(v) : L/b
    F=T*aT-1.5a
    cv=R*(dot(X,w.cp0)-1)+Lb/sqt*(0.75a/T-aT)
    cv>0 || throw(DomainError(cv))
    @inbounds for i in 1:n
        bi=m.b[i]
        if b==0
            w.u_TV[i]=RT*(w.h0[i]-1);w.lnphi[i]=zero(T)
            continue
        end
        Sk=2T*w.dAk[i]-3w.Ak[i]
        C=-Lb+1/(v+b)
        w.u_TV[i]=RT*(w.h0[i]-1)+(L*Sk+bi*F*C)/(b*sqt)
        w.lnphi[i]=log(RT/((v-b)*P))+bi/(v-b)-2w.Ak[i]*Lb/(sqt*RT)+
            a*bi*Lb/(b*sqt*RT)-a*bi/(b*sqt*(v+b)*RT)
    end
    (;T,P,rho,X,MW,v,cv_mass=cv/MW,u_TV=w.u_TV,lnphi=w.lnphi)
end

struct _SignedTrialRHS{R,W,P}
    reactor::R
    workspace::W
    reverse_plan::P
end
_SignedTrialRHS(reactor,work)=_SignedTrialRHS(reactor,work,Arrhenius._ReversibleRatePlan(reactor.gas.reaction))
function _signed_reactor_rhs(reactor,::Type{T}=Float64) where T
    _check_integer_trial_reactions(reactor.gas.reaction)
    reactor.energy===:adiabatic || error("adiabatic trial reactor required")
    if reactor isa IdealGasReactor
        reactor.constraint===:constant_volume || error("constant-volume trial reactor required")
        work=Arrhenius.ReactorWorkspace(reactor.gas,T)
    else
        work=RealGasKineticsWorkspace(reactor.gas,reactor.model,T)
    end
    _SignedTrialRHS(reactor,work)
end

function (rhs::_SignedTrialRHS{<:IdealGasReactor})(du,u,p,t)
    reactor,work=rhs.reactor,rhs.workspace
    gas=reactor.gas;T=u[end];density=reactor.density
    inverse_mw=zero(T)
    @inbounds for k in 1:gas.n_species
        inverse_mw+=u[k]/gas.MW[k]
        work.C[k]=u[k]*density/gas.MW[k]
    end
    @inbounds for k in 1:gas.n_species
        work.X[k]=u[k]/(gas.MW[k]*inverse_mw)
    end
    P=density*R*T*inverse_mw
    cal_cp_R!(work.cp_R,gas,T,P,work.X)
    cal_h_RT!(work.h_mole,gas,T,P,work.X)
    cal_s0_R!(work.entropy,gas,T,P,work.X)
    @inbounds for k in 1:gas.n_species
        work.h_mole[k]*=R*T;work.entropy[k]*=R
    end
    _integer_trial_wdot!(work.wdot,gas.reaction,T,work.C,work.entropy,work.h_mole,work.kinetics;
        rate_multipliers=reactor.rate_multipliers,reverse_plan=rhs.reverse_plan)
    capacity,source=zero(T),zero(T)
    @inbounds for k in 1:gas.n_species
        du[k]=work.wdot[k]*gas.MW[k]/density
        capacity+=u[k]*(work.cp_R[k]-1)/gas.MW[k]
        source+=work.wdot[k]*(work.h_mole[k]-R*T)
    end
    du[end]=-source/(density*R*capacity)
    nothing
end

function (rhs::_SignedTrialRHS{<:RedlichKwongReactor})(du,u,p,t)
    reactor,work=rhs.reactor,rhs.workspace
    gas=reactor.gas;T=u[end];n=gas.n_species
    total=zero(T)
    @inbounds for k in 1:n
        work.X[k]=u[k]/gas.MW[k]
        total+=work.X[k]
    end
    work.X./=total
    properties=_signed_trial_caloric!(work.thermo,reactor.model,T,reactor.density,work.X)
    @inbounds for k in 1:n
        work.C[k]=work.X[k]/properties.v
        work.activity[k]=work.X[k]*properties.P/(R*T)*exp(properties.lnphi[k])
        work.h0[k]=R*T*work.thermo.h0[k];work.s0[k]=R*work.thermo.s0[k]
    end
    _integer_trial_wdot!(work.wdot,gas.reaction,T,work.C,work.s0,work.h0,work.kinetics;
        pressure=properties.P,activity_concentrations=work.activity,rate_multipliers=reactor.rate_multipliers,reverse_plan=rhs.reverse_plan)
    source=zero(T)
    @inbounds for k in 1:n
        du[k]=work.wdot[k]*gas.MW[k]/reactor.density
        source+=work.wdot[k]*properties.u_TV[k]
    end
    du[end]=-source/(reactor.density*properties.cv_mass)
    nothing
end
