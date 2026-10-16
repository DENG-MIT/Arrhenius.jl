# Exact invariant-coordinate restriction for the source engine topology.
mutable struct EngineInvariantCoordinates
    indices::Vector{Int}
    removed::Vector{Int}
    absent_elements::Vector{String}
    full::Vector{Float64}
    derivative::Vector{Float64}
    jacobian::Matrix{Float64}
end

function engine_invariant_coordinates(problem)
    network=problem.f.model.network;gas=network.nodes.cylinder.initial.gas;ns=gas.n_species
    keys(network.nodes)==(:cylinder,:inlet,:injector,:outlet,:ambient) || error("source engine topology required")
    network.offsets==[1,0,0,0,0] || error("one source cylinder required")
    length(problem.u0)==ns+5 || error("source species/T/volume/three-ledger state required")
    inlet,fuel,outlet=network.flows
    inlet isa Valve && fuel isa MassFlowController && outlet isa Valve || error("source flow devices required")
    ((inlet.upstream,inlet.downstream),(fuel.upstream,fuel.downstream),(outlet.upstream,outlet.downstream))==
        ((:inlet,:cylinder),(:injector,:cylinder),(:cylinder,:outlet)) || error("source one-way flow connections required")
    inlet.time_function(0.)==inlet.time_function(last(problem.tspan)) || error("inlet gate is not frozen")
    fuel.mdot isa Real || error("fuel gate must be frozen within an interval")
    support=problem.u0[1:ns].!=0
    feeds=Vector{Float64}[]
    inlet.K*inlet.time_function(0.)>0 && push!(feeds,network.nodes.inlet.initial.mass_fractions)
    fuel.mdot>0 && push!(feeds,network.nodes.injector.initial.mass_fractions)
    for feed in feeds;support .|= feed.!=0;end
    absent=findall(e->!any(support .& (gas.ele_matrix[e,:].!=0)),axes(gas.ele_matrix,1))
    removed=findall(k->any(gas.ele_matrix[absent,k].!=0),1:ns)
    # Element cancellation or a small nonzero species never authorizes removal.
    all(iszero,problem.u0[removed]) || error("cannot remove a nonzero initial species")
    all(feed->all(iszero,feed[removed]),feeds) || error("removed species occurs in an active feed")
    indices=setdiff(collect(eachindex(problem.u0)),removed)
    n=length(problem.u0)
    EngineInvariantCoordinates(indices,removed,string.(gas.elements[absent]),zeros(n),zeros(n),zeros(n,n))
end

function engine_expand!(coordinates::EngineInvariantCoordinates,u)
    fill!(coordinates.full,0.)
    @inbounds for (j,k) in enumerate(coordinates.indices);coordinates.full[k]=u[j];end
    coordinates.full
end
engine_expand!(::Nothing,u)=u

struct EngineRestrictedRHS{F}
    f::F
    coordinates::EngineInvariantCoordinates
end
function (rhs::EngineRestrictedRHS)(du,u,p,t)
    c=rhs.coordinates
    rhs.f(c.derivative,engine_expand!(c,u),p,t)
    all(iszero,view(c.derivative,c.removed)) || error("the full RHS leaves the certified invariant subspace")
    @inbounds for (j,k) in enumerate(c.indices);du[j]=c.derivative[k];end
    nothing
end
function engine_restricted_jacobian(c,full_jacobian)
    function restricted_jacobian!(J,u,p,t)
        full_jacobian(c.jacobian,engine_expand!(c,u),p,t)
        all(iszero,view(c.jacobian,c.removed,c.indices)) ||
            error("the full Jacobian leaves the certified invariant subspace")
        @inbounds for j in eachindex(c.indices),i in eachindex(c.indices)
            J[i,j]=c.jacobian[c.indices[i],c.indices[j]]
        end
        nothing
    end
end
function engine_restricted_tgrad(c,full_tgrad)
    function restricted_tgrad!(du,u,p,t)
        full_tgrad(c.derivative,engine_expand!(c,u),p,t)
        all(iszero,view(c.derivative,c.removed)) || error("the full time derivative leaves the certified invariant subspace")
        @inbounds for (j,k) in enumerate(c.indices);du[j]=c.derivative[k];end
        nothing
    end
end

# The original error norm includes all full-state coordinates. Omitted
# coordinates have exactly zero error, so retain its original RMS denominator.
struct FullEngineRMS
    length::Int
end
# DiffEqBase's vector-atol residual calculation broadcasts this callable as
# a value alongside component arrays; it is a scalar object in that operation.
Base.broadcastable(norm::FullEngineRMS)=Ref(norm)
function (norm::FullEngineRMS)(u::AbstractArray,t)
    value=zero(eltype(u))
    @inbounds @fastmath for x in u;value+=abs2(x);end
    Base.FastMath.sqrt_fast(real(value)/norm.length)
end
(norm::FullEngineRMS)(x::Number,t)=abs(x)
(norm::FullEngineRMS)(f,u,t)=sqrt(sum(x->abs2(f(x)),u)/norm.length)

function engine_accept_restricted_weights!(weights,integrator,coordinates,full_atol)
    mass=engine_weight_mass(weights,engine_expand!(coordinates,integrator.u))
    weights.reference<=mass || error("accepted mass fell below its reference")
    push!(weights.references_used,weights.reference)
    push!(weights.accepted_masses,mass)
    push!(weights.accepted_times,integrator.t)
    weights.dynamic && (weights.reference=.5mass)
    engine_set_mass_atol!(weights,full_atol)
    @inbounds for (j,k) in enumerate(coordinates.indices);integrator.opts.abstol[j]=full_atol[k];end
    nothing
end
