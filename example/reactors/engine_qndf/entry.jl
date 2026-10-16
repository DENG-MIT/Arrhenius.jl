module NativeEngineQNDF

using Arrhenius, SciMLBase, OrdinaryDiffEqBDF, ForwardDiff, LinearAlgebra, SparseArrays
import Arrhenius.NPZ
import TOML

include(joinpath(@__DIR__,"..","ic_engine_setup.jl"))
for helper in ("derivatives.jl","signed_reactor.jl",
        "unseeded_derivatives.jl","scaling.jl",
        "mass_weights.jl","element_budget.jl",
        "invariant_coordinates.jl")
    include(joinpath(@__DIR__,helper))
end

mutable struct EngineDiagnosticLog{S}
    sink::S
    snapshots::Dict{String,Any}
end
EngineDiagnosticLog(sink=nothing)=EngineDiagnosticLog(sink,Dict{String,Any}())
function engine_emit!(log::EngineDiagnosticLog,name,data)
    snapshot=deepcopy(data)
    log.snapshots[name]=snapshot
    log.sink===nothing || log.sink(name,snapshot)
    nothing
end

struct EngineFileSink
    directory::String
    function EngineFileSink(directory)
        directory=abspath(directory)
        mkpath(directory)
        new(directory)
    end
end
function (sink::EngineFileSink)(name,data)
    occursin(r"^[A-Za-z0-9_-]+$",name) || throw(ArgumentError("invalid diagnostic name"))
    if all(value->value isa Number || value isa AbstractArray{<:Number},values(data))
        NPZ.npzwrite(joinpath(sink.directory,name*".npz"),data)
    else
        open(joinpath(sink.directory,name*".toml"),"w") do io
            TOML.print(io,data;sorted=true)
        end
    end
    nothing
end

# Preserve the original calculation failure even if its diagnostic sink fails.
function engine_failure_snapshot!(log,name,data)
    try
        engine_emit!(log,name,data)
    catch err
        log.snapshots["diagnostic-sink-exception"]=Dict(
            "exception_utf8"=>collect(codeunits(sprint(showerror,err))))
    end
    nothing
end

function engine_pass_checks!(log,checks)
    checks["checks_pass"]=true
    engine_emit!(log,"complete-checks",checks)
    checks
end

# Explicit reaction-order metadata selects a different Cantera mass-action
# kernel, even when its numbers equal stoichiometry. This adapter accepts the
# flattened source mechanism with elementary integer-order reactions only.
function engine_check_order_metadata(mechanism)
    document=Arrhenius.YAML.load_file(mechanism)
    reactions=get(document,"reactions",nothing)
    reactions isa AbstractVector || throw(ArgumentError("engine mechanism requires inline reactions"))
    for (index,reaction) in enumerate(reactions)
        reaction isa AbstractDict || throw(ArgumentError("engine reaction $index is not a mapping"))
        haskey(reaction,"orders") && throw(ArgumentError("engine reaction $index has explicit orders, which are unsupported by this signed-rate adapter"))
    end
    nothing
end

struct EngineCalculationFailure{E,D} <: Exception
    cause::E
    diagnostics::D
end
Base.showerror(io::IO,err::EngineCalculationFailure)=showerror(io,err.cause)

function engine_cli_failure_guard(calculate,directory)
    try
        calculate()
    catch err
        if err isa EngineCalculationFailure
            try
                log=EngineDiagnosticLog(EngineFileSink(directory))
                for (name,data) in err.diagnostics
                    engine_failure_snapshot!(log,name,data)
                end
                if haskey(log.snapshots,"diagnostic-sink-exception")
                    err.diagnostics["diagnostic-sink-exception"]=log.snapshots["diagnostic-sink-exception"]
                end
            catch sink_error
                err.diagnostics["diagnostic-sink-exception"]=Dict(
                    "exception_utf8"=>collect(codeunits(sprint(showerror,sink_error))))
            end
        end
        rethrow()
    end
end

include(joinpath(@__DIR__,"solve.jl"))

engine_bitwise_equal(a,b)=reinterpret(UInt64,a)==reinterpret(UInt64,b)
mutable struct EngineSourceAudit{D}
    diagnostics::D
    index::Int
    initial_state::Vector{Float64}
    boundary_state::Vector{Float64}
    last_endpoint::Vector{Float64}
    inactive::Vector{Int}
    full_coordinates_active::Bool
    first_fuel_rhs_seen::Bool
    prescribed_fuel_source_pass::Bool
end
EngineSourceAudit(log)=EngineSourceAudit(log,0,Float64[],Float64[],Float64[],Int[],false,false,false)

struct EngineSignedJacobian end
(::EngineSignedJacobian)(p)=engine_ad_unseeded_jacobian(p;rate_evaluator=trial_wdot!,clip_trials=false)
struct EngineCoordinates{A}
    audit::A
end
function (factory::EngineCoordinates)(problem)
    audit=factory.audit;index=audit.index
    c=engine_invariant_coordinates(problem)
    scaled,scale=engine_scaled_problem(problem,problem.jac)
    state=problem.u0
    stops=engine_switching_times(.16);start,stop=stops[index],stops[index+1]
    network=problem.f.model.network
    phase=mod(index-2,6)+1
    inlet=network.flows[1].time_function(0.)
    outlet=network.flows[3].time_function(0.)
    fuel=network.flows[2].mdot
    engine_emit!(audit.diagnostics,"boundary"*string(index),Dict(
        "initial_state"=>copy(state),"previous_endpoint"=>copy(audit.last_endpoint),
        "scale"=>scale,"scale_roundtrip"=>scaled.u0.*scale,"coordinates"=>c.indices,
        "interval"=>[start,stop],"inlet_gate"=>[inlet],"outlet_gate"=>[outlet],"fuel_rate"=>[fuel]))
    length(state)==105 || error("the exact 105-component source mechanism is required")
    engine_bitwise_equal(scaled.u0.*scale,state) || error("non-bitwise scale roundtrip")
    if index==1
        audit.initial_state=copy(state)
        audit.inactive=copy(c.removed)
        source=moving_wall_state(engine_network(network.nodes.cylinder.initial.gas))
        engine_bitwise_equal(state,source) || error("the original source initial state is required")
    else
        engine_bitwise_equal(state,audit.last_endpoint) || error("non-bitwise physical boundary")
    end
    length(c.indices)==(index<=3 ? 8 : 105) || error("source coordinate activation mismatch")
    index<=3 && !all(iszero,state[audit.inactive]) && error("nonzero absent-element species before opening")
    inlet==(phase in (1,6)) || error("source inlet branch mismatch")
    outlet==(phase in (5,6)) || error("source outlet branch mismatch")
    fuel==(phase==3 ? ENGINE_INJECTION_RATE : 0.) || error("source fuel branch mismatch")
    audit.boundary_state=copy(state)
    audit.full_coordinates_active=length(c.indices)==105
    c
end

struct EngineFuelRHS{F,A}
    f::F
    audit::A
    scale::Vector{Float64}
    index::Int
end
function (rhs::EngineFuelRHS)(du,u,p,t)
    audit=rhs.audit
    first=rhs.index==4 && !audit.first_fuel_rhs_seen
    if first
        audit.first_fuel_rhs_seen=true
        physical=u.*rhs.scale
        engine_emit!(audit.diagnostics,"first-fuel-rhs-input",Dict(
            "state"=>physical,"local_time"=>[t],"full_coordinates_active"=>[Int(audit.full_coordinates_active)]))
        audit.full_coordinates_active || error("fuel RHS ran before full activation")
        engine_bitwise_equal(physical,audit.boundary_state) || error("fuel-opening state changed")
        t==0. || error("first fuel RHS does not start at the existing source event")
    end
    rhs.f(du,u,p,t)
    if first
        derivative=copy(rhs.f.draw)
        gas=rhs.f.model.network.nodes.cylinder.initial.gas;n=gas.n_species
        fuel=only(findall(==("c12h26"),gas.species_names))
        rate=rhs.f.model.network.flows[2].mdot
        engine_emit!(audit.diagnostics,"first-fuel-rhs",Dict("state"=>u.*rhs.scale,
            "derivative"=>derivative,"prescribed_fuel_rate"=>[rate],"fuel_index"=>[fuel]))
        derivative[fuel]==rate==ENGINE_INJECTION_RATE || error("incorrect first fuel-species source")
        derivative[n+3]==rate || error("incorrect first fuel mass ledger")
        all(iszero,derivative[setdiff(audit.inactive,[fuel])]) || error("unexpected absent-element source")
        audit.prescribed_fuel_source_pass=true
    end
    nothing
end
struct EngineRHSObserver{A}
    audit::A
end
function (observer::EngineRHSObserver)(f,scale,index,start,stop)
    observer.audit.index=index
    EngineFuelRHS(f,observer.audit,scale,index)
end
function engine_finish_interval!(audit,index,states,record)
    engine_bitwise_equal(first(states),audit.boundary_state) || error("saved interval entry changed")
    index<=3 && !all(u->all(iszero,u[audit.inactive]),states) && error("accepted state left absent-element subspace")
    index>=4 && record["integrated_coordinate_count"]!=105 && error("full coordinates were deactivated")
    index==4 && !(audit.first_fuel_rhs_seen && audit.prescribed_fuel_source_pass) && error("fuel-opening audit incomplete")
    audit.last_endpoint=copy(last(states))
    nothing
end

include(joinpath(@__DIR__,"adapter.jl"))
include(joinpath(@__DIR__,"checks.jl"))

"Solve the eight-revolution n-dodecane engine with caller-owned QNDF integration."
function solve_ic_engine_qndf(mechanism;diagnostic_sink=nothing,progress=false)
    log=EngineDiagnosticLog(diagnostic_sink)
    audit=EngineSourceAudit(log)
    adapter=EngineQNDFAdapter(signed_jacobian_factory=EngineSignedJacobian(),
        rhs_observer=EngineRHSObserver(audit),coordinate_factory=EngineCoordinates(audit),diagnostics=log)
    try
        engine_check_order_metadata(mechanism)
        result=solve_ic_engine(mechanism;integrator=adapter,reltol=1e-12,progress)
        output=ic_engine_observables(result)
        summary=ic_engine_summary(result,output)
        checks=engine_qndf_checks(result,output,summary,adapter,audit)
        return (;result,output,summary,checks,records=adapter.records,accepted=adapter.saved,
            initial_state=copy(audit.initial_state),diagnostics=log.snapshots)
    catch err
        failure=Dict{String,Any}("exception_utf8"=>collect(codeunits(sprint(showerror,err))),
            "recorded_intervals"=>[length(adapter.records)])
        isempty(audit.last_endpoint) || (failure["last_completed_endpoint"]=copy(audit.last_endpoint))
        isempty(audit.initial_state) || (failure["source_initial_state"]=copy(audit.initial_state))
        engine_failure_snapshot!(log,"calculation-exception",failure)
        throw(EngineCalculationFailure(err,log.snapshots))
    end
end

export solve_ic_engine_qndf, EngineFileSink, EngineCalculationFailure
end
