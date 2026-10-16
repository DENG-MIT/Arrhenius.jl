# Shared engine/preflight solver entry. Qualified calls avoid the distinct
# solve! bindings exported by Arrhenius and SciML. Diagnostics are written
# before propagating an initialization or integration exception.
function engine_probe_solve(ode,algorithm,trace,diagnostics,name,reconstruct;time_offset=0.,kwargs...)
    integrator=nothing
    try
        integrator=SciMLBase.init(ode,algorithm;kwargs...)
        SciMLBase.solve!(integrator)
    catch err
        trace["exception_utf8"]=collect(codeunits(sprint(showerror,err)))
        if integrator!==nothing
            trace["current_time"]=[time_offset+integrator.t]
            trace["current_state"]=reconstruct(integrator.u)
            trace["saved_times"]=time_offset.+integrator.sol.t
            if !isempty(integrator.sol.u)
                trace["saved_states"]=reduce(hcat,[reconstruct(u) for u in integrator.sol.u])
            end
            for field in (:tmp,:u0,:u1)
                if hasproperty(integrator.cache,field)
                    value=getproperty(integrator.cache,field)
                    value isa Vector{Float64} && (trace["solver_cache_"*string(field)]=copy(value))
                end
            end
        end
        engine_failure_snapshot!(diagnostics,name,trace)
        rethrow()
    end
end
