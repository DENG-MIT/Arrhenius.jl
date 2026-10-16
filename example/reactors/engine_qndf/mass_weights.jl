# Conservative componentwise mass mapping for source fraction atol=1e-16.
mutable struct EngineMassWeights
    ns::Int
    scale::Vector{Float64}
    dynamic::Bool
    reference::Float64
    references_used::Vector{Float64}
    accepted_masses::Vector{Float64}
    accepted_times::Vector{Float64}
    mass_guard_rejections::Int
end
function EngineMassWeights(ns,scale,initial_mass;dynamic)
    EngineMassWeights(ns,scale,dynamic,0.5initial_mass,Float64[],Float64[],Float64[],0)
end
function engine_weight_mass(w::EngineMassWeights,u)
    mass=0.0
    @inbounds for k in 1:w.ns
        mass+=w.scale[k]*u[k]
    end
    mass
end
function engine_mass_guard!(w::EngineMassWeights,u)
    rejected=engine_weight_mass(w,u)<w.reference
    rejected && (w.mass_guard_rejections+=1)
    rejected
end
function engine_set_mass_atol!(w::EngineMassWeights,atol)
    @inbounds for k in 1:w.ns
        atol[k]=w.reference*1e-16/w.scale[k]
    end
    nothing
end
function engine_accept_mass_weights!(w::EngineMassWeights,integrator)
    mass=engine_weight_mass(w,integrator.u)
    w.reference<=mass || error("accepted mass fell below the error-weight reference")
    push!(w.references_used,w.reference)
    push!(w.accepted_masses,mass)
    push!(w.accepted_times,integrator.t)
    w.dynamic && (w.reference=0.5mass)
    engine_set_mass_atol!(w,integrator.opts.abstol)
    nothing
end
