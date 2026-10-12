struct EngineScaledRHS{F,M}
    f::F
    model::M
    scale::Vector{Float64}
    raw::Vector{Float64}
    draw::Vector{Float64}
end
function (f::EngineScaledRHS)(du,u,p,t)
    @. f.raw=f.scale*u
    f.f(f.draw,f.raw,p,t)
    @. du=f.draw/f.scale
    nothing
end
function engine_scaled_problem(problem,jacobian)
    ns=problem.f.model.network.nodes.cylinder.initial.gas.n_species
    gas=problem.f.model.network.nodes.cylinder.initial.gas
    scale=ones(length(problem.u0));Vref=ENGINE_CLEARANCE
    scale[1:ns].=gas.MW.*Vref
    scale[ns+1]=1000.;scale[ns+2]=Vref
    scale[ns+3]=sum(problem.u0[1:ns]);scale[ns+4:ns+5].=1000.
    f=EngineScaledRHS(problem.f,problem.f.model,scale,zero(problem.u0),zero(problem.u0))
    physical_jacobian=zeros(length(scale),length(scale))
    jac=(J,u,p,t)->begin
        @. f.raw=scale*u
        jacobian(physical_jacobian,f.raw,p,t)
        for j in eachindex(scale),i in eachindex(scale)
            J[i,j]=physical_jacobian[i,j]*scale[j]/scale[i]
        end
        nothing
    end
    tgrad=(du,u,p,t)->begin
        @. f.raw=scale*u
        problem.tgrad(f.draw,f.raw,p,t)
        @. du=f.draw/scale
        nothing
    end
    outside=(u,p,t)->begin
        @. f.raw=scale*u
        problem.isoutofdomain(f.raw,p,t)
    end
    return (;f,jac,tgrad,isoutofdomain=outside,u0=problem.u0./scale,tspan=problem.tspan,p=problem.p),scale
end

