"""Fixed-grid pseudo-arclength predictor-corrector for counterflow flames.
Supported parameters: fuel_mass_flux, oxidizer_mass_flux, pressure, inlet_flux_scale (proportional scaling of both inlet mass fluxes), or control_temperature for a controlled diffusion flame.
The bordered correction retains the banded flame Jacobian and uses a scalar
Schur complement. The initial flame must be converged; no grid refinement is
performed during tracing.
"""
struct CounterflowContinuationPoint
    parameter::Float64
    max_temperature::Float64
    max_strain_rate::Float64
    state::Matrix{Float64}
    residual_inf::Float64
    arclength::Float64
    turning_point::Bool
    tangent_state::Matrix{Float64}
    tangent_parameter::Float64
    next_step::Float64
    minimum_step::Float64
    maximum_step::Float64
    parameter_scale::Float64
    state_scale::Float64
    pressure::Float64
    fuel_mass_flux::Float64
    oxidizer_mass_flux::Float64
    inlet_flux_base::Tuple{Float64,Float64}
end

"""One fixed-grid solution produced by Cantera-style two-point control."""
struct TwoPointControlPoint
    control_temperature::Float64
    temperature_decrement::Float64
    max_temperature::Float64
    max_strain_rate::Float64
    fuel_mass_flux::Float64
    oxidizer_mass_flux::Float64
    state::Matrix{Float64}
    solve_seconds::Float64
    left_control_index::Int
    right_control_index::Int
    control_points::NTuple{4,Float64}
end

"""Trace a counterflow diffusion flame with Cantera's two-point temperature schedule.

At every step, choose control locations from the first/last crossing of the
current relative temperature, lower the selected temperatures by the current
decrement, and solve without grid refinement. The initial schedule and
decrement adaptation follow Cantera's diffusion_flame_continuation example.
"""
function two_point_control!(f::CounterflowDiffusionFlame; max_steps=1000,
        initial_temperature_decrement=20., max_temperature_decrement=100.,
        target_delta_Tmax=20., initial_control_spacing=.6,
        unstable_control_spacing=.95, strain_rate_tol=.1,
        temperature_decrement_schedule=nothing,
        max_consecutive_errors=3, max_time_steps=100, initial_time_step=1.0e-5, loglevel=0)
    (f.converged && isnothing(f.control_points)) ||
        throw(ArgumentError("two-point continuation requires a converged flame with control disabled"))
    max_steps>=1 && max_consecutive_errors>=1 ||
        throw(ArgumentError("max_steps and max_consecutive_errors must be positive"))
    all(x->isfinite(x)&&x>0,(initial_temperature_decrement,max_temperature_decrement,
        target_delta_Tmax,strain_rate_tol)) || throw(ArgumentError("continuation schedule values must be positive and finite"))
    0<initial_control_spacing<unstable_control_spacing<1 ||
        throw(ArgumentError("require 0 < initial_control_spacing < unstable_control_spacing < 1"))

    schedule=isnothing(temperature_decrement_schedule) ? nothing : Float64.(temperature_decrement_schedule)
    isnothing(schedule) || all(x->isfinite(x)&&x>0,schedule) ||
        throw(ArgumentError("temperature_decrement_schedule must contain positive finite values"))
    curve=TwoPointControlPoint[]
    step=min(Float64(initial_temperature_decrement),Float64(max_temperature_decrement))
    Tmax_old=maximum(temperature(f)); amax=_two_point_max_strain_rate(f,f.state)
    strain=amax; errors=0
    for i in 1:max_steps
        !isnothing(schedule) && i>length(schedule) && break
        !isnothing(schedule) && (step=schedule[i])
        spacing=strain>.98amax ? initial_control_spacing : unstable_control_spacing
        target=minimum(temperature(f))+spacing*(maximum(temperature(f))-minimum(temperature(f)))
        saved_state=copy(f.state); saved_control=f.control_points
        saved_mf=f.fuel_mass_flux; saved_mo=f.oxidizer_mass_flux; saved_converged=f.converged
        try
            set_two_point_control!(f;temperature=target,decrement=step)
            ctrl=f.control_points
            left_index=findfirst(==(ctrl[1]),f.grid)
            right_index=findfirst(==(ctrl[3]),f.grid)
            (!isnothing(left_index) && !isnothing(right_index)) ||
                error("two-point control coordinates must lie on the fixed grid")
            if ctrl[2]<f.fuel_temperature+100 || ctrl[4]<f.oxidizer_temperature+100
                f.state=saved_state; f.control_points=saved_control
                f.fuel_mass_flux=saved_mf; f.oxidizer_mass_flux=saved_mo; f.converged=saved_converged
                break
            end
            t0=time_ns()
            solve!(f;refine_grid=false,auto=false,max_time_steps,initial_time_step,loglevel)
            elapsed=(time_ns()-t0)/1e9
            Tmax=maximum(temperature(f)); strain=_two_point_max_strain_rate(f,f.state)
            push!(curve,TwoPointControlPoint(target,step,Tmax,strain,
                f.fuel_mass_flux,f.oxidizer_mass_flux,copy(f.state),elapsed,
                left_index,right_index,ctrl))
            errors=0
            delta=abs(Tmax-Tmax_old)
            if isnothing(schedule)
                if delta<.8target_delta_Tmax
                    step=min(step+3.,max_temperature_decrement)
                elseif delta>target_delta_Tmax
                    step*=.9target_delta_Tmax/delta
                end
            end
            amax=max(strain,amax); Tmax_old=Tmax
            loglevel>0 && println("two-point step $i: Tmax=$Tmax K strain=$strain 1/s dT=$step K")
        catch err
            f.state=saved_state; f.control_points=saved_control
            f.fuel_mass_flux=saved_mf; f.oxidizer_mass_flux=saved_mo; f.converged=saved_converged
            strain/amax<strain_rate_tol && break
            if !isnothing(schedule)
                loglevel>0 && println("two-point scheduled step $i failed at target=$target K dT=$step K: ",sprint(showerror,err))
                rethrow(err)
            end
            step*=.7; errors+=1
            loglevel>0 && println("two-point step $i failed: ",sprint(showerror,err),"; retry dT=$step K")
            errors>=max_consecutive_errors && break
        end
    end
    return curve
end

function _two_point_max_strain_rate(f::CounterflowDiffusionFlame,u)
    n=f.gas.n_species; N=length(f.grid)
    velocity=Vector{Float64}(undef,N)
    @inbounds for j in 1:N
        sumY_over_MW=0.
        for k in 1:n
            sumY_over_MW+=u[k+1,j]/f.gas.MW[k]
        end
        velocity[j]=u[end,j]*R*(1000u[1,j])*sumY_over_MW/f.pressure
    end
    amax=0.
    @inbounds for j in 1:N
        if j==1
            rate=(velocity[2]-velocity[1])/(f.grid[2]-f.grid[1])
        elseif j==N
            rate=(velocity[N]-velocity[N-1])/(f.grid[N]-f.grid[N-1])
        else
            rate=(velocity[j+1]-velocity[j-1])/(f.grid[j+1]-f.grid[j-1])
        end
        amax=max(amax,abs(rate))
    end
    return amax
end

function _continuation_strain_rate(f,u)
    n=f.gas.n_species
    v=[u[end,j]*R*(1000u[1,j])*sum(u[k+1,j]/f.gas.MW[k] for k in 1:n)/f.pressure
        for j in axes(u,2)]
    return maximum(abs.(diff(v)./diff(f.grid)))
end

function _continuation_parameter(f,name::Symbol,fluxbase=nothing)
    name in (:fuel_mass_flux,:oxidizer_mass_flux,:pressure,:inlet_flux_scale,:control_temperature) ||
        throw(ArgumentError("parameter must be :fuel_mass_flux, :oxidizer_mass_flux, :pressure, :inlet_flux_scale, or :control_temperature"))
    if name===:control_temperature
        f isa CounterflowDiffusionFlame && !isnothing(f.control_points) ||
            throw(ArgumentError(":control_temperature requires enabled diffusion two-point control"))
        return Float64((f.control_points[2]+f.control_points[4])/2)
    end
    isnothing(f.control_points) || throw(ArgumentError("disable two-point control before flux or pressure pseudo-arclength continuation"))
    if name===:inlet_flux_scale
        isnothing(fluxbase) && throw(ArgumentError("internal inlet flux baseline is missing"))
        scale=f.fuel_mass_flux/fluxbase[1]
        isapprox(f.oxidizer_mass_flux/fluxbase[2],scale;rtol=1e-10,atol=0) ||
            throw(ArgumentError("inlet mass fluxes do not share a common scale"))
        return Float64(scale)
    end
    Float64(getproperty(f,name))
end
function _continuation_set_parameter!(f,name::Symbol,value,fluxbase=nothing)
    isfinite(value) && value>0 || throw(ArgumentError("continuation parameter must remain positive and finite"))
    if name===:inlet_flux_scale
        isnothing(fluxbase) && throw(ArgumentError("internal inlet flux baseline is missing"))
        f.fuel_mass_flux=fluxbase[1]*value; f.oxidizer_mass_flux=fluxbase[2]*value
    elseif name===:control_temperature
        isnothing(f.control_points) && throw(ArgumentError("two-point control is not enabled"))
        isnothing(fluxbase) && throw(ArgumentError("control-temperature offsets are missing"))
        c=f.control_points
        f.control_points=(c[1],Float64(value)+fluxbase[1],c[3],Float64(value)+fluxbase[2])
    else
        setproperty!(f,name,Float64(value))
    end
    f.converged=false; value
end

function _continuation_border_solve!(f,w,u,p,parameter,r,fp,rhs,z,analytic=true,fluxbase=nothing)
    _continuation_set_parameter!(f,parameter,p,fluxbase)
    counterflow_residual!(r,f,u,w)
    _counterflow_jacobian!(f,u,w,r;analytic)
    bw=2size(u,1)-1
    _,piv=LinearAlgebra.LAPACK.gbtrf!(bw,bw,length(u),w.band)
    copyto!(rhs,vec(r)); rhs .*= -1
    LinearAlgebra.LAPACK.gbtrs!('N',bw,bw,length(u),w.band,piv,rhs)
    dp=sqrt(eps(Float64))*max(abs(p),1.0)
    _continuation_set_parameter!(f,parameter,p+dp,fluxbase)
    counterflow_residual!(fp,f,u,w); fp .= (fp .- r)./dp
    _continuation_set_parameter!(f,parameter,p,fluxbase)
    copyto!(z,vec(fp)); LinearAlgebra.LAPACK.gbtrs!('N',bw,bw,length(u),w.band,piv,z)
    return r
end

function _continuation_correct!(f,w,u0,p0,tu,tp,parameter,weights,qscale,ds;
        tolerance=1e-8,maxiters=18,analytic=true,fluxbase=nothing)
    u=copy(u0); p=p0; B,N=size(u); r=similar(u); fp=similar(u)
    rhs=Vector{Float64}(undef,length(u)); z=similar(rhs)
    pred=u0 .+ ds.*tu; ppred=p0+ds*tp
    cu=tu ./ reshape(weights.^2,B,1); cp=tp/qscale^2
    merit0=Inf
    for it in 1:maxiters
        _continuation_border_solve!(f,w,u,p,parameter,r,fp,rhs,z,analytic,fluxbase)
        g=sum(tu .* (u .- pred) ./ reshape(weights.^2,B,1))+tp*(p-ppred)/qscale^2
        rn=norm(r,Inf)
        rn<=tolerance && abs(g)<=tolerance*max(1.,ds) && return (u,p,rn,it)
        denom=cp-dot(vec(cu),z)
        abs(denom)>eps(Float64)*max(1.,abs(cp)) || return nothing
        dp=(-g-dot(vec(cu),rhs))/denom
        du=reshape(rhs .- z.*dp,B,N); accepted=false; alpha=1.
        for _ in 1:8
            un=u .+ alpha.*du; pn=p+alpha*dp
            if pn>0
                _continuation_set_parameter!(f,parameter,pn,fluxbase); counterflow_residual!(fp,f,un,w)
                gn=sum(tu .* (un .- pred) ./ reshape(weights.^2,B,1))+tp*(pn-ppred)/qscale^2
                merit=hypot(norm(fp,Inf),abs(gn))
                if isfinite(merit) && (merit<merit0 || merit<=(1-1e-4*alpha)*hypot(rn,abs(g)))
                    u=un; p=pn; merit0=merit; accepted=true; break
                end
            end
            alpha*=.5
        end
        accepted || return nothing
    end
    return nothing
end

function _continuation_tangent!(f,w,u,p,parameter,weights,qscale;analytic=true,fluxbase=nothing,direction=1)
    r=similar(u); fp=similar(u); rhs=Vector{Float64}(undef,length(u)); z=similar(rhs)
    _continuation_border_solve!(f,w,u,p,parameter,r,fp,rhs,z,analytic,fluxbase)
    du=direction.*reshape(-z,size(u)); tp=Float64(direction)
    nt=sqrt(sum((du./weights).^2)+(tp/qscale)^2); du./=nt; tp/=nt
    return du,tp
end

_continuation_base(f) = f isa CounterflowDiffusionFlame ? f : f.flow
_continuation_kind(f::CounterflowDiffusionFlame) = "diffusion"
_continuation_kind(f::CounterflowPremixedFlame) = "opposed_premixed"
_continuation_kind(f::CounterflowTwinPremixedFlame) = "twin_premixed"
_continuation_bytes(s) = collect(codeunits(string(s)))
_continuation_text(a,k) = String(vec(UInt8.(a[k])))
function _continuation_mechanism_signature(f)
    g=f.gas
    bytes2hex(SHA.sha256(repr((g.species_names,g.MW,g.thermo,g.reaction,g.trans);context=:limit=>false)))
end
function _continuation_boundary_signature(f)
    b=_continuation_base(f)
    repr((b.fuel_temperature,b.oxidizer_temperature,b.fuel_Y,b.oxidizer_Y,b.transport_model,
        b.soret_enabled,b.flux_gradient_basis,b.radiation_enabled,b.boundary_emissivities,b.fixed_temperature))
end

"""Save a continuation point, tangent, grid, and next arc step to an NPZ checkpoint."""
function save_continuation_checkpoint(path::AbstractString,f,point::CounterflowContinuationPoint;parameter,overwrite=false)
    ispath(path) && !overwrite && throw(ArgumentError("checkpoint exists; use overwrite=true"))
    param=Symbol(parameter); b=_continuation_base(f)
    param in (:fuel_mass_flux,:oxidizer_mass_flux,:pressure,:inlet_flux_scale) ||
        throw(ArgumentError("unsupported continuation checkpoint parameter"))
    isnothing(b.control_points) || throw(ArgumentError("two-point control checkpoint is unsupported"))
    d=Dict{String,Any}(
        "format_utf8"=>_continuation_bytes("arrhenius-counterflow-continuation-v1"),
        "kind_utf8"=>_continuation_bytes(_continuation_kind(f)),
        "parameter_utf8"=>_continuation_bytes(param),
        "mechanism_signature_utf8"=>_continuation_bytes(_continuation_mechanism_signature(f)),
        "species_names_utf8"=>_continuation_bytes(join(f.gas.species_names,"\n")),
        "boundary_signature_utf8"=>_continuation_bytes(_continuation_boundary_signature(f)),
        "grid"=>copy(b.grid),"state"=>copy(point.state),"tangent_state"=>copy(point.tangent_state),
        "parameter_value"=>[point.parameter],"tangent_parameter"=>[point.tangent_parameter],
        "arclength"=>[point.arclength],"next_step"=>[point.next_step],
        "minimum_step"=>[point.minimum_step],"maximum_step"=>[point.maximum_step],
        "parameter_scale"=>[point.parameter_scale],"state_scale"=>[point.state_scale],
        "pressure"=>[point.pressure],"fuel_mass_flux"=>[point.fuel_mass_flux],
        "oxidizer_mass_flux"=>[point.oxidizer_mass_flux],"inlet_flux_base"=>collect(point.inlet_flux_base))
    npzwrite(path,d); path
end

"""Restore a checkpoint after verifying flame kind, mechanism, and boundary conditions."""
function restore_continuation_checkpoint!(f,path::AbstractString)
    a=npzread(path); b=_continuation_base(f)
    _continuation_text(a,"format_utf8")=="arrhenius-counterflow-continuation-v1" || throw(ArgumentError("unknown continuation checkpoint format"))
    _continuation_text(a,"kind_utf8")==_continuation_kind(f) || throw(ArgumentError("continuation checkpoint flame kind mismatch"))
    _continuation_text(a,"mechanism_signature_utf8")==_continuation_mechanism_signature(f) || throw(ArgumentError("continuation checkpoint mechanism mismatch"))
    _continuation_text(a,"species_names_utf8")==join(f.gas.species_names,"\n") || throw(ArgumentError("continuation checkpoint species order mismatch"))
    _continuation_text(a,"boundary_signature_utf8")==_continuation_boundary_signature(f) || throw(ArgumentError("continuation checkpoint boundary mismatch"))
    isnothing(b.control_points) || throw(ArgumentError("disable two-point control before checkpoint restart"))
    z=vec(Float64.(a["grid"])); u=Matrix{Float64}(a["state"]); tu=Matrix{Float64}(a["tangent_state"]); n=f.gas.n_species
    length(z)>=5 && all(isfinite,z) && all(>(0),diff(z)) || throw(ArgumentError("invalid checkpoint grid"))
    size(u)==(n+4,length(z)) && all(isfinite,u) || throw(DimensionMismatch("invalid checkpoint state"))
    size(tu)==size(u) && all(isfinite,tu) || throw(DimensionMismatch("invalid checkpoint tangent"))
    p=only(a["parameter_value"]); tp=only(a["tangent_parameter"]); arc=only(a["arclength"])
    step=only(a["next_step"]); mind=only(a["minimum_step"]); maxd=only(a["maximum_step"])
    qs=only(a["parameter_scale"]); ss=only(a["state_scale"])
    all(isfinite,(p,tp,arc,step,mind,maxd,qs,ss)) && p>0 && arc>=0 && 0<mind<=step<=maxd && qs>0 && ss>0 ||
        throw(ArgumentError("invalid checkpoint continuation scalars"))
    fb=Tuple(Float64.(a["inlet_flux_base"])); length(fb)==2 && all(x->isfinite(x)&&x>=0,fb) ||
        throw(ArgumentError("invalid inlet-flux baseline"))
    param=Symbol(_continuation_text(a,"parameter_utf8"))
    param in (:fuel_mass_flux,:oxidizer_mass_flux,:pressure,:inlet_flux_scale) || throw(ArgumentError("invalid checkpoint parameter"))
    param===:inlet_flux_scale && !all(>(0),fb) && throw(ArgumentError("positive inlet-flux baseline required for proportional scaling"))
    pressure=only(a["pressure"]); mf=only(a["fuel_mass_flux"]); mo=only(a["oxidizer_mass_flux"])
    (isfinite(pressure)&&pressure>0&&isfinite(mf)&&mf>0&&isfinite(mo)&&mo>=0) || throw(ArgumentError("invalid checkpoint flame parameters"))
    b.grid=z; b.state=u; b.pressure=pressure; b.fuel_mass_flux=mf; b.oxidizer_mass_flux=mo; b.converged=true
    return (parameter=param,parameter_value=Float64(p),tangent_state=tu,tangent_parameter=Float64(tp),
        arclength=Float64(arc),next_step=Float64(step),min_ds=Float64(mind),max_ds=Float64(maxd),
        parameter_scale=Float64(qs),state_scale=Float64(ss),inlet_flux_base=fb)
end

function _continuation_point(f,p,u,rn,arc,turn,tu,tp,step,mind,maxd,qs,ss,fb)
    b=_continuation_base(f)
    strain=(f isa CounterflowDiffusionFlame && !isnothing(f.control_points)) ?
        _two_point_max_strain_rate(f,u) : _continuation_strain_rate(f,u)
    CounterflowContinuationPoint(Float64(p),maximum(temperature(f)),strain,copy(u),
        Float64(rn),Float64(arc),Bool(turn),copy(tu),Float64(tp),Float64(step),Float64(mind),Float64(maxd),Float64(qs),Float64(ss),
        b.pressure,b.fuel_mass_flux,b.oxidizer_mass_flux,fb)
end

"""Trace or restart a fixed-grid branch with a bordered pseudo-arclength Newton corrector."""
function pseudo_arclength!(f;parameter=nothing,restart=nothing,steps=40,ds=nothing,min_ds=nothing,max_ds=nothing,
        tolerance=1e-8,max_corrector_iters=18,parameter_scale=nothing,state_scale=nothing,analytic=true,direction=1)
    (f isa CounterflowDiffusionFlame || f isa CounterflowPremixedFlame || f isa CounterflowTwinPremixedFlame) ||
        throw(ArgumentError("pseudo-arclength is implemented for diffusion, opposed premixed, and twin premixed flames"))
    b=_continuation_base(f); saved=isnothing(restart) ? nothing : restore_continuation_checkpoint!(f,restart)
    parameter=isnothing(parameter) ? (isnothing(saved) ? :fuel_mass_flux : saved.parameter) : Symbol(parameter)
    isnothing(saved) || parameter===saved.parameter || throw(ArgumentError("restart parameter differs from checkpoint"))
    fb=isnothing(saved) ? (b.fuel_mass_flux,b.oxidizer_mass_flux) : saved.inlet_flux_base
    if !isnothing(saved)
        _continuation_set_parameter!(f,parameter,saved.parameter_value,parameter===:inlet_flux_scale ? fb : nothing)
    end
    direction in (-1,1) || throw(ArgumentError("direction must be -1 or 1"))
    steps>=1 || throw(ArgumentError("steps must be positive"))
    tolerance>0 && isfinite(tolerance) || throw(ArgumentError("positive finite tolerance required"))
    (f.converged || !isnothing(saved)) || throw(ArgumentError("initial flame must be converged"))
    p0=isnothing(saved) ? _continuation_parameter(f,parameter,parameter===:inlet_flux_scale ? fb : nothing) : saved.parameter_value
    controlbase=parameter===:control_temperature && isnothing(saved) ?
        (f.control_points[2]-p0,f.control_points[4]-p0) : nothing
    baseparam=parameter===:inlet_flux_scale ? fb : controlbase
    u0=copy(f.state)
    qs=isnothing(parameter_scale) ? (isnothing(saved) ? max(abs(p0),1.) : saved.parameter_scale) : Float64(parameter_scale)
    ss=isnothing(state_scale) ? (isnothing(saved) ? 1. : saved.state_scale) : Float64(state_scale)
    isfinite(qs) && qs>0 || throw(ArgumentError("parameter_scale must be positive and finite"))
    isfinite(ss) && ss>0 || throw(ArgumentError("state_scale must be positive and finite"))
    step=Float64(isnothing(ds) ? (isnothing(saved) ? .05 : saved.next_step) : ds)
    mind=Float64(isnothing(min_ds) ? (isnothing(saved) ? step/16 : saved.min_ds) : min_ds)
    maxd=Float64(isnothing(max_ds) ? (isnothing(saved) ? 2step : saved.max_ds) : max_ds)
    isfinite(step) && 0<mind<=step<=maxd || throw(ArgumentError("require 0<min_ds<=ds<=max_ds"))
    weights=ss .* [max(.05*maximum(abs,@view(u0[k,:])),1e-7) for k in axes(u0,1)]
    w=CounterflowWorkspace(f); r=similar(u0); counterflow_residual!(r,f,u0,w); rn0=norm(r,Inf)
    rn0<=tolerance || throw(ArgumentError("initial/restart residual $rn0 exceeds tolerance $tolerance"))
    if isnothing(saved)
        tu,tp=_continuation_tangent!(f,w,u0,p0,parameter,weights,qs;analytic,fluxbase=parameter in (:inlet_flux_scale,:control_temperature) ? baseparam : nothing,direction)
        arc=0.
    else
        tu,tp,arc=saved.tangent_state,saved.tangent_parameter,saved.arclength
    end
    f.state .= u0; _continuation_set_parameter!(f,parameter,p0,parameter in (:inlet_flux_scale,:control_temperature) ? baseparam : nothing); f.converged=true
    pts=CounterflowContinuationPoint[]; lastsign=sign(tp)
    push!(pts,_continuation_point(f,p0,u0,rn0,arc,false,tu,tp,step,mind,maxd,qs,ss,fb))
    for _ in 1:steps
        corrected=nothing
        while corrected===nothing
            corrected=_continuation_correct!(f,w,u0,p0,tu,tp,parameter,weights,qs,step;
                tolerance,maxiters=max_corrector_iters,analytic,fluxbase=parameter in (:inlet_flux_scale,:control_temperature) ? baseparam : nothing)
            if corrected===nothing
                step*=.5
                if step<mind
                    _continuation_set_parameter!(f,parameter,p0,parameter in (:inlet_flux_scale,:control_temperature) ? baseparam : nothing)
                    f.state .= u0; f.converged=true; error("pseudo-arclength corrector failed at min_ds=$mind")
                end
            end
        end
        u1,p1,rn,nit=corrected; arc+=step
        turn=sign(p1-p0)!=0 && sign(p1-p0)!=lastsign
        du=u1.-u0; dp=p1-p0
        orient=sum(du.*tu./reshape(weights.^2,size(u0,1),1))+dp*tp/qs^2
        orient<0 && (du .*= -1; dp *= -1)
        nt=sqrt(sum((du./weights).^2)+(dp/qs)^2); tu=du./nt; tp=dp/nt
        _continuation_set_parameter!(f,parameter,p1,parameter in (:inlet_flux_scale,:control_temperature) ? baseparam : nothing); f.state .= u1; f.converged=true
        step=min(maxd,step*(nit<=5 ? 1.25 : nit>=12 ? .7 : 1.0))
        push!(pts,_continuation_point(f,p1,u1,rn,arc,turn,tu,tp,step,mind,maxd,qs,ss,fb))
        u0=copy(u1); p0=p1; lastsign=sign(tp)
    end
    pts
end
export CounterflowContinuationPoint, TwoPointControlPoint, pseudo_arclength!, two_point_control!, save_continuation_checkpoint, restore_continuation_checkpoint!
