using Arrhenius

if !isdefined(@__MODULE__, :native_bdf)
    include(joinpath(dirname(dirname(pathof(Arrhenius))),"example","reactors","ode_solver.jl"))
end

_pulse_field(t)=190e-21*exp(-0.5*((t-24e-9)/3e-9)^2)

function _pulse_problem(rhs,u,tspan)
    jac=(J,x,p,t)->reactor_jacobian!(J,x,rhs,t)
    tgrad=(du,x,p,t)->(fill!(du,zero(eltype(du)));nothing)
    return (f=rhs,jac=jac,tgrad=tgrad,u0=copy(u),tspan=tspan,p=nothing)
end

function _pulse_solver_stats(solution)
    values=Pair{Symbol,Any}[]
    for name in propertynames(solution.stats)
        value=getproperty(solution.stats,name)
        value isa Number && push!(values,name=>value)
    end
    push!(values,:retcode=>string(solution.retcode))
    return (;values...)
end

function _pulse_matrix(rows,n)
    output=Matrix{Float64}(undef,length(rows),n)
    @inbounds for i in eachindex(rows)
        output[i,:].=rows[i]
    end
    return output
end

"""
    nanosecond_pulse_discharge(mechanism; reltol=1e-11, abstol=1e-23,
                               progress=nothing)

Run the native closed constant-pressure methane/air nanosecond pulse example.
The Boltzmann distribution is refreshed after each 1 ns chunk, matching the
source example's delayed field update. `progress`, when supplied, is called
once per completed chunk with a small named tuple. No files are written.
"""
function nanosecond_pulse_discharge(m::PlasmaMechanism;
        reltol=1e-11,abstol=1e-23,progress=nothing)
    isfinite(reltol) && reltol>0 || throw(ArgumentError("reltol must be positive and finite"))
    isfinite(abstol) && abstol>0 || throw(ArgumentError("abstol must be positive and finite"))

    state=PlasmaState(m;temperature=300.0,pressure=101325.0,
        mole_fractions=Dict("CH4"=>0.095,"O2"=>0.19,"N2"=>0.715,"e"=>1e-11))
    set_reduced_electric_field!(state,_pulse_field(0.0))
    update_eedf!(state)
    reactor=PlasmaEnergyReactor(state;volume=1.0)
    rhs=reactor_rhs(reactor)
    u=reactor_state(reactor)

    times=Float64[]
    temperatures=Float64[]
    densities=Float64[]
    X_rows=Vector{Vector{Float64}}()
    Y_rows=Vector{Vector{Float64}}()
    u_rows=Vector{Vector{Float64}}()
    cached_rows=Vector{Vector{Float64}}()
    chunk_indices=Int[]
    snapshots=Any[]
    chunk_stats=Any[]

    actual_time=0.0
    nominal_time=0.0
    chunks=0
    failure=nothing

    function snapshot(label,nominal)
        p=reactor_properties(rhs,u)
        row=(label=label,nominal_time=Float64(nominal),actual_time=actual_time,
            T=p.T,P=p.P,rho=p.rho,Te=p.Te,E=rhs.electric_field,
            EN=p.reduced_electric_field,mu=rhs.mobility,h_mass=p.h_mass,
            mass=p.mass,volume=p.volume,Y=copy(p.Y),X=copy(p.X),u=copy(u),
            eedf=copy(rhs.eedf.edge_eedf),energy_levels=copy(m.energy_levels))
        push!(snapshots,row)
        return row
    end

    initial=snapshot("initial",0.0)
    terminal=nothing
    try
        while nominal_time<90e-9
            chunk_end=min(nominal_time+1e-9,90e-9)
            requested=Float64[]
            next_time=actual_time
            while next_time<chunk_end
                next_time+=1e-10
                push!(requested,next_time)
            end

            stats=(retcode="empty",)
            if !isempty(requested)
                problem=_pulse_problem(rhs,u,(actual_time,requested[end]))
                solution=native_bdf(problem;linear_solver=:dense,reltol=reltol,
                    abstol=abstol,saveat=requested,save_start=false,
                    save_everystep=false)
                length(solution.t)==length(requested) ||
                    error("pulse solver did not return every requested output time")
                for i in eachindex(requested)
                    ti=Float64(solution.t[i])
                    ui=Vector{Float64}(solution.u[i])
                    p=reactor_properties(rhs,ui)
                    push!(times,ti)
                    push!(temperatures,p.T)
                    push!(densities,p.rho)
                    push!(X_rows,copy(p.X))
                    push!(Y_rows,copy(p.Y))
                    push!(u_rows,ui)
                    push!(chunk_indices,chunks)
                    push!(cached_rows,[p.Te,p.P,p.mass,p.h_mass,rhs.electric_field,
                        p.reduced_electric_field,rhs.mobility])
                end
                u=Vector{Float64}(solution.u[end])
                actual_time=Float64(solution.t[end])
                stats=_pulse_solver_stats(solution)
            end

            snapshot("chunk-"*lpad(string(chunks),3,'0'),nominal_time)
            push!(chunk_stats,(index=chunks,nominal_start=nominal_time,
                nominal_end=chunk_end,actual_end=actual_time,solver_stats=stats))
            update_eedf!(rhs,u;reduced_field=_pulse_field(nominal_time))
            nominal_time=chunk_end
            chunks+=1
            if progress!==nothing
                p=reactor_properties(rhs,u)
                progress((chunks=chunks,nominal_time=nominal_time,
                    actual_time=actual_time,T=p.T,Xe=p.X[m.electron_index]))
            end
        end
        terminal=snapshot("terminal",nominal_time)
    catch error
        failure=(type=string(typeof(error)),message=sprint(showerror,error),
            actual_time=actual_time,nominal_time=nominal_time,chunks=chunks)
        try
            terminal=snapshot("failed-terminal",nominal_time)
        catch snapshot_error
            terminal=(label="failed-terminal",nominal_time=nominal_time,
                actual_time=actual_time,snapshot_error=sprint(showerror,snapshot_error))
        end
    end

    ns=m.n_species
    completed=failure===nothing && nominal_time==90e-9
    return (completed=completed,error=failure,reltol=Float64(reltol),
        abstol=Float64(abstol),t=times,T=temperatures,
        X=_pulse_matrix(X_rows,ns),Y=_pulse_matrix(Y_rows,ns),rho=densities,
        u=_pulse_matrix(u_rows,ns+2),chunk_index=chunk_indices,
        cached=_pulse_matrix(cached_rows,7),plot_EN=_pulse_field.(times),
        initial=initial,terminal=terminal,snapshots=snapshots,
        chunk_stats=chunk_stats,chunks=chunks)
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    1<=length(ARGS)<=2 ||
        error("usage: nanosecond_pulse_discharge.jl MODEL_YAML [SPECIES_DATA_DIR]")
    data_paths=length(ARGS)==2 ? [ARGS[2]] : String[]
    mechanism=PlasmaMechanism(ARGS[1];data_paths)
    result=nanosecond_pulse_discharge(mechanism)
    if !result.completed || result.error!==nothing
        detail=result.error===nothing ? "calculation did not complete" : result.error.message
        error("nanosecond pulse discharge failed: $detail")
    end
    isempty(result.t) && error("nanosecond pulse discharge returned no output states")
    electron=mechanism.electron_index
    println((saved_points=length(result.t),time_s=result.t[end],
        temperature_K=result.T[end],electron_mole_fraction=result.X[end,electron]))
end
