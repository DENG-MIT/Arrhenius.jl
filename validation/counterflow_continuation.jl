# Cantera diffusion_flame_continuation.py at 726522be4e2a13454d8415b7ef799d621f665cf3.
module CounterflowContinuation

using Arrhenius, LinearAlgebra

export calculate

_strain(f) = begin
    z, v = f.grid, velocity(f)
    maximum(j == 1 ? abs((v[2]-v[1])/(z[2]-z[1])) :
        j == length(z) ? abs((v[end]-v[end-1])/(z[end]-z[end-1])) :
        abs((v[j+1]-v[j-1])/(z[j+1]-z[j-1])) for j in eachindex(z))
end

_snapshot(f) = (grid=copy(f.grid), state=copy(f.state), controls=f.control_points,
    fuel=f.fuel_mass_flux, oxidizer=f.oxidizer_mass_flux, converged=f.converged,
    anchor=f.anchor)
function _restore!(f, s)
    f.grid=s.grid; f.state=s.state; f.control_points=s.controls
    f.fuel_mass_flux=s.fuel; f.oxidizer_mass_flux=s.oxidizer
    f.converged=s.converged; f.anchor=s.anchor
    f
end

function _profile(f)
    (grid=copy(f.grid), state=copy(f.state), controls=f.control_points,
        T=copy(temperature(f)), Y=copy(mass_fractions(f)), velocity=copy(velocity(f)))
end

function _integral_and_width(f)
    q=heat_release_rate(f); z=f.grid; T=temperature(f)
    integral=sum(.5*(q[j]+q[j+1])*(z[j+1]-z[j]) for j in 1:length(z)-1)
    mid=.5*(minimum(T)+maximum(T)); hot=findall(>(mid),T)
    width=length(hot)>=2 ? z[last(hot)]-z[first(hot)] : 0.
    integral, width
end

function _row(f, step, decrement, errors, success, strain, amax, spacing)
    integral, width=_integral_and_width(f)
    (step=step, Tmax=maximum(temperature(f)), strain=strain, ratio=strain/amax,
        decrement=decrement, errors=errors, success=success, npoints=length(f.grid),
        fuel=f.fuel_mass_flux, oxidizer=f.oxidizer_mass_flux,
        integrated_heat_release=integral, flame_width=width, spacing=spacing)
end

"""Run the canonical two-point counterflow continuation on a fresh native flame.
`initial_slope`/`initial_curve`/`initial_prune` apply only to the first solve;
`max_points` bounds the grid size of every solve."""
function calculate(gas; slope=.1, curve=.2, prune=0., capture=false, maxsteps=1000,
        initial_slope=slope, initial_curve=curve, initial_prune=prune, max_points=1200)
    f=CounterflowDiffusionFlame(gas;fuel="H2:1",oxidizer="O2:1",mdot_fuel=.5,
        mdot_oxidizer=3.,T_fuel=300.,T_oxidizer=500.,P=1e5,width=.018)
    solve!(f;ratio=4.,slope=initial_slope,curve=initial_curve,prune=initial_prune,max_points,loglevel=0)
    profiles=capture ? NamedTuple[_profile(f)] : NamedTuple[]
    strain=amax=_strain(f); increment=20.; errors=0; source_success=false
    reason="step_cap"; records=NamedTuple[]
    for i in 1:maxsteps
        spacing=strain>.98*amax ? .6 : .95
        T=temperature(f); oldT=maximum(T)
        target=minimum(T)+spacing*(oldT-minimum(T)); backup=_snapshot(f)
        try
            set_two_point_control!(f;temperature=target)
        catch e
            e isa ArgumentError || rethrow()
            reason="control_selection_error"; _restore!(f,backup); break
        end
        c=f.control_points
        if c[2]-increment<f.fuel_temperature+100 || c[4]-increment<f.oxidizer_temperature+100
            reason="inlet_plus_100K"; source_success=true; _restore!(f,backup); break
        end
        set_two_point_control!(f;temperature=target,decrement=increment)
        converged=false
        try
            solve!(f;auto=false,ratio=4.,slope,curve,prune,max_points,max_time_steps=100,loglevel=0)
            converged=true
        catch e
            e isa ErrorException || rethrow()
            source_success=strain/amax<.1
            _restore!(f,backup)
            if source_success
                reason="strain_below_0.1_after_failure"
                break
            end
            increment*=.7; errors+=1
        end
        if converged
            successes=true
            increment_delta=abs(maximum(temperature(f))-oldT)
            if increment_delta<16
                increment=min(increment+3,100.)
            elseif increment_delta>20
                increment*=18/increment_delta
            end
            strain=_strain(f); amax=max(amax,strain)
            errors=0
            if capture; push!(profiles,_profile(f)); end
        else
            successes=false
        end
        push!(records,_row(f,i,increment,errors,successes,strain,amax,spacing))
        if errors>=3
            reason="three_successive_errors"; break
        end
    end
    final_ratio=strain/amax
    (;flame=f, data=records, profiles=profiles, maximum_strain=amax,
        final_ratio=final_ratio, source_success=source_success,
        termination_reason=reason)
end

end
