# Shared complete calculation for the ionized free and burner examples.
using Arrhenius

"Run both ionized stages with a loaded mechanism and fresh flame; return seconds, points, profiles."
function run_ion_source_sequence(gas,case;discretization=:conservative,
        after_stage::F=nothing,clock_ns::C=time_ns,loglevel=0) where {F,C}
    case in ("free","burner") || throw(ArgumentError("case must be free or burner"))
    (gas.n_species,gas.n_reactions)==(56,331) && gas.trans.model==:ionized_gas ||
        throw(ArgumentError("supply the source gri30_ion mechanism"))
    seconds=zeros(2);points=zeros(Int,2);snapshots=Dict{String,Any}();f=nothing
    for (i,stage) in enumerate(("frozen","field"))
        started=clock_ns()
        if i==1
            Tin=case=="free" ? 300. : 600.
            X=mole_fractions(gas,"CH4:1,O2:2,N2:7.52")
            if case=="free"
                f=FreeFlame(gas;T=Tin,P=one_atm,X,width=.05,discretization)
            else
                mdot=.15*one_atm*sum(X.*gas.MW)/(Arrhenius.R*Tin)
                f=BurnerFlame(gas;T=Tin,P=one_atm,X,width=.05,mdot,discretization)
            end
        end
        set_electric_field!(f,i==2)
        solve!(f;ratio=3.,slope=.05,curve=.1,loglevel)
        snapshots[stage]=Dict("grid"=>copy(f.grid),"T"=>temperature(f),
            "Y"=>mass_fractions(f),"X"=>mole_fractions(f),"E"=>electric_field(f),
            "velocity"=>velocity(f),"rho"=>density(f),"qdot"=>heat_release_rate(f),
            "inlet_Y"=>copy(f.inlet_Y),"P"=>[f.pressure],"state"=>copy(f.state))
        seconds[i]=(clock_ns()-started)/1e9
        f.converged || error("$case $stage did not converge")
        all(a->all(isfinite,a),values(snapshots[stage])) || error("nonfinite ionized flame output")
        points[i]=length(f.grid)
        after_stage===nothing || after_stage(f,stage)
    end
    return seconds,points,snapshots
end
