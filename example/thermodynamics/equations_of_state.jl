# Native CO2 counterpart of Cantera's thermo/equations_of_state.py.
# Usage: julia --project=EOS_ENV example/thermodynamics/equations_of_state.jl CO2.yaml CO2.json
# EOS_ENV needs Arrhenius and Clapeyron 0.6.28. The YAML's first phase must be
# CO2-Ideal with its kinetic sidecar; the second phase is CO2-RK. Generate the
# full Helmholtz parameter JSON with mechanism/export_helmholtz.py.
using Arrhenius, Clapeyron

function _co2_reference_first!(properties)
    for row in 1:3
        reference = properties[row,1]
        for column in axes(properties,2)
            properties[row,column] -= reference
        end
    end
    return properties
end

function _co2_ideal_properties(gas,T,pressures)
    gas.species_names == ["CO2"] || throw(ArgumentError("expected a pure CO2 ideal-gas phase"))
    properties = Matrix{Float64}(undef,5,length(pressures))
    density = Vector{Float64}(undef,length(pressures))
    X = [1.0]
    for (i,p) in enumerate(pressures)
        properties[:,i] .= (cal_hmass_mean(gas,T,p,X),cal_umass_mean(gas,T,p,X),
            cal_smass_mean(gas,T,p,X),cal_cpmass_mean(gas,T,p,X),cal_cvmass_mean(gas,T,p,X)) ./ 1000
        density[i] = p*only(gas.MW)/(Arrhenius.R*T)
    end
    return (;properties=_co2_reference_first!(properties),density)
end

function _co2_rk_properties(model,T,pressures)
    model.species_names == ["CO2"] || throw(ArgumentError("expected a pure CO2 RK phase"))
    properties = Matrix{Float64}(undef,5,length(pressures))
    density,pressure = zeros(length(pressures)),zeros(length(pressures))
    X = [1.0]
    for (i,p) in enumerate(pressures)
        state = redlich_kwong_state(model;T,P=p,X)
        properties[:,i] .= (state.h_mass,state.u_mass,state.s_mass,state.cp_mass,state.cv_mass) ./ 1000
        density[i],pressure[i] = state.rho,state.P
    end
    return (;properties=_co2_reference_first!(properties),density,pressure)
end

function _co2_helmholtz_properties(model,T,pressures)
    psat = T < model.properties.Tc ? first(saturation_pressure(model,T)) : NaN
    T >= model.properties.Tc || (isfinite(psat) && psat > 0) ||
        error("CO2 saturation calculation did not converge")
    properties = Matrix{Float64}(undef,5,length(pressures))
    density,pressure = zeros(length(pressures)),zeros(length(pressures))
    phase = Vector{Symbol}(undef,length(pressures))
    mw = Clapeyron.molecular_weight(model) # kg/mol, using the fitted EOS value
    for (i,p) in enumerate(pressures)
        p == psat && throw(ArgumentError("a saturated state needs a vapor fraction"))
        phase[i] = isnan(psat) ? :unknown : p < psat ? :vapour : :liquid
        v = volume(model,p,T;phase=phase[i],threaded=false) # m³/mol
        # Share the native volume root across all five thermodynamic properties.
        properties[:,i] .= (Clapeyron.VT0.enthalpy(model,v,T),
            Clapeyron.VT0.internal_energy(model,v,T),Clapeyron.VT0.entropy(model,v,T),
            Clapeyron.VT0.isobaric_heat_capacity(model,v,T),
            Clapeyron.VT0.isochoric_heat_capacity(model,v,T)) ./ (1000mw)
        density[i],pressure[i] = mw/v,Clapeyron.VT0.pressure(model,v,T)
    end
    return (;properties=_co2_reference_first!(properties),density,pressure,phase,saturation_pressure=psat)
end

"""
    equations_of_state_calculation(ideal, rk, helmholtz; T=300.0,
        pressures=1e5 .* range(1,100;length=1000))

Evaluate the complete CO2 pressure sweep with native ideal-gas, Redlich–Kwong
and Helmholtz equations of state. Each model returns a 5×N property matrix
(h, u, s, cp, cv), in kJ/kg or kJ/kg/K, and density in kg/m³. Enthalpy,
internal energy and entropy are relative to that model's first pressure point.
The RK model selects its largest physical root; the Helmholtz model selects
stable vapor/liquid roots using its own calculated saturation pressure.
"""
function equations_of_state_calculation(ideal,rk,helmholtz;T=300.0,
        pressures=1e5 .* range(1,100;length=1000))
    T = Float64(T)
    p = Float64.(collect(pressures))
    isfinite(T) && T > 0 || throw(ArgumentError("positive finite temperature required"))
    !isempty(p) && all(x -> isfinite(x) && x > 0,p) ||
        throw(ArgumentError("positive finite pressures required"))
    return (;T,pressure=p,ideal=_co2_ideal_properties(ideal,T,p),
        rk=_co2_rk_properties(rk,T,p),helmholtz=_co2_helmholtz_properties(helmholtz,T,p))
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) == 2 || error("supply CO2.yaml and CO2.json")
    ideal = CreateSolution(ARGS[1])
    rk = RedlichKwongThermo(ARGS[1];phase="CO2-RK")
    helmholtz = SingleFluid(read(ARGS[2],String);coolprop_userlocations=false)
    result = equations_of_state_calculation(ideal,rk,helmholtz)
    println("Computed ",length(result.pressure)," pressures with all three equations of state.")
    println("Helmholtz saturation pressure: ",result.helmholtz.saturation_pressure," Pa")
    for model in (:ideal,:rk,:helmholtz)
        println(model," final (h,u,s,cp,cv): ",getproperty(result,model).properties[:,end])
    end
end
