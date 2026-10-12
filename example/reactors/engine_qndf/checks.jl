# Complete-source summary and accepted-history conservation checks.
function engine_qndf_summary_checks(calculation)
    s = calculation.summary
    return all(isfinite,reduce(vcat,calculation.result.states)) &&
        s["end_time_s"] == 0.16 && s["integration_points"] > s["display_points"] > 2880 &&
        s["maximum_output_interval_s"] <= (1/(360*50))*(1+1e-12) &&
        s["maximum_output_temperature_change_K"] <= 20+1e-8 &&
        s["volume_identity_error_m3"] < 1e-10 && s["mass_balance_relative_drift"] < 2e-7 &&
        s["energy_balance_relative_drift"] < 2e-6 && s["minimum_mass_fraction"] >= -1e-13 &&
        maximum(values(s["relative_quadrature_change"])) < 1e-4 &&
        s["work_quadrature_ledger_relative_error"] < 1e-4
end
function engine_qndf_checks(result,output,summary,adapter,audit)
    data=adapter.saved;gas=result.gas;n=gas.n_species
    initial_mass=sum(audit.initial_state[1:n]);mass_error=0.;energy_error=0.
    terms=zeros(4);coarse_terms=zeros(4)
    length(adapter.records)==25 || error("complete source requires 25 intervals")
    for index in 1:25
        key="segment_"*lpad(index,2,'0')
        mass_error=max(mass_error,maximum(abs,data[key*"_mass_balance"].-data["segment_01_mass_balance"][1])/initial_mass)
        energy_scale=max(maximum(abs,data[key*"_internal_energy"]),maximum(abs,data[key*"_pressure_work"]),1.)
        energy_error=max(energy_error,maximum(abs,data[key*"_energy_balance"].-data["segment_01_energy_balance"][1])/energy_scale)
        interval=Dict(replace(k,key*"_output_"=>"")=>v for (k,v) in data if startswith(k,key*"_output_"))
        terms.+=engine_quadrature_terms(interval,gas;omit_initial=false,integral=engine_quadratic_integral)
        coarse=unique!(vcat(collect(1:2:length(interval["time"])),length(interval["time"])))
        coarse_terms.+=engine_quadrature_terms(interval,gas;indices=coarse,omit_initial=false,integral=engine_quadratic_integral)
    end
    work=data["segment_25_pressure_work"][end]-data["segment_01_pressure_work"][1]
    convergence=abs.(coarse_terms./terms.-1)
    efficiency_convergence=abs((coarse_terms[2]/coarse_terms[1])/(terms[2]/terms[1])-1)
    co_convergence=abs((coarse_terms[3]/coarse_terms[4])/(terms[3]/terms[4])-1)
    work_error=abs(terms[2]/work-1)
    totals=Dict(k=>sum(r[k] for r in adapter.records) for k in
        ("accepted_steps","rejected_steps","rhs_evaluations","jacobian_evaluations",
         "nonlinear_iterations","nonlinear_failures","mass_guard_rejections"))
    checks=Dict("checks_pass"=>false,"mass_history_relative_drift"=>mass_error,
        "energy_history_relative_drift"=>energy_error,"quadrature_terms"=>terms,
        "coarse_quadrature_terms"=>coarse_terms,"quadrature_relative_change"=>convergence,
        "efficiency_quadrature_change"=>efficiency_convergence,"CO_quadrature_change"=>co_convergence,
        "work_quadrature_ledger_relative_error"=>work_error,"ledger_work"=>work,
        "CO_ppm"=>1e6*terms[3]/terms[4],"solver_totals"=>totals,
        "prescribed_fuel_source_pass"=>audit.prescribed_fuel_source_pass)
    # Save all values before the final gates, including failures.
    engine_emit!(audit.diagnostics,"complete-checks",checks)
    mass_error<2e-7 && energy_error<2e-6 || error("accepted-history conservation failed")
    terms[3]>0 && terms[4]>0 && coarse_terms[3]>0 && coarse_terms[4]>0 || error("nonzero complete-source CO terms required")
    maximum(convergence)<1e-4 || error("accepted-state quadrature did not converge")
    efficiency_convergence<1e-4 && co_convergence<1e-4 || error("integral ratios did not converge")
    work_error<1e-4 || error("work integral disagrees with its ledger")
    engine_qndf_summary_checks((;result,output,summary)) || error("original complete-source summary gates failed")
    engine_pass_checks!(audit.diagnostics,checks)
end
