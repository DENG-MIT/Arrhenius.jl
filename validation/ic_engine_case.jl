using Arrhenius
using NPZ
using Test
using TOML
include(joinpath(@__DIR__,"..","example","reactors","ic_engine_solver.jl"))
include(joinpath(@__DIR__,"..","example","reactors","ic_engine_setup.jl"))

struct EngineOutputControlProbe{M}
    model::M
end
function (rhs::EngineOutputControlProbe)(du,u,p,t)
    du[1] = 9e6exp(-t/1e-4)
end
function validate_engine_output_controls()
    # An exact finite heating pulse exercises output limits independently of
    # the engine's chemical mechanism and its much tighter accuracy settings.
    model = (network=(nodes=(cylinder=(initial=(gas=(n_species=0,),),),),),)
    problem = (f=EngineOutputControlProbe(model),u0=[300.0],tspan=(0.0,0.002),p=nothing,
        jac=(J,u,p,t)->fill!(J,0.0),tgrad=(du,u,p,t)->(du[1]=-9e10exp(-t/1e-4)),
        isoutofdomain=(u,p,t)->!all(isfinite,u) || u[1]<=0)
    solution = native_engine_sdirk(problem;reltol=1e-6,abstol=1e-6,
        save_everystep=true,dense=false)
    temperatures = first.(solution.u)
    @testset "native engine temperature and crank-angle output limits" begin
        @test maximum(abs,diff(temperatures)) <= 20.0+1e-8
        @test maximum(diff(solution.t)) <= 1/(360*ENGINE_FREQUENCY)*(1+1e-12)
        @test temperatures[end] ≈ 1200-900exp(-20) atol=0.002
        selected = engine_output_indices(solution.t,temperatures)
        @test maximum(abs,diff(temperatures[selected])) <= 20.0+1e-8
        @test maximum(diff(solution.t[selected])) <= 1/(360*ENGINE_FREQUENCY)*(1+1e-12)
    end
end

function validate_ic_engine_case(directory;run=true)
    reference = npzread(joinpath(directory,"ic_engine.reference.npz"))
    mechanism = joinpath(directory,"dodecane_IG.yaml")
    if run
        result = solve_ic_engine(mechanism;integrator=native_engine_sdirk,times=reference["time"],progress=true)
        native = ic_engine_observables(result)
        native["states"] = reduce(hcat,result.states)
        npzwrite(joinpath(directory,"ic_engine.native.npz"),native)
        write_ic_engine_csv(joinpath(directory,"ic_engine.native.csv"),native,result.gas)
    else
        native = npzread(joinpath(directory,"ic_engine.native.npz"))
    end
    gas = CreateSolution(mechanism)
    t = reference["time"]
    temperature_error = maximum(abs,native["temperature"]-reference["temperature"])
    pressure_error = maximum(abs,native["pressure"]./reference["pressure"].-1)
    volume_error = maximum(abs,native["volume"]-reference["volume"])
    mass_error = maximum(abs,native["mass"]./reference["mass"].-1)
    y_error = maximum(abs,native["Y"]-reference["Y"])
    entropy_error = maximum(abs,native["entropy_mass"]-reference["entropy_mass"])
    volume_identity = maximum(abs,native["volume"]-engine_volume.(t))
    mass_balance = maximum(abs,native["mass_balance"].-native["mass_balance"][1])/native["mass"][1]
    energy_scale = max(maximum(abs,native["integrated_work"]),maximum(abs,native["internal_energy"]),1.0)
    energy_balance = maximum(abs,native["energy_balance"].-native["energy_balance"][1])/energy_scale
    switches = engine_switching_times(last(t))[2:end-1]
    continuous = [all(abs(time-switch)>1e-11 for switch in switches) for time in t]
    rate_errors = Dict(key=>maximum(abs,native[key][continuous]-reference[key][continuous]) /
        max(maximum(abs,reference[key][continuous]),1e-30)
        for key in ("mdot_in","mdot_out","mdot_fuel","work_rate","heat_release_rate"))
    native_integrals,reference_integrals = ic_engine_integrals(native,gas),ic_engine_integrals(reference,gas)
    integral_errors = Dict(String(key)=>abs(getproperty(native_integrals,key)/getproperty(reference_integrals,key)-1)
                          for key in keys(native_integrals))
    @testset "Cantera 4 complete eight-revolution engine" begin
        @test native["time"] == t
        @test last(t) >= 8/ENGINE_FREQUENCY
        @test native["crank_angle"] ≈ reference["crank_angle"] rtol=0 atol=1e-13
        @test temperature_error < 0.5
        @test pressure_error < 2e-4
        @test volume_error < 2e-10
        @test mass_error < 5e-5
        @test y_error < 1e-4
        @test entropy_error < 0.5
        @test volume_identity < 1e-10
        @test mass_balance < 2e-7
        @test energy_balance < 2e-6
        @test rate_errors["mdot_in"] < 2e-3
        @test rate_errors["mdot_out"] < 2e-3
        @test rate_errors["mdot_fuel"] < 1e-13
        @test rate_errors["work_rate"] < 2e-3
        @test rate_errors["heat_release_rate"] < 5e-3
        @test integral_errors["heat_J"] < 5e-4
        @test integral_errors["work_J"] < 5e-4
        @test integral_errors["efficiency"] < 5e-4
        @test integral_errors["CO_ppm"] < 0.01
        @test minimum(native["Y"]) >= -1e-13
    end
    validate_engine_output_controls()
    metrics = Dict("points"=>length(t),"temperature_error_K"=>temperature_error,
        "pressure_relative_error"=>pressure_error,"volume_error_m3"=>volume_error,
        "mass_relative_error"=>mass_error,"mass_fraction_error"=>y_error,
        "entropy_error_J_per_kg_K"=>entropy_error,"volume_identity_error_m3"=>volume_identity,
        "minimum_mass_fraction"=>minimum(native["Y"]),
        "native_reltol"=>1e-13,"native_species_atol_kg"=>1e-26,
        "native_temperature_atol_K"=>1e-10,"native_volume_atol_m3"=>1e-20,
        "mass_balance_relative_drift"=>mass_balance,"energy_balance_relative_drift"=>energy_balance,
        "relative_rate_errors"=>rate_errors,"relative_integral_errors"=>integral_errors,
        "native_integrals"=>Dict(String(k)=>v for (k,v) in pairs(native_integrals)),
        "reference_integrals"=>Dict(String(k)=>v for (k,v) in pairs(reference_integrals)))
    open(joinpath(directory,"ic_engine_validation.toml"),"w") do io
        TOML.print(io,metrics;sorted=true)
    end
    println(metrics)
    return metrics
end

function validate_ic_engine_standalone(csv_path;accepted_reference=nothing)
    summary = TOML.parsefile(csv_path*".toml")
    @testset "native engine standalone trajectory and quadrature" begin
        @test summary["end_time_s"] == 8/ENGINE_FREQUENCY
        @test summary["integration_points"] > summary["display_points"] > 2880
        @test summary["maximum_output_interval_s"] <= 1/(360*ENGINE_FREQUENCY)*(1+1e-12)
        @test summary["maximum_output_temperature_change_K"] <= 20.0+1e-8
        @test summary["volume_identity_error_m3"] < 1e-10
        @test summary["mass_balance_relative_drift"] < 2e-7
        @test summary["energy_balance_relative_drift"] < 2e-6
        @test summary["minimum_mass_fraction"] >= -1e-13
        @test maximum(values(summary["relative_quadrature_change"])) < 1e-4
        @test summary["work_quadrature_ledger_relative_error"] < 1e-4
    end
    if accepted_reference !== nothing
        reference = npzread(joinpath(accepted_reference,"ic_engine.reference.npz"))
        gas = CreateSolution(joinpath(accepted_reference,"dodecane_IG.yaml"))
        reference_integrals = engine_integral_values(engine_quadrature_terms(reference,gas;omit_initial=false))
        errors = Dict(String(k)=>abs(summary["integrals"][String(k)]/v-1) for (k,v) in pairs(reference_integrals))
        @testset "Cantera 4 accepted-step engine integrals" begin
            @test last(reference["time"]) == 8/ENGINE_FREQUENCY
            @test errors["heat_J"] < 1e-4
            @test errors["work_J"] < 1e-5
            @test errors["efficiency"] < 1e-4
            @test errors["CO_ppm"] < 1e-4
        end
        summary["accepted_reference_relative_errors"] = errors
        summary["accepted_reference_integrals"] = Dict(String(k)=>v for (k,v) in pairs(reference_integrals))
        open(csv_path*".validation.toml","w") do io
            TOML.print(io,summary;sorted=true)
        end
    end
    return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("provide the Cantera 4 engine reference directory")
    validate_ic_engine_case(ARGS[1];run=!("--cached" in ARGS))
    if "--standalone" in ARGS
        csv_path = ARGS[findfirst(==("--standalone"),ARGS)+1]
        accepted_reference = "--accepted-reference" in ARGS ? ARGS[findfirst(==("--accepted-reference"),ARGS)+1] : nothing
        validate_ic_engine_standalone(csv_path;accepted_reference)
    end
end
