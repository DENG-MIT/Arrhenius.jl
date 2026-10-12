using Arrhenius, LinearAlgebra, SHA, TOML, Libdl
import Arrhenius.NPZ
include(joinpath(@__DIR__,"..","example","reactors","network_solver.jl"))
include(joinpath(@__DIR__,"..","example","reactors","moving_wall_setup.jl"))
include(joinpath(@__DIR__,"numerical_threads.jl"))

const CUSTOM2_SCALARS=("time","temperature","pressure","volume","mass","density","velocity",
    "internal_energy_mass","enthalpy_mass","entropy_mass","gibbs_mass","cp_mass","cv_mass",
    "mean_molecular_weight","gas_internal_energy","wall_kinetic_energy","mechanical_energy",
    "environment_pressure","pressure_work_rate","energy_input","pressure_work_output",
    "mass_balance","energy_balance")

"Pack already observed rows; shared by the real exporter and a no-chemistry fixture."
function custom2_pack(times,states,species,elements,molecular_weights,rows)
    length(times)==length(rows)==size(states,2) || throw(DimensionMismatch("custom2 row count"))
    size(states,1)==16 && length(species)==10 || throw(DimensionMismatch("custom2 state/species count"))
    length(molecular_weights)==10 || throw(DimensionMismatch("custom2 molecular weights"))
    data=Dict{String,Any}("state_full"=>states,"species_names_utf8"=>collect(codeunits(join(species,'\n'))),
        "element_names_utf8"=>collect(codeunits(join(elements,'\n'))),"molecular_weights"=>molecular_weights)
    for key in CUSTOM2_SCALARS
        data[key]=[getproperty(row,Symbol(key)) for row in rows]
    end
    data["time"]=collect(times)
    for key in ("Y","X","species_mass","elements")
        data[key]=permutedims(reduce(hcat,[getproperty(row,Symbol(key)) for row in rows]))
    end
    all(size(data[k])==(length(times),10) for k in ("Y","X","species_mass")) || throw(DimensionMismatch("custom2 species rows"))
    size(data["elements"])==(length(times),length(elements)) || throw(DimensionMismatch("custom2 element rows"))
    return data
end

function custom2_arrays(result)
    model=result.model;gas=model.network.nodes.reactor.initial.gas;rhs=moving_wall_rhs(model)
    rows=map(zip(result.solution.t,result.solution.u)) do (t,u)
        d=moving_wall_diagnostics(rhs,u,t);s=d.nodes.reactor;Y=s.mass_fractions;X=Y2X(gas,Y)
        T,P,V=s.temperature,s.pressure,s.volume;ambient=d.nodes.environment.pressure
        (;time=t,temperature=T,pressure=P,volume=V,mass=s.mass,density=s.density,velocity=only(d.wall_velocities),
            internal_energy_mass=s.internal_energy,enthalpy_mass=s.enthalpy,entropy_mass=cal_smass_mean(gas,T,P,X),
            gibbs_mass=cal_gmass_mean(gas,T,P,X),cp_mass=cal_cpmass_mean(gas,T,P,X),cv_mass=cal_cvmass_mean(gas,T,P,X),
            mean_molecular_weight=sum(X.*gas.MW),gas_internal_energy=d.total_internal_energy,
            wall_kinetic_energy=d.wall_kinetic_energy,mechanical_energy=d.total_internal_energy+d.wall_kinetic_energy+ambient*V,
            environment_pressure=ambient,pressure_work_rate=d.pressure_work_rate,
            energy_input=d.energy_input,pressure_work_output=d.pressure_work_output,
            mass_balance=d.mass_balance,energy_balance=d.energy_balance,Y,X,species_mass=s.mass.*Y,
            elements=[d.element_inventories[e] for e in gas.elements])
    end
    return custom2_pack(result.solution.t,reduce(hcat,result.solution.u),gas.species_names,gas.elements,gas.MW,rows)
end

function custom2_write(output,data,record)
    mkpath(output)
    NPZ.npzwrite(joinpath(output,"native.npz"),data)
    record["native_sha256"]=bytes2hex(sha256(read(joinpath(output,"native.npz"))))
    open(joinpath(output,"native.toml"),"w") do io;TOML.print(io,record;sorted=true);end
end

function custom2_inputs(input)
    paths=[joinpath(input,name) for name in ("h2o2.yaml","h2o2.yaml.npz")]
    append!(paths,[@__FILE__,joinpath(@__DIR__,"numerical_threads.jl")])
    append!(paths,[joinpath(@__DIR__,"..","example","reactors",name) for name in ("custom2.jl","network_solver.jl","moving_wall_setup.jl")])
    for (folder,_,files) in walkdir(joinpath(pkgdir(Arrhenius),"src")),name in files
        endswith(name,".jl") && push!(paths,joinpath(folder,name))
    end
    project=Base.active_project()
    isnothing(project) || append!(paths,[project,joinpath(dirname(project),"Manifest.toml")])
    return Dict(abspath(path)=>bytes2hex(sha256(read(path))) for path in paths)
end

function custom2_runtime(;enforce)
    settings=benchmark_julia_thread_settings(;enforce)
    libraries=Dict(basename(path)=>bytes2hex(sha256(read(path))) for path in Libdl.dllist()
        if isfile(path) && occursin(r"blas|mkl|accelerate|iomp|libgomp"i,basename(path)))
    any(path->occursin(r"cantera|coolprop|python"i,path),Libdl.dllist()) && error("foreign solver loaded")
    return Dict("settings"=>settings,"numerical_library_sha256"=>libraries,"blas_config"=>string(BLAS.get_config()))
end

function custom2_main(input,output)
    ispath(output) && error("use a new custom2 output directory")
    mkpath(output)
    before=custom2_inputs(input);runtime=custom2_runtime(;enforce=true)
    record=Dict{String,Any}("completed"=>false,"inputs_before"=>before,"runtime_before"=>runtime,
        "julia_version"=>string(VERSION),"solver"=>"unchanged public native_network_bdf / QNDF",
        "reltol"=>1e-10,"physical_abstol"=>1e-18,"mass_ledger_abstol_kg"=>1e-15,
        "energy_work_ledger_abstol_J"=>1e-5,"initial_dt_s"=>1e-10,
        "scope"=>"Full custom2 public helper; 101 output points through 0.5 s; output-grid states, not all accepted steps.")
    open(joinpath(output,"started.toml"),"w") do io;TOML.print(io,record;sorted=true);end
    try
        # Use the existing helper and its default integrator settings unchanged.
        result=run_moving_wall_example("custom2",input,joinpath(output,"native.csv");integrator=native_network_bdf)
        data=custom2_arrays(result)
        record["inputs_after"]=custom2_inputs(input)
        record["runtime_after"]=custom2_runtime(;enforce=false)
        record["inputs_unchanged"]=record["inputs_after"]==before
        record["completed"]=true
        custom2_write(output,data,record)
        record["inputs_unchanged"] || error("custom2 sources or inputs changed")
    catch err
        record["failure"]=sprint(showerror,err,catch_backtrace())
        open(joinpath(output,"failure.toml"),"w") do io;TOML.print(io,record;sorted=true);end
        rethrow()
    end
    println("Complete custom2 public helper history saved.")
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("usage: custom2_case.jl INPUT_DIR OUTPUT_DIR")
    custom2_main(abspath(ARGS[1]),abspath(ARGS[2]))
end
