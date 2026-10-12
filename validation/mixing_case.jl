using Arrhenius, LinearAlgebra, SHA, TOML
import Arrhenius.NPZ
include(joinpath(@__DIR__,"..","example","reactors","mixing_solver.jl"))
include(joinpath(@__DIR__,"numerical_threads.jl"))

function main(input,output)
    mkpath(output)
    paths=[joinpath(input,name) for name in ("gri30.yaml","gri30.yaml.npz","air.yaml","air.yaml.npz")]
    append!(paths,[joinpath(@__DIR__,"..","example","reactors",name) for name in ("mixing_solver.jl","network_cases.jl")])
    append!(paths,[@__FILE__,joinpath(@__DIR__,"numerical_threads.jl")])
    for (folder,_,files) in walkdir(joinpath(pkgdir(Arrhenius),"src")),name in files
        endswith(name,".jl") && push!(paths,joinpath(folder,name))
    end
    filehash(path)=bytes2hex(sha256(read(path)))
    before=Dict(path=>filehash(path) for path in paths)
    runtime_before=benchmark_julia_thread_settings(;enforce=true)
    gas=CreateSolution(joinpath(input,"gri30.yaml"))
    air=CreateSolution(joinpath(input,"air.yaml"))
    result=solve_mixing_network(gas,air;save_history=true)
    network=result.network
    initial=network_state(network)
    rhs=network_rhs(network)
    diagnostic=network_diagnostics(rhs,result.state)
    state=diagnostic.nodes.mixer
    X=Y2X(gas,state.mass_fractions)
    derivative=similar(initial)
    rhs(derivative,result.state,nothing,0.)
    payload=Dict{String,Any}("state"=>result.state,"initial_state"=>initial,
        "Y"=>state.mass_fractions,"X"=>X,"flows"=>diagnostic.mass_flow_rates,
        "temperature"=>[state.temperature],"pressure"=>[state.pressure],
        "density"=>[state.density],"mass"=>[state.mass],
        "enthalpy"=>[state.enthalpy],"internal_energy"=>[state.internal_energy],
        "entropy"=>[cal_smass_mean(gas,state.temperature,state.pressure,X)],
        "gibbs"=>[cal_gmass_mean(gas,state.temperature,state.pressure,X)],
        "cp"=>[cal_cpmass_mean(gas,state.temperature,state.pressure,X)],
        "cv"=>[cal_cvmass_mean(gas,state.temperature,state.pressure,X)],
        "mean_molecular_weight"=>[sum(X.*gas.MW)],
        "chemical_potentials_RT"=>cal_g(gas,state.temperature,state.pressure,X)./(R*state.temperature),
        "species_rate"=>derivative[1:end-1],"temperature_rate"=>[derivative[end]],
        "total_energy_rate"=>[diagnostic.total_energy_rate],
        "external_energy_rate"=>[diagnostic.external_energy_rate],
        "element_rates"=>collect(values(diagnostic.element_rates)))
    payload["iterates"]=result.states
    NPZ.npzwrite(joinpath(output,"native.npz"),payload)
    after=Dict(path=>filehash(path) for path in paths)
    runtime_after=benchmark_julia_thread_settings(;enforce=false)
    report=Dict("inputs_before"=>before,"inputs_after"=>after,"inputs_unchanged"=>before==after,
        "runtime_before"=>runtime_before,"runtime_after"=>runtime_after,"julia_version"=>string(VERSION),
        "native_sha256"=>filehash(joinpath(output,"native.npz")),"converged"=>result.converged,"iterations"=>result.iterations,
        "fixed_scale_residual_per_s"=>result.residual,"message"=>result.message,
        "history"=>[Dict(string(k)=>getfield(row,k) for k in fieldnames(typeof(row))) for row in result.history],
        "current_state_residual_per_s"=>max(maximum(abs,derivative[1:end-1])/state.mass,abs(derivative[end])/state.temperature),
        "minimum_Y"=>minimum(state.mass_fractions),"chemistry_enabled"=>network.nodes.mixer.chemistry,
        "species_counts"=>[gas.n_species,air.n_species],"blas_threads"=>BLAS.get_num_threads(),
        "candidate_sha256"=>filehash(joinpath(@__DIR__,"..","example","reactors","mixing_solver.jl")))
    open(joinpath(output,"native.toml"),"w") do io
        TOML.print(io,report;sorted=true)
    end
    println("Mixer: ",result.iterations," iterations; residual ",report["current_state_residual_per_s"]," s^-1")
    before==after || error("source or inputs changed during calculation")
    result.converged || error("native stationary solve did not converge; outputs preserved")
end
length(ARGS)==2 || error("usage: mixing_case.jl INPUT_DIR OUTPUT_DIR")
main(ARGS[1],ARGS[2])
