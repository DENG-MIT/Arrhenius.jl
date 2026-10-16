# Loading and serialization contract for the actual public exporter; no chemistry.
# julia --project=CALLER_ENV validation/ic_engine_qndf_export_contract.jl [OUTPUT_DIR]
using Test, TOML
module EngineCaseExporter end
const exporter_path=joinpath(@__DIR__,"ic_engine_qndf_case.jl")
mktempdir() do unrelated
    cd(unrelated) do
        before=readdir()
        Base.include(EngineCaseExporter,exporter_path)
        @test readdir()==before
    end
end
const E=EngineCaseExporter

function mocked_engine_calculation()
    n=100;state=zeros(105);state[1]=1.;state[101]=300.;state[102]=1e-4
    names=["species_"*string(k) for k in 1:n]
    times=E.NativeEngineQNDF.engine_switching_times(.16)
    accepted=Dict{String,Any}();records=Dict{String,Any}[]
    scalar_names=("time","temperature","pressure","volume","mass","entropy_mass",
        "mdot_in","mdot_out","mdot_fuel","work_rate","heat_release_rate")
    for index in 1:25
        key="segment_"*lpad(index,2,'0');t=times[index:index+1]
        accepted[key*"_time"]=copy(t)
        accepted[key*"_state"]=hcat(state,state)
        for name in scalar_names
            accepted[key*"_output_"*name]=name=="time" ? copy(t) : ones(2)
        end
        y=zeros(2,n);y[:,1].=1
        accepted[key*"_output_Y"]=y
        push!(records,Dict("start"=>t[1],"stop"=>t[2],"full_state_count"=>105,
            "rms_norm_denominator"=>105,"integrated_coordinate_count"=>(index<=3 ? 8 : 105),
            "element_flux_relative_error"=>0.,"element_flux_quadrature_change"=>0.,
            "maximum_reference_to_mass_ratio"=>.5))
    end
    checks=Dict{String,Any}("checks_pass"=>true,"quadrature_terms"=>ones(4),
        "coarse_quadrature_terms"=>ones(4),"ledger_work"=>1.,
        "mass_history_relative_drift"=>0.,"energy_history_relative_drift"=>0.,
        "quadrature_relative_change"=>zeros(4),"efficiency_quadrature_change"=>0.,
        "CO_quadrature_change"=>0.,"work_quadrature_ledger_relative_error"=>0.,
        "prescribed_fuel_source_pass"=>true)
    summary=Dict{String,Any}("end_time_s"=>.16,"minimum_mass_fraction"=>0.,
        "volume_identity_error_m3"=>0.,"maximum_output_interval_s"=>1e-5,
        "maximum_output_temperature_change_K"=>1.,
        "integrals"=>Dict("heat_J"=>1.,"work_J"=>1.,"efficiency"=>1.,"CO_ppm"=>1.))
    result=(;gas=(;species_names=names),states=[copy(state) for _ in times])
    output=Dict{String,Any}("time"=>times,"Y"=>reduce(vcat,[permutedims(state[1:n]) for _ in times]))
    (;accepted,initial_state=state,result,checks,summary,records,output)
end

function exporter_schema_contract(directory)
    mkpath(directory)
    # Arbitrary bytes are sufficient here: only the real fingerprint and archive
    # paths are tested. No mechanism parser or chemistry calculation is called.
    mechanism=joinpath(directory,"mock.yaml")
    write(mechanism,"schema fixture\n");write(mechanism*".npz","sidecar fixture\n")
    calculation=mocked_engine_calculation()
    before=E.engine_case_inputs(mechanism)
    runtime_before=E.engine_case_runtime(;enforce=true)
    after=E.engine_case_inputs(mechanism)
    runtime_after=E.engine_case_runtime(;enforce=false)
    output=joinpath(directory,"mock.npz")
    record=E.write_engine_case(calculation,output,before,after,runtime_before,runtime_after)
    data=E.NPZ.npzread(output)
    @test record["checks_pass"] && record["inputs_unchanged"]
    @test record["inputs_before"]==record["inputs_after"]==before
    @test split(String(data["species_names_utf8"]),'\n')==calculation.result.gas.species_names
    @test size(data["exact_source_initial_state"])==(105,)
    @test size(data["states"])==(105,26)
    @test size(data["quadrature_terms"])==size(data["coarse_quadrature_terms"])==(4,)
    @test size(data["ledger_work"])==size(data["checks_pass"])==(1,)
    for index in 1:25
        key="segment_"*lpad(index,2,'0')
        @test size(data[key*"_time"])==(2,)
        @test size(data[key*"_state"])==(105,2)
        @test size(data[key*"_output_Y"])==(2,100)
        for name in ("time","temperature","pressure","volume","mass","entropy_mass",
                "mdot_in","mdot_out","mdot_fuel","work_rate","heat_release_rate")
            @test size(data[key*"_output_"*name])==(2,)
        end
    end
    # verify_mechanisms consumes these two top-level record keys and compares
    # sidecar source_sha256 to the actual mechanism bytes.
    @test record["mechanism_sha256"]==E.engine_case_hashfile(mechanism)
    @test record["sidecar_sha256"]==E.engine_case_hashfile(mechanism*".npz")
    @test length(before["thread_guard_sha256"])==64
    @test all(haskey(before,key) for key in ("source_sha256","helper_sha256","exporter_sha256",
        "caller_project_sha256","caller_manifest_sha256"))
    @test runtime_before["settings"]["julia_threads"]==runtime_after["settings"]["julia_threads"]==1
    @test runtime_before["settings"]["blas_threads"]==runtime_after["settings"]["blas_threads"]==1
    @test all(value==1 for (key,value) in runtime_after["settings"] if startswith(key,"mkl_threads:"))
    @test runtime_after["loaded_numerical_library_sha256"] isa Dict
    for key in ("mechanism_sha256","sidecar_sha256","source_sha256","helper_sha256",
            "exporter_sha256","thread_guard_sha256","caller_project_sha256","caller_manifest_sha256")
        changed=deepcopy(after);changed[key]=Dict("changed"=>"fixture")
        failed=joinpath(directory,"changed-"*key*".npz")
        @test_throws ErrorException E.write_engine_case(calculation,failed,before,changed,runtime_before,runtime_after)
        @test isfile(failed)
        @test !TOML.parsefile(replace(failed,".npz"=>".toml"))["inputs_unchanged"]
    end
    invalid_runtime=Dict("checks_pass"=>false,"capture_error"=>"mocked lazy backend mismatch")
    failed=joinpath(directory,"backend-failure.npz")
    @test_throws ErrorException E.write_engine_case(calculation,failed,before,after,runtime_before,invalid_runtime)
    @test isfile(failed)
    @test !TOML.parsefile(replace(failed,".npz"=>".toml"))["checks_pass"]
    path=joinpath(directory,"mutated.yaml");write(path,"one");write(path*".npz","sidecar")
    original=E.engine_case_inputs(path);write(path,"two")
    @test original["mechanism_sha256"]!=E.engine_case_inputs(path)["mechanism_sha256"]
    nothing
end

@testset "Public engine exporter schema and provenance without chemistry" begin
    if isempty(ARGS)
        mktempdir(exporter_schema_contract)
    else
        exporter_schema_contract(abspath(ARGS[1]))
    end
end
