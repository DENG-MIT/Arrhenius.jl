# Focused timing-harness contracts. Imports caller packages, never calls an ODE.
using Test
include(joinpath(@__DIR__,"ic_engine_timing.jl"))

function timing_fixture()
    records=[Dict{String,Any}("index"=>i,"integrated_coordinate_count"=>i<=3 ? 8 : 105,
        "rhs_evaluations"=>3i,"seconds_including_first_specialization"=>.1i) for i in 1:25]
    summary=Dict{String,Any}("end_time_s"=>.16,"integration_points"=>3001,"display_points"=>2881,
        "maximum_output_interval_s"=>1/18000,"maximum_output_temperature_change_K"=>20.,
        "volume_identity_error_m3"=>0.,"mass_balance_relative_drift"=>0.,"energy_balance_relative_drift"=>0.,
        "minimum_mass_fraction"=>0.,"relative_quadrature_change"=>Dict("heat"=>0.),
        "work_quadrature_ledger_relative_error"=>0.)
    checks=Dict{String,Any}("checks_pass"=>true,"quadrature_terms"=>ones(4),"coarse_quadrature_terms"=>ones(4),"ledger_work"=>1.)
    result=(;gas=(;species_names=["X"]),states=[[1.,300.],[.9,400.]],times=[0.,.16],
        integrals=(;heat=1.),quadrature_integrals=(;heat=1.),coarse_integrals=(;heat=1.),integration_points=3001)
    (;records,summary,checks,result,accepted=Dict{String,Any}("segment_01_state"=>ones(2,2)),
        output=Dict("temperature"=>[300.,400.]),initial_state=[1.,300.],diagnostics=Dict("boundary"=>[1.,0.]))
end

@testset "independently pinned first-call library additions" begin
    mktempdir() do directory
        first=joinpath(directory,"libopenblas.so");second=joinpath(directory,"libmkl_avx2.so.2")
        unlisted=joinpath(directory,"libmkl_unlisted.so")
        write(first,"baseline");write(second,"prior dispatch");write(unlisted,"not allowed")
        receipt=joinpath(directory,"prior.toml")
        proof=Dict("checks_pass"=>true,"runtime_after"=>Dict("loaded_numerical_library_sha256"=>
            Dict(basename(first)=>engine_hash(first),basename(second)=>engine_hash(second))))
        open(io->TOML.print(io,proof),receipt,"w")
        receipt_hash=engine_hash(receipt)
        pins=engine_library_pins(receipt,receipt_hash;mapped=[realpath(first)])
        @test Set(keys(pins))==Set(realpath.([first,second]))
        @test !haskey(pins,realpath(unlisted))
        @test_throws ErrorException engine_library_pins(receipt,"wrong hash";mapped=[realpath(first)])
        @test_throws ErrorException engine_library_pins(receipt,receipt_hash;mapped=realpath.([first,unlisted]))
        write(second,"unexpected new bytes")
        @test_throws ErrorException engine_library_pins(receipt,receipt_hash;mapped=[realpath(first)])
    end
end

@testset "complete engine deterministic replay" begin
    a=timing_fixture();b=deepcopy(a)
    @test engine_require_complete(a)
    @test engine_exact_equal(engine_replay_payload(a),engine_replay_payload(b))
    b.records[1]["seconds_including_first_specialization"]+=1
    @test engine_exact_equal(engine_replay_payload(a),engine_replay_payload(b))
    b.records[1]["rhs_evaluations"]+=1
    @test_throws ErrorException engine_require_replay(engine_replay_payload(b),engine_replay_payload(a),"counter mutation")
    b=deepcopy(a);b.records[1]["other_seconds"]=1.
    @test !engine_exact_equal(engine_replay_payload(a),engine_replay_payload(b))
    b=deepcopy(a);b.diagnostics["boundary"][2]=-0.
    @test !engine_exact_equal(engine_replay_payload(a),engine_replay_payload(b))
    b=deepcopy(a);b.accepted["segment_01_state"][1]+=eps()
    @test !engine_exact_equal(engine_replay_payload(a),engine_replay_payload(b))
    b=deepcopy(a);pop!(b.records)
    @test_throws ErrorException engine_require_complete(b)
    b=deepcopy(a);b.summary["end_time_s"]=.12
    @test_throws ErrorException engine_require_complete(b)
    b=deepcopy(a);b.records[4]["integrated_coordinate_count"]=8
    @test_throws ErrorException engine_require_complete(b)
    mktempdir() do directory
        path=joinpath(directory,"retained.npz")
        payload=engine_capture_completed(a,path)
        altered=deepcopy(payload);altered["checks"]["ledger_work"]=2.
        @test_throws ErrorException engine_require_replay(altered,payload,"after save")
        @test isfile(path) && isfile(path*".replay.jls")
        @test engine_exact_equal(deserialize(path*".replay.jls"),payload)
    end
end

@testset "source and loaded-library membership guards" begin
    mktempdir() do directory
        path=joinpath(directory,"first.so");second=joinpath(directory,"second.so")
        write(path,"first");write(second,"second")
        source=engine_tree_hashes(directory)
        write(joinpath(directory,"added.jl"),"# added")
        @test source!=engine_tree_hashes(directory)
        pins=Dict(path=>engine_hash(path),second=>engine_hash(second))
        snapshot(paths)=Dict("environment"=>Dict("threads"=>"1"),"settings"=>Dict("blas"=>1),
            "mapped"=>paths,"mapped_sha256"=>Dict(p=>engine_hash(p) for p in paths))
        a=snapshot([path]);b=snapshot([path,second])
        @test engine_require_runtime(a,a,pins)
        @test engine_require_runtime(a,b,pins;allow_first_load=true)
        @test_throws ErrorException engine_require_runtime(a,b,pins)
        @test_throws ErrorException engine_require_runtime(b,a,pins;allow_first_load=true)
        @test_throws ErrorException engine_require_runtime(a,b,Dict(path=>pins[path]);allow_first_load=true)
        b=deepcopy(a);b["settings"]["blas"]=2
        @test_throws ErrorException engine_require_runtime(a,b,pins)
        b=deepcopy(a);b["environment"]["threads"]="2"
        @test_throws ErrorException engine_require_runtime(a,b,pins)
        write(path,"changed")
        @test_throws ErrorException engine_require_runtime(a,a,pins)
    end
end
