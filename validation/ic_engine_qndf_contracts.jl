# Run in the caller's SciML environment. This loads the real CLI without
# executing its engine calculation, then checks the common solver and sinks.
using Test, TOML
module EnginePublicEntry end
module EnginePublicCLI end
const reactors=normpath(joinpath(@__DIR__,"..","example","reactors"))
mktempdir() do unrelated
    cd(unrelated) do
        before=readdir()
        Base.include(EnginePublicEntry,joinpath(reactors,"ic_engine_qndf_solver.jl"))
        Base.include(EnginePublicCLI,joinpath(reactors,"ic_engine.jl"))
        @test readdir()==before
    end
end
const E=EnginePublicEntry.NativeEngineQNDF
Base.include(E,joinpath(@__DIR__,"ic_engine_qndf_solver_contract.jl"))

@testset "Engine common solver and diagnostic delivery" begin
    @test E.engine_qndf_solver_contract()
    log=E.EngineDiagnosticLog()
    source=Dict("state"=>[1.,2.])
    E.engine_emit!(log,"state",source)
    source["state"][1]=3.
    @test log.snapshots["state"]["state"]==[1.,2.]
    @test log.sink===nothing
    failing=E.EngineDiagnosticLog((name,data)->error("intentional sink failure"))
    @test_throws ErrorException E.engine_emit!(failing,"preserved",Dict("state"=>[8.]))
    @test failing.snapshots["preserved"]["state"]==[8.]
    @test E.engine_failure_snapshot!(failing,"original-error",Dict("exception_utf8"=>collect(codeunits("original"))))===nothing
    @test String(failing.snapshots["original-error"]["exception_utf8"])=="original"
    @test occursin("intentional sink failure",String(failing.snapshots["diagnostic-sink-exception"]["exception_utf8"]))
    mktempdir() do folder
        sink=E.EngineFileSink(folder)
        @test_throws ArgumentError sink("../escape",Dict("state"=>[1.]))
        @test isempty(readdir(folder))
        sink("safe",Dict("state"=>[2.]))
        @test E.NPZ.npzread(joinpath(folder,"safe.npz"))["state"]==[2.]
        diagnostic=E.EngineDiagnosticLog(sink)
        checks=Dict{String,Any}("checks_pass"=>false,"quadrature_terms"=>[1.,2.,3.,4.],
            "solver_totals"=>Dict("accepted_steps"=>5))
        # Exercise the actual pre-gate sink and post-success publisher used
        # by the full calculation, including non-array check dictionaries.
        E.engine_emit!(diagnostic,"complete-checks",checks)
        failed=TOML.parsefile(joinpath(folder,"complete-checks.toml"))
        @test !failed["checks_pass"]
        @test failed["solver_totals"]["accepted_steps"]==5
        @test failed["quadrature_terms"]==[1.,2.,3.,4.]
        @test E.engine_pass_checks!(diagnostic,checks)===checks
        @test checks["checks_pass"]
        @test diagnostic.snapshots["complete-checks"]["checks_pass"]
        @test TOML.parsefile(joinpath(folder,"complete-checks.toml"))["checks_pass"]
    end
end

@testset "Engine explicit-order eligibility" begin
    mktempdir() do folder
        path=joinpath(folder,"mechanism.yaml")
        write(path,"reactions:\n- equation: A + B => C\n  rate-constant: {A: 1.0, b: 0, Ea: 0}\n")
        @test E.engine_check_order_metadata(path)===nothing
        for orders in ("{A: 1, B: 1}","{}","{A: 0.5}","{A: -1}")
            write(path,"reactions:\n- equation: A + B => C\n  orders: $orders\n")
            @test_throws ArgumentError E.engine_check_order_metadata(path)
            diagnostic=joinpath(folder,"diagnostic")
            failure=try
                E.solve_ic_engine_qndf(path;diagnostic_sink=E.EngineFileSink(diagnostic))
            catch err
                err
            end
            @test failure isa E.EngineCalculationFailure
            @test failure.cause isa ArgumentError
            @test occursin("explicit orders",sprint(showerror,failure))
            @test haskey(failure.diagnostics,"calculation-exception")
            saved=E.NPZ.npzread(joinpath(diagnostic,"calculation-exception.npz"))
            @test occursin("explicit orders",String(saved["exception_utf8"]))
            @test !haskey(saved,"source_initial_state")
            @test !haskey(saved,"last_completed_endpoint")
            @test saved["recorded_intervals"]==[0]
        end
        for document in ("reactions: all\n","phases: []\n","reactions: [3]\n")
            write(path,document)
            @test_throws ArgumentError E.engine_check_order_metadata(path)
        end
    end
end

@testset "Engine real CLI loading and source scope" begin
    @test isdefined(EnginePublicCLI,:solve_ic_engine_qndf)
    @test isdefined(EnginePublicCLI,:write_ic_engine_csv)
    @test E.NPZ===E.Arrhenius.NPZ
    @test realpath(String(first(methods(E.engine_network)).file))==realpath(joinpath(reactors,"ic_engine_setup.jl"))
    @test realpath(String(first(methods(EnginePublicCLI.engine_network)).file))==realpath(joinpath(reactors,"ic_engine_setup.jl"))
    @test !isdefined(E,:engine_timing_main)
    @test !isdefined(E,:run_full_private_engine)
    @test !isdefined(E,:engine_remaining_probe)
    @test length(E.engine_switching_times(.16))==26
    @test last(E.engine_switching_times(.16))==.16
    @test_throws MethodError E.solve_ic_engine_qndf("unused.yaml";initial_state=zeros(105))
    @test_throws MethodError E.solve_ic_engine_qndf("unused.yaml";prefix="unused.npz")
    @test_throws MethodError E.solve_ic_engine_qndf("unused.yaml";times=[0.,.001])
    @test_throws MethodError E.solve_ic_engine_qndf("unused.yaml";reltol=1e-5)
    for name in ("entry.jl","adapter.jl","solve.jl","checks.jl")
        source=read(joinpath(reactors,"engine_qndf",name),String)
        @test !occursin(r"C:/|C:\\\\|/home/|/Users/|npzread|TOML.parsefile|cd\(",source)
    end
    cli=read(joinpath(reactors,"ic_engine.jl"),String)
    entry=read(joinpath(reactors,"engine_qndf","entry.jl"),String)
    @test !occursin("pkgdir(Arrhenius)",entry)
    @test occursin("import Arrhenius.NPZ",entry)
    @test !occursin(r"(?m)^(using|import) NPZ",entry)
    @test occursin("solve_ic_engine_qndf(ARGS[1];progress=true)",cli)
    @test occursin("NativeEngineQNDF.engine_cli_failure_guard",cli)
    @test occursin("write_ic_engine_csv",cli)
    @test occursin("TOML.print(io,summary;sorted=true)",cli)
end

@testset "Engine actual CLI failure guard without chemistry" begin
    mktempdir() do folder
        diagnostics=Dict{String,Any}("complete-checks"=>Dict("checks_pass"=>false,
            "solver_totals"=>Dict("accepted_steps"=>0)))
        primary=E.EngineCalculationFailure(ErrorException("intentional calculation failure"),diagnostics)
        directory=joinpath(folder,"saved")
        caught=try
            E.engine_cli_failure_guard(directory) do
                throw(primary)
            end
        catch err
            err
        end
        @test caught===primary
        @test !TOML.parsefile(joinpath(directory,"complete-checks.toml"))["checks_pass"]
        @test sprint(showerror,caught)=="intentional calculation failure"
        blocked=joinpath(folder,"a-file")
        write(blocked,"not a directory")
        caught=try
            E.engine_cli_failure_guard(joinpath(blocked,"child")) do
                throw(primary)
            end
        catch err
            err
        end
        @test caught===primary
        @test haskey(primary.diagnostics,"diagnostic-sink-exception")
        @test sprint(showerror,caught)=="intentional calculation failure"
        untouched=joinpath(folder,"unused-success-directory")
        value=E.engine_cli_failure_guard(()->42,untouched)
        @test value==42
        @test !ispath(untouched)
    end
end
