using Test
function engine_qndf_solver_contract()
    norm=FullEngineRMS(105)
    @test Base.broadcastable(norm) isa Base.RefValue{FullEngineRMS}
    @test ((f,x)->f(x,0.)).(norm,[-2.,0.,3.])==[2.,0.,3.]
    @test norm(ones(8),0.)≈sqrt(8/105) rtol=8eps()
    @test norm(abs,ones(8),0.)≈sqrt(8/105) rtol=8eps()
    base=only(m for (id,m) in Base.loaded_modules if id.name=="DiffEqBase")
    # The installed in-place vector-atol path that failed before the first
    # accepted engine step must broadcast the scalar callable successfully.
    u0=collect(1.:8.);u1=reverse(u0);residual_error=fill(2e-14,8);atol=fill(1e-16,8);out=zeros(8)
    base.calculate_residuals!(out,residual_error,u0,u1,atol,1e-12,norm,0.)
    @test out≈residual_error./(atol.+max.(abs.(u0),abs.(u1)).*1e-12) rtol=8eps()
    # Scalar-atol and difference-residual variants also call the scalar norm.
    base.calculate_residuals!(out,residual_error,u0,u1,1e-16,1e-12,norm,0.)
    @test out≈residual_error./(atol.+max.(abs.(u0),abs.(u1)).*1e-12) rtol=8eps()
    base.calculate_residuals!(out,u0,u1,atol,1e-12,norm,0.)
    @test out≈(u1.-u0)./(atol.+max.(abs.(u0),abs.(u1)).*1e-12) rtol=8eps()
    f=(du,u,p,t)->(du.=-u;nothing)
    jac=(J,u,p,t)->begin
        fill!(J,0.)
        for i in axes(J,1);J[i,i]=-1.;end
        nothing
    end
    initial=vcat([1.,2.],zeros(6))
    ode=ODEProblem(ODEFunction(f;jac),initial,(0.,.02))
    mktempdir() do folder
        diagnostics=EngineDiagnosticLog(EngineFileSink(folder))
        reconstruct=u->vcat(copy(u),zeros(97))
        trace=Dict{String,Any}("initial_state"=>reconstruct(initial))
        path=joinpath(folder,"normal.npz")
        options=(;internalnorm=norm,reltol=1e-12,abstol=fill(1e-16,8),
            dt=1e-4,dtmax=.01,dense=false,save_everystep=true)
        solution=engine_probe_solve(ode,QNDF(),trace,diagnostics,"normal",reconstruct;options...)
        @test SciMLBase.successful_retcode(solution)
        @test solution.u[end]≈initial.*exp(-.02) rtol=1e-10 atol=1e-16
        @test all(iszero,solution.u[end][3:end])
        @test !isfile(path)
        # Exercise the same exception serializer after initialization/advance.
        failed_trace=copy(trace);failed_path=joinpath(folder,"accepted_failure.npz")
        callback=DiscreteCallback((u,t,integrator)->true,
            integrator->error("intentional accepted-step contract failure");save_positions=(false,false))
        @test_throws ErrorException engine_probe_solve(ode,QNDF(),failed_trace,diagnostics,"accepted_failure",reconstruct;
            options...,callback,time_offset=.001)
        saved=NPZ.npzread(failed_path)
        @test saved["initial_state"]==reconstruct(initial)
        @test length(saved["current_state"])==105
        @test size(saved["saved_states"],1)==105
        @test size(saved["saved_states"],2)==length(saved["saved_times"])
        @test only(saved["current_time"])>.001
        @test occursin("intentional accepted-step contract failure",String(saved["exception_utf8"]))
        # An initialization error must retain the exact initial state even
        # when no integrator exists to supply current or accepted states.
        broken=ODEProblem((du,u,p,t)->error("intentional initialization contract failure"),initial,(0.,.02))
        init_trace=copy(trace);init_path=joinpath(folder,"initialization_failure.npz")
        @test_throws ErrorException engine_probe_solve(broken,QNDF(),init_trace,diagnostics,"initialization_failure",reconstruct;options...)
        init_saved=NPZ.npzread(init_path)
        @test init_saved["initial_state"]==reconstruct(initial)
        @test !haskey(init_saved,"current_state")
        @test occursin("intentional initialization contract failure",String(init_saved["exception_utf8"]))
    end
    true
end
