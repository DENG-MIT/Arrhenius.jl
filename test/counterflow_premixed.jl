using Test, Arrhenius, LinearAlgebra
if !isdefined(Arrhenius,:CounterflowPremixedFlame)
    Base.include(Arrhenius,joinpath(@__DIR__,"..","src","PremixedCounterflowFlames.jl"))
end

@testset "premixed axisymmetric boundaries" begin
    gas=CreateSolution(joinpath(dirname(pathof(Arrhenius)),"..","mechanism","h2o2.yaml"))
    composition="H2:1.6,O2:1,AR:7"
    kwargs=(;reactants=composition,mdot_reactants=.12,mdot_products=.06,
        T_reactants=373.,P=.05*one_atm,width=.2)
    @test_throws ArgumentError Arrhenius.CounterflowPremixedFlame(gas;kwargs...,products="H2O:1")
    @test_throws ArgumentError Arrhenius.CounterflowTwinPremixedFlame(gas;reactants=composition,mdot=-.1)
    @test_throws ArgumentError Arrhenius.ImpingingJet(gas;reactants=composition,mdot=.06,T_surface=0.)
    f=Arrhenius.CounterflowPremixedFlame(gas;kwargs...)
    @test !Arrhenius._flame_auto_species_jacobian!(nothing,f,nothing,nothing,nothing)
    solve!(f;slope=.1,curve=.2)
    @test f.converged
    @test 2000<maximum(temperature(f))<2100
    @test velocity(f)[1]>0 && velocity(f)[end]<0
    @test f.state[end,[1,end]] ≈ [.12,-.06] atol=1e-9
    @test maximum(abs.(sum(mass_fractions(f);dims=1).-1))<1e-10
    r=similar(f.state); Arrhenius.counterflow_residual!(r,f)
    @test norm(r,Inf)<1e-8
    # Per-kind nonzero solve counters: solver_stats(f) counts only this flame's
    # last solve! (Cantera keeps Sim1D.solver_stats on the instance).
    ps=Arrhenius.solver_stats(f)
    @test ps isa CounterflowStats && Arrhenius.solver_stats(f)===ps
    @test ps.jacobians>0 && ps.newton_iterations>=ps.jacobians
    @test !isempty(ps.grids) && last(ps.grids)==length(f.grid)
    fc=pseudo_arclength!(f;parameter=:fuel_mass_flux,steps=1,ds=.005,tolerance=1e-7)
    @test length(fc)==2 && all(p->p.residual_inf<1e-7,fc)
    premixed_cp_dir=mktempdir(); premixed_cp=joinpath(premixed_cp_dir,"premixed.npz")
    save_continuation_checkpoint(premixed_cp,f,fc[1];parameter=:fuel_mass_flux)
    f_resume=Arrhenius.CounterflowPremixedFlame(gas;kwargs...)
    fc_resume=pseudo_arclength!(f_resume;restart=premixed_cp,steps=1,tolerance=1e-7)
    rm(premixed_cp_dir;recursive=true,force=true)
    @test fc_resume[1].state==fc[1].state
    @test fc_resume[2].parameter≈fc[2].parameter rtol=1e-12
    @test fc_resume[2].state≈fc[2].state rtol=1e-9 atol=1e-10

    # Inert wall and symmetry are checked directly with signed trial traces.
    twin=Arrhenius.CounterflowTwinPremixedFlame(gas;reactants=composition,mdot=.12,T=373.,P=.05*one_atm,width=.2)
    wall=Arrhenius.ImpingingJet(gas;reactants=composition,mdot=.06,T_inlet=373.,T_surface=500.,P=.05*one_atm,width=.2)
    for flow in (twin,wall)
        rr=similar(flow.state)
        Arrhenius.counterflow_residual!(rr,flow)
        @test rr[end,end]==0
        @test maximum(abs.(rr[2:gas.n_species+1,end]))<1e-12
        @test_throws ArgumentError Arrhenius.set_mass_flux!(flow;products=.01)
        Arrhenius.set_mass_flux!(flow;reactants=.07)
        @test flow.fuel_mass_flux==.07 && !flow.converged
        @test Arrhenius._premixed_reinitialize(flow,flow.grid).fuel_mass_flux==.07
    end
    rt=similar(twin.state); Arrhenius.counterflow_residual!(rt,twin)
    @test rt[1,end]==0
    @test rt[end-2,end] ≈ twin.state[end-2,end]-twin.state[end-2,end-1]
    # Symm1D applies its species equations independently of Soret transport.
    mechanism=joinpath(dirname(pathof(Arrhenius)),"..","mechanism","h2o2.yaml")
    transport=MultiTransportData(mechanism*".multicomponent.npz",gas;mechanism)
    twin.transport_model=:multicomponent
    twin.multicomponent_data=transport
    twin.soret_enabled=true
    tw=Arrhenius.CounterflowWorkspace(twin); rs=similar(twin.state)
    Arrhenius.counterflow_residual!(rs,twin,twin.state,tw)
    iex=tw.excess_right
    for k in 1:gas.n_species
        k==iex || @test rs[k+1,end] ≈ twin.state[k+1,end]-twin.state[k+1,end-1]
    end
    twin.transport_model=:mixture_averaged
    twin.multicomponent_data=nothing
    twin.soret_enabled=false
    hydrogen=findfirst(==("H2"),gas.species_names)
    wall.state[hydrogen+1,end]+=1e-4
    wall_workspace=Arrhenius.CounterflowWorkspace(wall)
    rw=similar(wall.state); Arrhenius.counterflow_residual!(rw,wall,wall.state,wall_workspace)
    @test rw[hydrogen+1,end] ≈ -6.969555716579648e-7 rtol=1e-10
    @test rw[hydrogen+1,end] ≈ wall_workspace.properties.flux[hydrogen,end]+
        wall.state[end,end]*wall.state[hydrogen+1,end] atol=1e-14
    @test rw[1,end] ≈ 0 atol=1e-14
    @test rw[end-2,end]==0
    @test_throws ArgumentError solve!(wall;prune=.5)
    @test_throws ArgumentError solve!(wall;grid_min=0.)

    # Twin-flame kind: counters of its own solve! (the mdot was set to .07 above).
    ts=Arrhenius.solver_stats(twin)
    @test ts.jacobians==0 && isempty(ts.grids)          # fresh flame starts at zero
    solve!(twin;slope=.1,curve=.2)
    @test twin.converged
    @test ts.jacobians>0 && ts.newton_iterations>=ts.jacobians
    @test !isempty(ts.grids) && last(ts.grids)==length(twin.grid)
    tc=pseudo_arclength!(twin;parameter=:fuel_mass_flux,steps=1,ds=.005,tolerance=1e-7)
    @test length(tc)==2 && all(p->p.residual_inf<1e-7,tc)
    twin_cp_dir=mktempdir(); twin_cp=joinpath(twin_cp_dir,"twin.npz")
    save_continuation_checkpoint(twin_cp,twin,tc[1];parameter=:fuel_mass_flux)
    twin_resume=Arrhenius.CounterflowTwinPremixedFlame(gas;reactants=composition,mdot=.12,T=373.,P=.05*one_atm,width=.2)
    tc_resume=pseudo_arclength!(twin_resume;restart=twin_cp,steps=1,tolerance=1e-7)
    rm(twin_cp_dir;recursive=true,force=true)
    @test tc_resume[1].state==tc[1].state
    @test tc_resume[2].parameter≈tc[2].parameter rtol=1e-12
    @test tc_resume[2].state≈tc[2].state rtol=1e-9 atol=1e-10

    cold=Arrhenius.ImpingingJet(gas;reactants="AR:1",mdot=.06,T_inlet=373.,
        T_surface=500.,P=.05*one_atm,width=.2,initial_products=:inlet)
    solve!(cold)
    @test cold.converged
    @test Arrhenius.extinct(cold)
    @test temperature(cold)[[1,end]] ≈ [373.,500.] atol=1e-6
    @test minimum(temperature(cold)) >= 373. - 1e-6
    @test maximum(temperature(cold)) <= 500. + 1e-6
    @test abs(velocity(cold)[end])<1e-10
    @test maximum(abs.(mass_fractions(cold)[findfirst(==("AR"),gas.species_names),:].-1))<1e-10
    # Impinging-jet kind: nonzero counters, and a warm solve! REPLACES them.
    js=Arrhenius.solver_stats(cold)
    @test js.jacobians>0 && js.newton_iterations>=js.jacobians
    @test last(js.grids)==length(cold.grid)
    cold_jac=js.jacobians
    solve!(cold)
    @test cold.converged && js.jacobians<cold_jac
end

@testset "premixed steady retry counts successful steps" begin
    # Wiring-layer regression for the solve! cadence: Cantera retries the
    # steady Newton solve only after m_steps={10} successful pseudo-time
    # steps (SteadyStateSystem::timeStep increments n on success only), so a
    # failed step advances the attempt count without counting toward the
    # retry. COUNTERFLOW_STATS.steady_retries records (successful steps,
    # attempts) at every retry decision, counted since the previous retry of
    # the same transient loop. An attempt-counted wiring instead retries at
    # attempt multiples and records fewer than ten successes whenever
    # failures intervene (measured 9, 7 and 6 on this solve in the prior
    # attempt-cadence intervention).
    gas=CreateSolution(joinpath(dirname(pathof(Arrhenius)),"..","mechanism","h2o2.yaml"))
    f=Arrhenius.CounterflowPremixedFlame(gas;reactants="H2:1.6,O2:1,AR:7",
        mdot_reactants=.12,mdot_products=.06,T_reactants=373.,P=.05*one_atm,width=.2)
    Arrhenius.reset!(Arrhenius.COUNTERFLOW_STATS)
    solve!(f;slope=.1,curve=.2)
    @test f.converged
    retries=Arrhenius.COUNTERFLOW_STATS.steady_retries
    # the transient path with its steady retries must actually be exercised
    @test !isempty(retries)
    @test all(first.(retries).==10)
    # failures intervene between retries on this solve (windows wider than
    # their ten accepted steps), so the count above discriminates the two
    # cadences instead of passing vacuously. The attempt counter is retained
    # only for this telemetry.
    @test Arrhenius.COUNTERFLOW_STATS.time_step_failures>0
    @test any(last.(retries).>10)
end

@testset "twin excess index is shared by residual, growth, and Jacobian" begin
    gas=CreateSolution(joinpath(dirname(pathof(Arrhenius)),"..","mechanism","h2o2.yaml"))
    twin=CounterflowTwinPremixedFlame(gas;reactants="H2:1.6,O2:1,AR:7",
        mdot=.12,T=373.,P=.05*one_atm,width=.2)
    u=copy(twin.state); n=gas.n_species; N=length(twin.grid)
    w=Arrhenius.CounterflowWorkspace(twin); r=similar(u); growth=similar(u)
    # Make two major species nearly tied at the symmetry plane and record the
    # index from a full residual evaluation.
    u[2:n+1,end].=0
    selected=1; competitor=2
    u[selected+1,end]=.500000000001
    u[competitor+1,end]=.499999999999
    Arrhenius.counterflow_residual!(r,twin,u,w)
    frozen=w.excess_right
    @test frozen==selected

    # The trial state crosses the argmax, but partial evaluations, including
    # the growth residual, must retain the index selected above.
    x=copy(u)
    x[selected+1,end]-=4e-12
    x[competitor+1,end]+=4e-12
    @test argmax(@view x[2:n+1,end])==competitor
    Arrhenius.counterflow_residual!(r,twin,x,w;update_excess=false)
    @test w.excess_right==frozen
    Arrhenius._counterflow_growth_residual!(growth,twin,x,w)
    @test w.excess_right==frozen
    @test r[frozen+1,end]==sum(@view(x[2:n+1,end]))-1
    @test growth[frozen+1,end]==1-sum(@view(x[2:n+1,end]))
    for k in 1:n
        k==frozen && continue
        @test r[k+1,end]==x[k+1,end]-x[k+1,end-1]
        @test growth[k+1,end]==x[k+1,end]-x[k+1,end-1]
    end

    # Colored finite differences may cross the near-tie too; the Jacobian
    # closure row must remain the same frozen normalization row.
    jac=Arrhenius._counterflow_jacobian!(twin,x,w,r;analytic=false)
    @test w.excess_right==frozen
    B=size(x,1); stride=4B-2
    for k in 1:n
        col=(N-1)*B+k+1; row=(N-1)*B+frozen+1
        packed=stride+1+row-col
        h=1e-7*max(abs(x[k+1,end]),k+1<=n+1 ? 1e-5 : .01)
        trial=copy(x); trial[k+1,end]+=h
        fd=similar(x)
        Arrhenius.counterflow_residual!(fd,twin,trial,w;update_excess=false)
        expected=(fd[frozen+1,end]-r[frozen+1,end])/h
        @test jac[packed,col]≈expected atol=2e-8
    end
end
