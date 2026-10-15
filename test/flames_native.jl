using Arrhenius, LinearAlgebra, Test

@testset "native premixed flames and transport" begin
    mechanism = joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    gas = CreateSolution(mechanism)
    X = mole_fractions(gas,"H2:1.1,O2:1,AR:5")
    transport = TransportWorkspace(gas)
    @test_throws ArgumentError mixture_transport!(transport,gas,one_atm,300.,zeros(gas.n_species))
    @test_throws ArgumentError mixture_transport!(transport,gas,one_atm,300.,X;basis=:invalid)
    @test_throws ArgumentError FreeFlame(gas;X,T=NaN)
    @test_throws ArgumentError FreeFlame(gas;X,width=-.1)
    @test_throws ArgumentError FreeFlame(gas;X,grid=[0.,.1,.1,.2,.3])
    @test_throws ArgumentError BurnerFlame(gas;X,mdot=0.)
    f = FreeFlame(gas;X,discretization=:finite_difference)
    @test_throws ArgumentError set_transport!(f,:multicomponent)
    @test_throws ArgumentError set_transport!(f,:mixture_averaged;soret=true)
    @test_throws ArgumentError solve!(f;initial_time_step=0.)
    @test_throws ArgumentError solve!(f;slope=0.)
    @test_throws ArgumentError solve!(f;prune=.1)
    solve!(f)
    @test f.converged
    # Cantera 4.0.0a2, commit 726522be4e2a13454d8415b7ef799d621f665cf3,
    # exact bundled 9-species mechanism, independent native initialization.
    @test flame_speed(f) ≈ .7055995685698858 rtol=.01
    @test maximum(temperature(f)) ≈ 1852.7139776420147 rtol=.01
    @test maximum(abs.(sum(mass_fractions(f);dims=1) .- 1)) < 1e-10
    @test minimum(mass_fractions(f)) > -1e-12
    @test maximum(f.state[end,:])-minimum(f.state[end,:]) < 1e-10
    @test norm(flame_residual!(similar(f.state),f),Inf) < 1e-8

    # solver_stats: per-flame deterministic counters. Cantera keeps Sim1D
    # solver_stats on the instance; each flame counts only its own solves and
    # the exported PREMIXED_STATS mirrors the accumulated family total.
    counts(s)=(s.jacobians,s.newton_calls,s.newton_iterations,s.time_steps,
        s.time_step_failures,s.steady_failures,s.regrids,copy(s.grids))
    stats=solver_stats(f)
    @test stats isa PremixedStats && solver_stats(f)===stats
    @test solver_stats(FreeFlame(gas;X))!==stats         # per-flame, not a shared global
    fcounts=counts(stats)
    @test fcounts[1]>0                                   # f's own solve reached its counters
    fresh=FreeFlame(gas;X,discretization=:finite_difference)
    @test counts(solver_stats(fresh))==(0,0,0,0,0,0,0,Int[])
    solve!(fresh)
    run1=counts(solver_stats(fresh))
    @test run1[1]>0 && run1[2]>0 && run1[3]>0            # jacobians, newton_calls, newton_iterations
    @test run1[3]>=run1[1]                                # every Jacobian serves a Newton iteration
    @test last(run1[8])==length(fresh.grid)               # final steady attempt kept the final grid
    @test counts(stats)==fcounts                          # solving fresh left f's counters untouched
    @test PREMIXED_STATS.jacobians>=run1[1]               # global mirrors the family total
    again=FreeFlame(gas;X,discretization=:finite_difference)
    solve!(again)
    @test counts(solver_stats(again))==run1               # equivalent cold solve is deterministic
    solve!(fresh;ratio=2.,slope=.02,curve=.05)            # a second solve! replaces the counters
    c2=counts(solver_stats(fresh))
    @test 0<c2[1]<run1[1]                                 # tighter criteria recount from zero
    @test c2[8]!=run1[8]                                  # the replaced grids differ from run1
    @test counts(stats)==fcounts
    # Burner flames count through the same shared helpers, per flame.
    burner=BurnerFlame(gas;X,mdot=.1,discretization=:finite_difference)
    bs=solver_stats(burner)
    @test bs!==stats && counts(bs)==(0,0,0,0,0,0,0,Int[])
    solve!(burner)
    @test burner.converged && bs.jacobians>0 && bs.newton_iterations>=bs.jacobians
    @test last(bs.grids)==length(burner.grid) && counts(stats)==fcounts
    # Fixed-temperature burner: same per-flame accounting on the imposed-profile path.
    fixed=BurnerFlame(gas;X,mdot=.1,discretization=:finite_difference)
    set_temperature_profile!(fixed,[0.,.25,1.],[300.,1600.,300.])
    fxs=solver_stats(fixed)
    solve!(fixed)
    @test fixed.converged && fxs.jacobians>0 && fxs.newton_iterations>=fxs.jacobians
    @test counts(stats)==fcounts

    # Cantera's fourth refine criterion: prune=0 must reproduce the default.
    explicit = solve!(FreeFlame(gas;X,discretization=:finite_difference);prune=0.)
    @test explicit.grid==f.grid && explicit.state==f.state
    pruned = solve!(FreeFlame(gas;X,discretization=:finite_difference);prune=.02)
    @test pruned.converged
    @test length(pruned.grid) < length(f.grid)
    @test flame_speed(pruned) ≈ .7055995685698858 rtol=.01

    @test maximum(abs.(density(f).*velocity(f) .- f.state[end,:])) < 1e-12
    @test maximum(heat_release_rate(f)) > 0
    @test maximum(abs.(sum(mole_fractions(f);dims=1) .- 1)) < 1e-12
    mktempdir() do directory
        snapshot = joinpath(directory,"flame.npz")
        table = joinpath(directory,"flame.csv")
        save_flame(snapshot,f)
        @test_throws ArgumentError save_flame(snapshot,f)
        restored = FreeFlame(gas;X)
        restore_flame!(restored,snapshot)
        @test restored.converged
        @test restored.state == f.state
        @test restored.grid == f.grid
        @test flame_speed(restored) == flame_speed(f)
        @test restored.discretization == :finite_difference
        # A restart must retain its equations, including when the destination
        # was constructed with another discretization. Old archives retain the
        # finite-difference semantics under which they were written.
        arrays = Arrhenius.npzread(snapshot)
        arrays["discretization_utf8"] = collect(codeunits("conservative"))
        Arrhenius.npzwrite(snapshot,arrays)
        restore_flame!(restored,snapshot)
        @test restored.discretization == :conservative
        delete!(arrays,"discretization_utf8")
        Arrhenius.npzwrite(snapshot,arrays)
        restore_flame!(restored,snapshot)
        @test restored.discretization == :finite_difference
        arrays["discretization_utf8"] = collect(codeunits("unknown"))
        Arrhenius.npzwrite(snapshot,arrays)
        @test_throws ArgumentError restore_flame!(restored,snapshot)
        save_flame(table,f;basis=:mole)
        lines = readlines(table)
        @test length(lines) == length(f.grid)+1
        @test startswith(lines[1],"z_m,velocity_m_s,temperature_K")
    end

    # Check the colored band Jacobian against a separate one-column difference,
    # holding transport fixed, as required by this modified Newton method.
    w = Arrhenius.FlameWorkspace(f)
    r = flame_residual!(similar(f.state),f,f.state,w)
    band = copy(Arrhenius._flame_jacobian(f,f.state,w,r))
    B,N = size(f.state)
    kl = 2B-1
    for (k,j) in ((1,f.anchor),(2,f.anchor+1),(B,1))
        u = copy(f.state)
        step = 1e-7*max(abs(u[k,j]),k==1 ? .1 : 1e-5)
        u[k,j] += step
        # Refresh the base state before the independent perturbation.
        flame_residual!(similar(r),f,f.state,w)
        rp = flame_residual!(similar(r),f,u,w;update_transport=false)
        column = (j-1)*B+k
        from_band = zeros(length(u))
        for row in max(1,column-kl):min(length(u),column+kl)
            from_band[row] = band[2kl+1+row-column,column]
        end
        @test from_band ≈ vec((rp-r)/step) rtol=1e-9 atol=1e-9
    end

    data = MultiTransportData(mechanism*".multicomponent.npz",gas;mechanism)
    wm = MultiTransportWorkspace(data)
    lambda = multicomponent_transport!(wm,data,gas,one_atm,900.,X)
    D = copy(wm.diffusion)
    DT = copy(wm.thermal_diffusion)
    @test lambda > 0
    @test abs(sum(DT)) < 1e-17
    @test multicomponent_transport!(wm,data,gas,2one_atm,900.,X) ≈ lambda rtol=1e-12
    @test wm.diffusion ≈ D/2 rtol=1e-12
    @test wm.thermal_diffusion ≈ DT rtol=1e-12
    @test multicomponent_thermal_conductivity!(wm,data,gas,one_atm,900.,X) == lambda
    @test wm.diffusion ≈ D/2 rtol=1e-12 # conductivity-only call leaves diffusion intact
    @test_throws ArgumentError multicomponent_thermal_conductivity!(wm,data,one_atm,900.,2X,wm.cp_R)
    flux = zeros(gas.n_species)
    gradX = collect(1.:gas.n_species); gradX .-= sum(gradX)/length(gradX)
    multicomponent_fluxes!(flux,wm,data,2one_atm,900.,X,gradX,10.)
    @test abs(sum(flux)) < 1e-13
    set_transport!(f,:multicomponent;data)
    @test !f.converged
    solve!(f)
    @test flame_speed(f) ≈ .7204474366630716 rtol=.01
    @test norm(flame_residual!(similar(f.state),f),Inf) < 1e-8

    burner = BurnerFlame(gas;X="H2:1.5,O2:1,AR:7",T=373.,P=.05one_atm,mdot=.06,width=.5,discretization=:finite_difference)
    solve!(burner;slope=.05,curve=.1)
    @test burner.converged
    @test maximum(temperature(burner)) ≈ 1857.234327060171 rtol=.01
    @test norm(flame_residual!(similar(burner.state),burner),Inf) < 1e-8
    @test minimum(mass_fractions(burner)) > -1e-12
    # A supplied mesh is part of the user's boundary-value problem. Its domain
    # and coordinates must survive construction, including a shifted origin.
    supplied_grid = .1 .+ .5 .* [0,.1,.2,.3,.5,.7,1]
    supplied = BurnerFlame(gas;X="H2:1.5,O2:1,AR:7",T=373.,P=.05one_atm,mdot=.06,grid=supplied_grid,discretization=:finite_difference)
    @test supplied.grid == supplied_grid
    @test supplied.state[1,1] == .373
    @test supplied.state[end,:] == fill(.06,length(supplied_grid))
    @test_throws ArgumentError set_temperature_profile!(burner,[0.,1.],[300.,1750.])
    @test_throws ArgumentError set_temperature_profile!(burner,[0.,.5],[373.,1750.])
    positions = [0.,.005,.01,.02,.05,.1,1.]
    temperatures = [373.,650.,1000.,1350.,1650.,1750.,1750.]
    set_temperature_profile!(burner,positions,temperatures)
    solve!(burner;slope=.05,curve=.1)
    @test burner.converged
    @test maximum(temperature(burner)) ≈ 1750.
    @test maximum(abs.(burner.state[end,:] .- .06)) < 1e-10
    @test norm(flame_residual!(similar(burner.state),burner),Inf) < 1e-8
    @test temperature(burner) ≈ [Arrhenius._interpolate_profile(
        .5positions,temperatures,z) for z in burner.grid]
end
