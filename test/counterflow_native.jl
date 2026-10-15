using Test, Arrhenius, LinearAlgebra, NPZ
if !isdefined(Arrhenius,:CounterflowDiffusionFlame)
    Base.include(Arrhenius,joinpath(@__DIR__,"..","src","CounterflowFlames.jl"))
end
@testset "native axisymmetric counterflow" begin
    mechanism=joinpath(dirname(pathof(Arrhenius)),"..","mechanism","h2o2.yaml")
    if isfile(mechanism*".npz")
        gas=CreateSolution(mechanism)
        kwargs=(fuel="H2:1,AR:1",oxidizer="O2:.2,AR:.8",mdot_fuel=.24,mdot_oxidizer=.72)
        @test_throws ArgumentError Arrhenius.CounterflowDiffusionFlame(gas;kwargs...,P=-1)
        @test_throws ArgumentError Arrhenius.CounterflowDiffusionFlame(gas;kwargs...,grid=[0.,.01,.005,.015,.02])
        @test_throws ArgumentError Arrhenius.CounterflowDiffusionFlame(gas;kwargs...,boundary_emissivities=(0.,1.1))
        f=Arrhenius.CounterflowDiffusionFlame(gas;kwargs...)
        solve!(f;slope=.2,curve=.3)
        @test f.converged
        @test 2300<maximum(temperature(f))<2600
        @test temperature(f)[[1,end]] ≈ [300.,300.] atol=1e-6
        @test f.state[end,1] ≈ .24 atol=1e-8
        @test f.state[end,end] ≈ -.72 atol=1e-8
        @test velocity(f)[1]>0 && velocity(f)[end]<0
        @test minimum(mass_fractions(f))>-1e-7
        @test maximum(abs.(sum(mass_fractions(f);dims=1).-1))<1e-10
        @test maximum(abs.(Arrhenius.pressure_curvature(f).-Arrhenius.pressure_curvature(f)[1]))<1e-7
        r=similar(f.state)
        Arrhenius.counterflow_residual!(r,f)
        @test norm(r,Inf)<1e-8
        # solver_stats: per-flame counters; the exported global mirrors the family.
        cstats=solver_stats(f)
        @test cstats isa CounterflowStats && solver_stats(f)===cstats
        @test cstats.jacobians>0 && !isempty(cstats.grids)
        @test cstats.newton_iterations>=cstats.jacobians
        @test last(cstats.grids)==length(f.grid)
        @test COUNTERFLOW_STATS.jacobians>=cstats.jacobians
        # residual_evals counts every counterflow_residual! call at entry,
        # manual ones included: three extra calls add exactly three.
        re=cstats.residual_evals
        r3=similar(f.state)
        for _ in 1:3
            counterflow_residual!(r3,f)
        end
        @test cstats.residual_evals==re+3
        if haskey(ENV,"COUNTERFLOW_REFERENCE_DIR")
            ref=npzread(joinpath(ENV["COUNTERFLOW_REFERENCE_DIR"],"counterflow-h2-reference.npz"))
            interp(values)=[Arrhenius._interpolate_profile(ref["grid"],vec(values),z) for z in f.grid]
            @test maximum(abs.(temperature(f).-interp(ref["T"])))<1
            @test maximum(abs.(velocity(f).-interp(ref["velocity"])))<1e-3
            @test Arrhenius.pressure_curvature(f)[1] ≈ ref["Lambda"][1] rtol=1e-4
            @test maximum(abs.(mass_fractions(f).-reduce(vcat,permutedims(interp(ref["Y"][k,:])) for k in 1:gas.n_species)))<1e-4
        end
        Tbefore=maximum(temperature(f)); grid=copy(f.grid)
        @test all(iszero,Arrhenius.radiative_heat_loss(f))
        Arrhenius.set_radiation!(f,true;boundary_emissivities=(.3,.7))
        solve!(f;refine_grid=false)
        @test f.converged && f.grid==grid
        @test maximum(temperature(f))<Tbefore
        source=Arrhenius.radiation_source(f)
        @test maximum(source.heat_loss)>0
        @test source.heat_loss[end]==0
        @test maximum(source.planck_absorption)>0

        @testset "two-point control equations and compatibility" begin
            cgas=CreateSolution(mechanism)
            cf=CounterflowDiffusionFlame(cgas;fuel="H2:1",oxidizer="O2:1",mdot_fuel=.5,
                mdot_oxidizer=3.,T_fuel=300.,T_oxidizer=500.,P=1e5,
                grid=collect(range(0.,.018;length=7)))
            cf.state[1,:]=[.3,1.,1.8,2.6,2.,1.2,.5]
            # The Cantera-parity initial guess stores rho*u of the blended
            # boundary state; pin exact inlet fluxes for the control checks.
            cf.state[end,1]=.5; cf.state[end,end]=-3.
            state0=copy(cf.state); grid0=copy(cf.grid); n=cgas.n_species
            @test_throws ArgumentError set_two_point_control!(cf;temperature=Inf)
            @test_throws ArgumentError set_two_point_control!(cf;temperature=150.)
            @test_throws ArgumentError set_two_point_control!(cf;temperature=1500.,decrement=-1.)
            @test_throws ArgumentError set_two_point_control!(cf;temperature=1500.,decrement=1400.)
            @test cf.state==state0 && cf.grid==grid0 && isnothing(cf.control_points)
            cf.fixed_temperature=temperature(cf)
            @test_throws ArgumentError set_two_point_control!(cf;temperature=1500.)
            @test cf.state==state0 && isnothing(cf.control_points)
            empty!(cf.fixed_temperature)
            set_two_point_control!(cf;temperature=1500.,decrement=20.)
            # Direct internal-helper calls accumulate without resetting, and
            # solve! resets at entry: the two-point/fixed-temperature
            # compatibility check runs AFTER that reset, so the throwing
            # solve! zeroes the counters the helper call had filled.
            st=copy(cf.state)
            Arrhenius._counterflow_newton!(cf,Arrhenius.CounterflowWorkspace(cf))
            cs=solver_stats(cf)
            @test cs.newton_calls==1
            rr2=similar(cf.state); re2=cs.residual_evals
            for _ in 1:3
                counterflow_residual!(rr2,cf)
            end
            @test cs.residual_evals==re2+3
            copyto!(cf.state,st)
            cf.fixed_temperature=temperature(cf)
            @test_throws ArgumentError solve!(cf;refine_grid=false)
            @test all(iszero,getfield(cs,n) for n in fieldnames(CounterflowStats)
                if fieldtype(CounterflowStats,n)==Int)
            @test isempty(cs.grids) && isempty(cs.calls)
            empty!(cf.fixed_temperature)
            zL,tL,zR,tR=cf.control_points; iL=findfirst(==(zL),cf.grid); iR=findfirst(==(zR),cf.grid)
            @test iL<iR
            @test tL==1000*state0[1,iL]-20
            @test tR==1000*state0[1,iR]-20
            rr=similar(cf.state); counterflow_residual!(rr,cf)
            @test rr[end-1,iL]≈.02 atol=1e-14
            @test rr[n+2,iR]≈.02 atol=1e-14
            @test rr[end,end]≈0 atol=1e-14
            @test all(isfinite,rr)
            controlled_state=copy(cf.state)
            qdot=heat_release_rate(cf)
            @test length(qdot)==length(cf.grid) && all(isfinite,qdot)
            @test cf.state==controlled_state
            disable_two_point_control!(cf)
            @test cf.state==state0 && cf.grid==grid0 && isnothing(cf.control_points)
            @test cf.fuel_mass_flux==.5 && cf.oxidizer_mass_flux==3.

            set_two_point_control!(cf;temperature=1500.,decrement=20.)
            trace=findfirst(==("HO2"),cgas.species_names)+1
            for j in 2:6
                δ=(iseven(j) ? 0. : -1e-12)-cf.state[trace,j]
                cf.state[trace,j]+=δ; cf.state[cf.dependent_species+1,j]-=δ
            end
            w=CounterflowWorkspace(cf); base=similar(cf.state)
            counterflow_residual!(base,cf,cf.state,w)
            band=copy(Arrhenius._counterflow_jacobian!(cf,cf.state,w,base;analytic=false))
            hybrid=copy(Arrhenius._counterflow_jacobian!(cf,cf.state,w,base))
            B,N=size(cf.state); bw=2B-1; J=zeros(length(cf.state),length(cf.state))
            for col in 1:length(cf.state), row in max(1,col-bw):min(length(cf.state),col+bw)
                J[row,col]=band[2bw+1+row-col,col]
            end
            ref=zeros(size(J)); wp=CounterflowWorkspace(cf); rp=similar(cf.state); trial=copy(cf.state)
            counterflow_residual!(rp,cf,cf.state,wp)
            for j in 1:N, k in 1:B
                h=1e-7*max(abs(cf.state[k,j]),k==1 ? .1 : k<=n+1 ? 1e-5 : .01)
                trial .= cf.state; trial[k,j]+=h
                counterflow_residual!(rp,cf,trial,wp;update_transport=false)
                ref[:,(j-1)*B+k].=vec((rp.-base)./h)
            end
            @test norm(J-ref,Inf)/max(1.,norm(ref,Inf))<1e-7
            @test maximum(abs.(J.-ref)./(1e-5.+abs.(ref)))<1e-4
            # Analytic local chemistry replaces only the chemistry finite differences;
            # their truncation error bounds the hybrid difference.
            H=zeros(size(J))
            for col in 1:length(cf.state), row in max(1,col-bw):min(length(cf.state),col+bw)
                H[row,col]=hybrid[2bw+1+row-col,col]
            end
            @test norm(H-ref,Inf)/max(1.,norm(ref,Inf))<1e-5
            significant=abs.(ref).>1e-3*norm(ref,Inf)
            @test maximum(abs.(H.-ref)[significant]./abs.(ref)[significant])<1e-3
            @test J[(N-1)*B+B,(N-1)*B+n+2]≈1 atol=1e-6
            @test J[B-1,B]==0
            @testset "fresh correction convergence guards" begin
                scales=fill(NaN,length(cf.state))
                Arrhenius._counterflow_residual_scales!(scales,band,cf.state,bw)
                @test scales ≈ max.(1.,abs.(J)*abs.(vec(cf.state))) rtol=1e-14
                Arrhenius._counterflow_residual_scales!(scales,band,zero(cf.state),bw)
                @test all(==(1.),scales)

                guards=zeros(size(cf.state)); tol=1e-8
                @test Arrhenius._counterflow_constraints_converged(cf,guards,tol)
                for (row,col) in ((1,1),(1,N),(B-1,3),(B,3),(n+2,3),
                                  (cf.dependent_species+1,3))
                    guards[row,col]=2tol
                    @test !Arrhenius._counterflow_constraints_converged(cf,guards,tol)
                    guards[row,col]=0.
                end
                guards[cf.dependent_species+1,3]=2e-9
                @test !Arrhenius._counterflow_constraints_converged(cf,guards,tol)
                guards .= 0.; guards[B-2,3]=1.
                @test Arrhenius._counterflow_constraints_converged(cf,guards,tol)
                guards .= 0.
                for value in (NaN,Inf)
                    guards[n+2,3]=value
                    @test !Arrhenius._counterflow_constraints_converged(cf,guards,tol)
                end
            end
            cp=cf.control_points
            @test Arrhenius._refine_flame!(cf;ratio=4.,slope=.1,curve=.2,max_points=100)
            @test cf.control_points==cp
            @test cp[1] in cf.grid && cp[3] in cf.grid
            @testset "optional pruning" begin
                flat=CounterflowDiffusionFlame(cgas;fuel="H2:1",oxidizer="O2:1",mdot_fuel=.5,mdot_oxidizer=3.,
                    grid=collect(range(0.,.018;length=9)))
                flat.state .= flat.state[:,1]
                flat_grid=copy(flat.grid); flat_state=copy(flat.state)
                @test !Arrhenius._refine_flame!(flat;prune=0.)
                @test flat.grid==flat_grid && flat.state==flat_state
                @test !Arrhenius._refine_flame!(flat;prune=.05)
                @test flat.grid==flat_grid && flat.state==flat_state
                flat.state[1,:]=[.3,.3,.3,.3,1.,1.,.3,.3,.3]
                flat_state=copy(flat.state)
                @test Arrhenius._refine_flame!(flat;slope=1.,curve=1.,prune=.5)
                @test length(flat.grid)<length(flat_grid)
                @test flat.grid[[1,end]]==flat_grid[[1,end]] && all(diff(flat.grid).>0)
                @test all(z->z in flat_grid,flat.grid)
                @test flat.state==flat_state[:,[findfirst(==(z),flat_grid) for z in flat.grid]]
                removed=findall(z->!(z in flat.grid),flat_grid)
                @test all(diff(removed).>1)
                Arrhenius._refine_flame!(cf;ratio=4.,slope=.1,curve=.2,prune=.05,max_points=150)
                @test cf.control_points==cp && cp[1] in cf.grid && cp[3] in cf.grid
            end
            premixed=CounterflowPremixedFlame(cgas;reactants="H2:2,O2:1,AR:7",
                mdot_reactants=.12,mdot_products=.06)
            @test isnothing(premixed.flow.control_points)
            @test size(premixed.flow.state,1)==n+4
        end
    else
        @test_skip false
    end
end

@testset "nakamura counterflow damping regression" begin
    # Solver-decision regression for the nakamura_2017 slowdown: the counterflow
    # damped Newton must follow Cantera MultiNewton::dampStep (<=7 trials, sqrt(2)
    # decay, then Jacobian refresh). With the original 24-halving search the
    # steady Newton stalled on stiff trace species (quadratic reverse rates such
    # as N2H2 + M <=> 2 NH + M) and this case needed 550 Jacobians and 560
    # pseudo-time steps; Cantera needs 77 and 70 (nakcf task counts).
    mechanism=joinpath(@__DIR__,"data","nakamura_2017.yaml")  # TODO: pinned CollectionOfMechanisms artifact
    if isfile(mechanism*".npz")
        gas=CreateSolution(mechanism)
        f=CounterflowDiffusionFlame(gas;fuel="NH3:0.7,H2:0.3",oxidizer="O2:0.21,N2:0.79",
            mdot_fuel=.2,mdot_oxidizer=.6,T_fuel=300.,T_oxidizer=300.,P=one_atm,width=.02)
        Arrhenius.reset!(Arrhenius.COUNTERFLOW_STATS)
        solve!(f;ratio=3.,slope=.03,curve=.06,prune=.02)
        @test f.converged
        @test maximum(temperature(f)) ≈ 2121.46 rtol=1e-3
        stats=Arrhenius.COUNTERFLOW_STATS
        @test stats.jacobians < 250
        @test stats.time_steps < 250
        @test stats.recoveries == 0
        r=similar(f.state); Arrhenius.counterflow_residual!(r,f)
        # Steady convergence may end on Cantera weighted-correction criterion rather
        # than the 1e-8 residual floor; 1e-5 separates converged from stalled.
        @test norm(r,Inf)<1e-5
        @test minimum(mass_fractions(f))>-1e-7
        @test maximum(abs.(sum(mass_fractions(f);dims=1).-1))<1e-10
    else
        @test_skip false
    end
end

@testset "solver_stats dispatch for every counterflow kind" begin
    # Constructing each public kind must yield its own stats object (Cantera
    # keeps solver_stats per instance); dispatch needs no solve.
    gas=CreateSolution(joinpath(dirname(pathof(Arrhenius)),"..","mechanism","h2o2.yaml"))
    diffusion=CounterflowDiffusionFlame(gas;fuel="H2:1,AR:1",oxidizer="O2:.2,AR:.8",
        mdot_fuel=.24,mdot_oxidizer=.72,P=one_atm,width=.02)
    composition="H2:1.6,O2:1,AR:7"; P=.05*one_atm
    premixed=CounterflowPremixedFlame(gas;reactants=composition,mdot_reactants=.12,
        mdot_products=.06,T_reactants=373.,P,width=.2)
    twin=CounterflowTwinPremixedFlame(gas;reactants=composition,mdot=.12,T=373.,P,width=.2)
    jet=ImpingingJet(gas;reactants=composition,mdot=.06,T_inlet=373.,T_surface=500.,P,width=.2)
    flames=(diffusion,premixed,twin,jet)
    @test all(f -> solver_stats(f) isa CounterflowStats, flames)
    stats=solver_stats.(flames)
    @test all(stats[i]!==stats[j] for i in 1:4 for j in i+1:4)
    @test all(s -> s.jacobians==0 && isempty(s.grids), stats)  # fresh flames start at zero
    # One small h2o2 solve per kind: per-flame counters come back nonzero.
    for f in flames
        s=solver_stats(f)
        solve!(f;ratio=8.,slope=.3,curve=.5)
        @test f.converged
        @test s.jacobians>0 && s.newton_iterations>0
        @test s.newton_iterations>=s.jacobians
        @test last(s.grids)==length(f.grid)
    end
end
