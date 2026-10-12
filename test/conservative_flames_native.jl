using Arrhenius,LinearAlgebra,Test

@testset "conservative premixed species and enthalpy balances" begin
    mechanism=joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    gas=CreateSolution(mechanism)
    data=MultiTransportData(mechanism*".multicomponent.npz",gas)
    X="H2:1.1,O2:1,AR:5"
    @test_throws ArgumentError FreeFlame(gas;X,discretization=:invalid)
    @test FreeFlame(gas;X).discretization == :conservative
    @test BurnerFlame(gas;X,mdot=.03).discretization == :conservative
    @test Arrhenius._flame_face_centering(0.) == 1
    @test 0 < Arrhenius._flame_face_centering(10.) < Arrhenius._flame_face_centering(.1) < 1
    @test Arrhenius._flame_face_centering(1e8) ≈ 2e-8 atol=1e-15
    atomic=gas.ele_matrix'\gas.MW
    E=Diagonal(atomic)*gas.ele_matrix*Diagonal(1 ./gas.MW)
    for (model,basis,soret) in ((:mixture_averaged,:mole,false),(:mixture_averaged,:mass,false),
            (:mixture_averaged,:mole,true),(:mixture_averaged,:mass,true),
            (:multicomponent,:mole,false),(:multicomponent,:mole,true))
        f=FreeFlame(gas;X,width=.06,discretization=:conservative,
            transport_model=model,flux_gradient_basis=basis,soret,multicomponent_data=data)
        solve!(f)
        w=Arrhenius.FlameWorkspace(f)
        r=flame_residual!(similar(f.state),f,f.state,w)
        c=w.conservative;mdot=f.state[end,1]
        @test f.converged && norm(r,Inf)<1e-8
        @test minimum(mass_fractions(f)) > -1e-12
        @test maximum(abs,E*(mass_fractions(f)[:,end]-f.inlet_Y)) < 1e-6
        @test maximum(abs,E*c.species_flux .- mdot.*(E*f.inlet_Y))/mdot < 1e-6
        @test maximum(abs,sum(w.flux;dims=1))/mdot < 1e-12
        @test (maximum(c.enthalpy_flux)-minimum(c.enthalpy_flux))/maximum(abs,c.enthalpy_flux) < 1e-5
        # A sufficiently long adiabatic domain approaches the independently
        # computed constant-enthalpy equilibrium temperature for every transport.
        @test temperature(f)[end] ≈ equilibrate(gas;T=300.,P=one_atm,X,mode=:HP).T atol=.05
        if model==:mixture_averaged && basis==:mole && !soret
            previous=copy(f.state);previous[2,end]-=1e-5;previous[4,end]+=1e-5
            dt=1e-3
            rp=flame_residual!(similar(r),f,f.state,w;previous,dt)
            hprevious=cal_h_RT(gas,1000previous[1,end],f.pressure,mole_fractions(gas,X)).*R.*(1000previous[1,end])
            oldh=dot(previous[2:end-1,end],hprevious./gas.MW)
            delta=-Arrhenius._flame_timescale/(1000*c.cp[end]*dt)*(c.enthalpy[end]-oldh)
            @test rp[1,end]-r[1,end] ≈ delta rtol=1e-10 atol=1e-13
            @test abs(delta)>1e-8 # Composition storage contributes with no T change.
            # The final half-cell retains the same nearest-neighbor Jacobian
            # stencil, including its energy/species pseudo-time coupling.
            band=copy(Arrhenius._flame_jacobian(f,f.state,w,rp;previous,dt,analytic=false))
            B,N=size(f.state);kl=2B-1
            for (k,j) in ((1,N),(2,N-1),(B,1))
                u=copy(f.state);step=1e-7*max(abs(u[k,j]),k==1 ? .1 : 1e-5);u[k,j]+=step
                flame_residual!(similar(r),f,f.state,w;previous,dt)
                perturbed=flame_residual!(similar(r),f,u,w;previous,dt,update_transport=false)
                column=(j-1)*B+k;from_band=zeros(length(u))
                for row in max(1,column-kl):min(length(u),column+kl)
                    from_band[row]=band[2kl+1+row-column,column]
                end
                @test from_band ≈ vec((perturbed-rp)/step) rtol=1e-8 atol=1e-9
            end
        end
    end
    f=BurnerFlame(gas;T=373.,P=.05one_atm,mdot=.06,X="H2:1.5,O2:1,AR:7",width=.5,discretization=:conservative)
    set_temperature_profile!(f,[0.,.005,.01,.02,.05,.1,1.],[373.,650.,1000.,1350.,1650.,1750.,1750.])
    solve!(f;slope=.05,curve=.1)
    @test temperature(f) ≈ f.imposed_temperature atol=1e-10
    @test maximum(abs,E*(mass_fractions(f)[:,end]-f.inlet_Y)) < 1e-6
    @test norm(flame_residual!(similar(f.state),f),Inf)<1e-8

    kink=BurnerFlame(gas;T=373.,P=.05one_atm,mdot=.03,X="H2:1.5,O2:1,AR:7",
        grid=collect(0.:.001:.02),discretization=:conservative,soret=true,multicomponent_data=data)
    set_temperature_profile!(kink,[0.,.5,1.],[373.,1000.,1500.])
    @test all(z->z in kink.grid,kink.profile_positions)
    hspecies=findfirst(==("H"),gas.species_names)+1
    inert=findfirst(==("AR"),gas.species_names)+1
    for j in eachindex(kink.grid)
        kink.state[2:end-1,j].=kink.inlet_Y
        y=.001*abs(kink.grid[j]-.01)
        kink.state[hspecies,j]+=y;kink.state[inert,j]-=y
    end
    # A genuine derivative jump at the prescribed knot must not force an
    # endless sequence of smaller cells. Slope and ratio remain active.
    @test !Arrhenius._refine_flame!(kink;slope=1.,curve=.1,ratio=3.)
    for j in eachindex(kink.grid)
        y=1e-5*exp(-((kink.grid[j]-.012)/.002)^2)
        kink.state[hspecies,j]+=y;kink.state[inert,j]-=y
    end
    oldgrid=copy(kink.grid)
    @test Arrhenius._refine_flame!(kink;slope=1.,curve=.1,ratio=3.)
    @test any(z->.011<z<.014 && !(z in oldgrid),kink.grid)
end
