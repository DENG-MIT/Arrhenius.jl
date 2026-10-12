using Arrhenius, ForwardDiff, LinearAlgebra, Test

function flame_test_dual_fields(value,::Type{D}) where D
    names=propertynames(value)
    return NamedTuple{names}(map(names) do name
        field=getproperty(value,name)
        field isa Array{Float64} ? D.(field) : field
    end)
end

function flame_test_automatic_jacobian(f,w;previous=nothing,dt=Inf)
    shape=size(f.state)
    return ForwardDiff.jacobian(vec(f.state)) do inputs
        D=eltype(inputs)
        names=propertynames(w)
        storage=NamedTuple{names}(map(names) do name
            field=getproperty(w,name)
            name==:kinetics ? KineticsWorkspace(f.gas.reaction,D) :
                name==:conservative ? flame_test_dual_fields(field,D) :
                field isa Array{Float64} ? D.(field) : field
        end)
        state=reshape(inputs,shape)
        vec(flame_residual!(similar(state),f,state,storage;
            previous,dt,update_transport=false))
    end
end

@testset "trial correction metric finite trace and units" begin
    state=zeros(5,5); state[1,:].=1.2; state[2,:].=.1; state[3,:].=-1e-12; state[5,:].=.07
    initial=copy(state); step=fill(.003,size(state))
    for transient in (false,true)
        weights=Arrhenius._flame_correction_weights(state,transient)
        atol=transient ? 1e-11 : 1e-9
        physical=copy(state); physical[1,:].*=1000
        pstep=copy(step); pstep[1,:].*=1000
        expected=sqrt(sum((pstep[k,j]/(1e-4*sum(abs,physical[k,:])/5+atol))^2 for k in 1:5,j in 1:5)/25)
        @test Arrhenius._flame_correction_norm(step,weights) ≈ expected rtol=1e-14
        @test all(isfinite,weights) && all(>(0),weights) && isequal(state,initial)
    end
end

@testset "conservative flame local Jacobian" begin
    mechanism=joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    gas=CreateSolution(mechanism)
    data=MultiTransportData(mechanism*".multicomponent.npz",gas)
    for (kind,model,basis,soret,branch,pseudo) in (
            (:free,:mixture_averaged,:mole,false,:positive,false),
            (:free,:mixture_averaged,:mole,false,:zero,true),
            (:fixed,:mixture_averaged,:mass,true,:negative,true),
            (:fixed,:mixture_averaged,:mass,true,:mixed,true),
            (:free,:multicomponent,:mole,true,:positive,true),
            (:free,:multicomponent,:mole,true,:mixed,true))
        grid=collect(range(0.,.03;length=5))
        f=kind==:free ? FreeFlame(gas;X="H2:1.1,O2:1,AR:5",grid) :
            BurnerFlame(gas;X="H2:1.1,O2:1,AR:5",grid,mdot=.07)
        set_transport!(f,model;data,soret,flux_gradient_basis=basis)
        n,N=gas.n_species,length(grid);B=n+2
        for j in 1:N
            f.state[1,j]=.53+.21j
            for k in 1:n
                f.state[k+1,j]=.02+.003k+.002sin(k+j)
            end
            f.state[2:n+1,j]./=sum(@view(f.state[2:n+1,j]))
            f.state[end,j]=.07+.004j
        end
        trace=findfirst(==("H2O2"),gas.species_names)+1
        branch==:zero && (f.state[trace,:].=0)
        branch==:negative && (f.state[trace,:].=-2e-8)
        branch==:mixed && (f.state[trace,1:2:end].=-2e-8)
        if kind==:fixed
            f.state[1,1]=f.inlet_temperature/1000
            set_temperature_profile!(f,grid,temperature(f);relative=false)
        end
        previous=pseudo ? .99 .* f.state : nothing;dt=pseudo ? .003 : Inf
        w=Arrhenius.FlameWorkspace(f)
        r=flame_residual!(similar(f.state),f,f.state,w;previous,dt)
        @test Arrhenius._conservative_flame_jacobian!(f,f.state,w,r;previous,dt)
        J=zeros(B*N,B*N)
        for col in 1:B*N,row in max(1,col-2B+1):min(B*N,col+2B-1)
            J[row,col]=w.band[4B-1+row-col,col]
        end
        original=copy(f.state)
        # The Jacobian uses the positive one-sided derivative at max(Y,0)'s
        # kink. AD at an exact tie selects a different branch; use its limit.
        branch==:zero && (f.state[trace,:].=1e-20)
        Jad=flame_test_automatic_jacobian(f,w;previous,dt)
        f.state.=original
        scale=max.(abs.(Jad),maximum(abs,Jad;dims=2).*1e-8,1e-12)
        @test maximum(abs.(J-Jad)./scale)<3e-6
        direction=reshape([.03sin(k) for k in eachindex(f.state)],size(f.state))
        branch==:zero && (direction[trace,:].=abs.(@view(direction[trace,:])))
        branch==:negative && (direction[trace,:].=0)
        branch==:mixed && (direction[trace,1:2:end].=0)
        h=1e-5
        r1=flame_residual!(similar(r),f,original+h.*direction,w;
            previous,dt,update_transport=false)
        r2=flame_residual!(similar(r),f,original+(h/2).*direction,w;
            previous,dt,update_transport=false)
        numerical=vec((4r2-r1-3r)/h)
        exact=J*vec(direction)
        @test maximum(abs.(exact-numerical)./max.(abs.(exact),1e-6))<1e-4
    end
end

@testset "multicomponent trace flux cancellation" begin
    mechanism=joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    gas=CreateSolution(mechanism)
    data=MultiTransportData(mechanism*".multicomponent.npz",gas)
    f=FreeFlame(gas;X="H2:1.1,O2:1,AR:5",grid=collect(range(0.,.01;length=5)))
    set_transport!(f,:multicomponent;data)
    k=findfirst(==("H2O2"),gas.species_names)
    q=findfirst(==("H2"),gas.species_names)
    oxygen=findfirst(==("O2"),gas.species_names)
    argon=findfirst(==("AR"),gas.species_names)
    n=gas.n_species
    for j in axes(f.state,2)
        f.state[1,j]=1.2;f.state[end,j]=0. # isolate diffusion from advection
        f.state[2:n+1,j].=1e-22
        f.state[k+1,j]=(1+j)*1e-20
        f.state[q+1,j]=.12+.01j;f.state[oxygen+1,j]=.25-.008j
        f.state[argon+1,j]=1-sum(@view(f.state[2:n+1,j]))+f.state[argon+1,j]
    end
    w=Arrhenius.FlameWorkspace(f)
    Arrhenius._flame_properties!(w,f,f.state)
    transport=MultiTransportWorkspace(data)
    setprecision(BigFloat,256) do
        mw=BigFloat.(gas.MW)
        for j in 1:4
            # Assemble the original Float64 physical prefactors independently;
            # the oracle never reads the centered cache under test.
            y=0.5 .*(f.state[2:n+1,j]+f.state[2:n+1,j+1])
            ysum=0.;inverse=0.
            for l in 1:n;ysum+=y[l];inverse+=y[l]/gas.MW[l];end
            meanMW=ysum/inverse
            x=[(y[l]/ysum)*meanMW/gas.MW[l] for l in 1:n]
            T=500*(f.state[1,j]+f.state[1,j+1])
            multicomponent_transport!(transport,data,gas,f.pressure,T,x)
            p=BigFloat[f.pressure/(Arrhenius.R*T*meanMW)*gas.MW[k]*gas.MW[l]*
                transport.diffusion[k,l] for l in 1:n]
            left=BigFloat.(f.state[2:n+1,j])./mw;left./=sum(left)
            right=BigFloat.(f.state[2:n+1,j+1])./mw
            denominator=sum(right);right./=denominator
            # dX/dY_H2 at fixed other mass fractions, with normalization exact
            # in the high-precision oracle. Transport coefficients stay frozen.
            derivative=-right/(mw[q]*denominator)
            derivative[q]+=1/(mw[q]*denominator)
            dz=BigFloat(f.grid[j+1]-f.grid[j])
            expected_flux=Float64(dot(p,right-left)/dz)
            expected_derivative=Float64(dot(p,derivative)/dz)
            Arrhenius._flame_face_derivatives!(w.conservative.derivatives,f,f.state,w,j)
            @test w.flux[k,j]≈expected_flux rtol=1e-10 atol=1e-30
            @test w.conservative.derivatives.flux_right[k+1,q+1]≈expected_derivative rtol=1e-10 atol=1e-30
        end
    end
end

@testset "correction merit selector" begin
    mechanism=joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    gas=CreateSolution(mechanism)
    grid=collect(range(0.,.03;length=5))
    free=FreeFlame(gas;X="H2:1.1,O2:1,AR:5",grid)
    burner=BurnerFlame(gas;X="H2:1.1,O2:1,AR:5",grid,mdot=.07)
    imposed=BurnerFlame(gas;X="H2:1.1,O2:1,AR:5",grid,mdot=.07)
    set_temperature_profile!(imposed,grid,temperature(imposed);relative=false)
    @test Arrhenius._flame_correction_enabled(free)
    @test Arrhenius._flame_correction_enabled(burner)
    @test !Arrhenius._flame_correction_enabled(imposed)
    legacy=FreeFlame(gas;X="H2:1.1,O2:1,AR:5",grid,discretization=:finite_difference)
    @test !Arrhenius._flame_correction_enabled(legacy)
end
