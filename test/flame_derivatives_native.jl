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
        f=kind==:free ? FreeFlame(gas;X="H2:1.1,O2:1,AR:5",grid,discretization=:conservative) :
            BurnerFlame(gas;X="H2:1.1,O2:1,AR:5",grid,mdot=.07,discretization=:conservative)
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
    f=FreeFlame(gas;X="H2:1.1,O2:1,AR:5",grid=collect(range(0.,.01;length=5)),discretization=:conservative)
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
    free=FreeFlame(gas;X="H2:1.1,O2:1,AR:5",grid,discretization=:conservative)
    burner=BurnerFlame(gas;X="H2:1.1,O2:1,AR:5",grid,mdot=.07,discretization=:conservative)
    imposed=BurnerFlame(gas;X="H2:1.1,O2:1,AR:5",grid,mdot=.07,discretization=:conservative)
    set_temperature_profile!(imposed,grid,temperature(imposed);relative=false)
    @test Arrhenius._flame_correction_enabled(free)
    @test Arrhenius._flame_correction_enabled(burner)
    @test !Arrhenius._flame_correction_enabled(imposed)
    legacy=FreeFlame(gas;X="H2:1.1,O2:1,AR:5",grid,discretization=:finite_difference)
    @test !Arrhenius._flame_correction_enabled(legacy)
end


@testset "auto finite-difference flame interior species Jacobian" begin
    mechanism=joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    gas=CreateSolution(mechanism)
    grid=collect(range(0.,.01;length=5))
    f=FreeFlame(gas;X="H2:1.1,O2:1,AR:5",grid,discretization=:finite_difference)
    n=gas.n_species;B=n+2
    for j in eachindex(grid)
        f.state[1,j]=.7+.08j
        f.state[2:n+1,j].=.01 .+ .001 .* (1:n) .+ .0002j
        f.state[2:n+1,j]./=sum(@view(f.state[2:n+1,j]))
        f.state[end,j]=.08+.003j
    end
    previous=.997 .* f.state
    dt=.004
    w=Arrhenius.FlameWorkspace(f)
    r=flame_residual!(similar(f.state),f,f.state,w;previous,dt)
    stats=solver_stats(f)
    node_start=stats.property_nodes; residual_start=stats.residual_evals
    band=copy(Arrhenius._flame_jacobian(f,f.state,w,r;previous,dt,analytic=true))
    optimized_nodes=stats.property_nodes-node_start
    optimized_calls=stats.residual_evals-residual_start
    node_start=stats.property_nodes; residual_start=stats.residual_evals
    legacy=copy(Arrhenius._flame_jacobian(f,f.state,w,r;previous,dt,analytic=true,local_property_updates=false))
    legacy_nodes=stats.property_nodes-node_start
    legacy_calls=stats.residual_evals-residual_start
    @test band==legacy
    @test optimized_nodes==length(grid)+2n
    @test legacy_nodes>optimized_nodes
    @test optimized_calls==legacy_calls
    node_start=stats.property_nodes; residual_start=stats.residual_evals
    fd=copy(Arrhenius._flame_jacobian(f,f.state,w,r;previous,dt,analytic=false))
    fd_nodes=stats.property_nodes-node_start
    fd_calls=stats.residual_evals-residual_start
    node_start=stats.property_nodes; residual_start=stats.residual_evals
    fd_legacy=copy(Arrhenius._flame_jacobian(f,f.state,w,r;previous,dt,analytic=false,local_property_updates=false))
    fd_legacy_nodes=stats.property_nodes-node_start
    @test fd==fd_legacy
    @test fd_nodes==(B-1)*length(grid)
    @test fd_legacy_nodes==B*length(grid)
    @test fd_calls==stats.residual_evals-residual_start
    J=zeros(B*length(grid),B*length(grid))
    for col in 1:B*length(grid),row in max(1,col-2B+1):min(B*length(grid),col+2B-1)
        J[row,col]=band[4B-1+row-col,col]
    end
    Jad=flame_test_automatic_jacobian(f,w;previous,dt)
    maxerr=0.0;worst=(0,0,0,0,0.0,0.0)
    family=zeros(3)
    for point in 2:length(grid)-1,component in 2:B-1
        col=(point-1)*B+component
        for rowpoint in 2:length(grid)-1, rowcomponent in 1:B
            row=(rowpoint-1)*B+rowcomponent
            err=abs(J[row,col]-Jad[row,col])/max(abs(Jad[row,col]),1e-8)
            fam=rowcomponent==1 ? 1 : rowcomponent==B ? 3 : 2
            family[fam]=max(family[fam],err)
            if err>maxerr
                maxerr=err;worst=(rowcomponent,rowpoint,point,component,J[row,col],Jad[row,col])
            end
        end
    end
    @info "analytic Jacobian derivative diagnostic" maxerr worst family
    @test maxerr < 2e-7

end

@testset "auto analytic burner interior species Jacobian" begin
    mechanism=joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    gas=CreateSolution(mechanism)
    grid=collect(range(0.,.01;length=5))
    for prescribed_temperature in (false,true)
        f=BurnerFlame(gas;X="H2:1.1,O2:1,AR:5",grid,mdot=.08,
            discretization=:finite_difference)
        prescribed_temperature && set_temperature_profile!(f,[grid[1],grid[end]],[f.inlet_temperature,900.];relative=false)
        n=gas.n_species;B=n+2
        for j in eachindex(grid)
            f.state[1,j]=.7+.04j
            f.state[2:n+1,j].=.01 .+ .001 .* (1:n) .+ .0002j
            f.state[2:n+1,j]./=sum(@view(f.state[2:n+1,j]))
            f.state[end,j]=.08
        end
        w=Arrhenius.FlameWorkspace(f)
        r=flame_residual!(similar(f.state),f,f.state,w)
        stats=solver_stats(f); node_start=stats.property_nodes
        analytic=copy(Arrhenius._flame_jacobian(f,f.state,w,r;analytic=true))
        optimized_nodes=stats.property_nodes-node_start
        legacy=copy(Arrhenius._flame_jacobian(f,f.state,w,r;analytic=true,local_property_updates=false))
        legacy_nodes=stats.property_nodes-node_start-optimized_nodes
        @test analytic==legacy
        @test optimized_nodes==length(grid)+4n
        @test legacy_nodes>optimized_nodes
        automatic=flame_test_automatic_jacobian(f,w)
        kl=2B-1
        worst=0.0;worst_entry=(0,0,0.0,0.0)
        for point in 2:length(grid)-1,component in 2:B-1
            col=(point-1)*B+component
            for rowpoint in max(2,point-1):min(length(grid)-1,point+1),rowcomponent in 1:B
                row=(rowpoint-1)*B+rowcomponent
                a=analytic[2kl+1+row-col,col]
                d=automatic[row,col]
                error=abs(a-d)/max(abs(d),1e-6)
                if error>worst; worst=error; worst_entry=(row,col,a,d); end
            end
        end
        @info "burner analytic species Jacobian versus ForwardDiff" prescribed_temperature worst worst_entry
        @test worst < 2e-4
        # The native burner inlet balance includes the face flux, so keep
        # the adjacent species-column boundary entry finite-differenced.
        col=B+2; inlet_row=2
        @test analytic[2kl+1+inlet_row-col,col] ≈ automatic[inlet_row,col] rtol=2e-4
        @test automatic[inlet_row,col] ≈ 0.02936003949483076 rtol=1e-8
        if prescribed_temperature
            for point in 2:length(grid)-1,component in 2:B-1
                col=(point-1)*B+component
                for rowpoint in 2:length(grid)-1
                    row=(rowpoint-1)*B+1
                    @test analytic[2kl+1+row-col,col] == 0.0
                end
            end
        end
    end
end

@testset "Cantera auto Jacobian fixture: free flame source species columns" begin
    root=joinpath(@__DIR__,"data","cantera")
    ref=Arrhenius.npzread(joinpath(root,"free_auto_stagni.npz"))
    expected=Arrhenius.npzread(joinpath(root,"free_auto_stagni_source.npz"))["source_jacobian"]
    mechanism=joinpath(@__DIR__,"..","mechanism","stagni_2023.yaml")
    gas=CreateSolution(mechanism); grid=vec(ref["grid"])
    f=FreeFlame(gas;T=300.0,P=one_atm,X="NH3:0.7,H2:0.3,O2:0.5625,N2:2.115",
        grid,width=grid[end]-grid[1],transport_model=:mixture_averaged,
        discretization=:finite_difference)
    n,N=gas.n_species,length(grid)
    f.anchor=4
    f.state[1,:].=vec(ref["values/T"])./1000
    for (k,name) in enumerate(gas.species_names)
        f.state[k+1,:].=vec(ref["values/"*name])
    end
    for j in 1:N
        rho=Arrhenius._flame_density(f,1000f.state[1,j],@view(f.state[2:n+1,j])).density
        f.state[end,j]=rho*ref["values/velocity"][j]
    end
    w=Arrhenius.FlameWorkspace(f)
    r=similar(f.state); flame_residual!(r,f,f.state,w)
    boundary_jac=Arrhenius._flame_jacobian(f,f.state,w,r;analytic=true)
    @test Arrhenius._flame_cantera_species_perturbation(f.state,2,1) ==
        1e-5*abs(f.state[2,1])+1e-10
    B=n+2; halfband=2B-1; h_index=findfirst(==("H"),gas.species_names)
    left_row=B; left_col=B+h_index+1
    right_row=(N-1)*B+B; right_col=(N-2)*B+h_index+1
    left_entry=boundary_jac[2halfband+1+left_row-left_col,left_col]
    right_entry=boundary_jac[2halfband+1+right_row-right_col,right_col]
    @test left_entry ≈ 24.0162922725 rtol=2e-8
    @test right_entry ≈ -23.1400344324 rtol=2e-8
    # Full Cantera auto matrix transformed to native (T/1000,Y,m) coordinates.
    # Regenerate with validation/generate_free_auto_jacobian_fixture.py.
    expected_full=Arrhenius.npzread(joinpath(root,"free_auto_stagni_jacobian_reference.npz"))["jacobian"]
    dense_jac=zeros(N*B,N*B)
    for col in 1:N*B, row in max(1,col-halfband):min(N*B,col+halfband)
        dense_jac[row,col]=boundary_jac[2halfband+1+row-col,col]
    end
    full_jacobian_relative_error=norm(dense_jac-expected_full)/norm(expected_full)
    @info "Cantera full auto Jacobian fixture comparison" full_jacobian_relative_error
    for (name,cols) in (("interior species",[j*B+k for j in 1:N-2 for k in 2:B-1]),
            ("boundary species",[j*B+k for j in (0,N-1) for k in 2:B-1]),
            ("temperature",[j*B+1 for j in 0:N-1]),("mass flux",[j*B+B for j in 0:N-1]))
        rel=norm(dense_jac[:,cols]-expected_full[:,cols])/norm(expected_full[:,cols])
        @info "Cantera Jacobian block" name rel
    end
    @test full_jacobian_relative_error < 2e-10
    d=Arrhenius._ConservativeFlameDerivatives(gas,N)
    G=zeros(n,n+1)
    maxrel=0.0
    for j in 2:N-1
        ok=Arrhenius._flame_source_jacobian!(G,d,gas,f.pressure,1000f.state[1,j],
            @view(f.state[2:n+1,j]),@view(w.X[:,j]),@view(w.h[:,j]),w.rate_caches[j])
        @test ok
        diff=G[:,2:end]-expected[:,:,j]
        maxrel=max(maxrel,norm(diff)/max(norm(expected[:,:,j]),1e-300))
    end
    @info "Cantera source derivative fixture comparison" maxrel
    @test maxrel < 1e-12
end
@testset "fixed-temperature finite-difference Newton merit" begin
    mechanism=joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    gas=CreateSolution(mechanism)
    burner=BurnerFlame(gas;X="H2:1.1,O2:1,AR:5",mdot=.07,
        grid=collect(range(0.,.03;length=5)),discretization=:finite_difference)
    set_temperature_profile!(burner,[0.,.25,1.],[300.,1600.,300.])
    solve!(burner;refine_grid=false)
    stats=solver_stats(burner)
    @info "fixed_t_h2o2_counts" jacobians=stats.jacobians newton_iterations=stats.newton_iterations
    @test stats.jacobians == 6
    @test stats.newton_iterations == 110
    # Same discrete 5-node state from origin/main ad987042, used to guard solution changes.
    expected_values=parse.(Float64,split(strip(read(joinpath(@__DIR__,"data","h2o2_fixed_t_main_state.csv"),String)),","))
    expected_state=reshape(expected_values,size(burner.state))
    state_difference=maximum(abs.(burner.state .- expected_state))
    @info "fixed_t_h2o2_state_comparison" max_abs_difference=state_difference
    @test state_difference <= 1e-7
    @test burner.converged
    @test maximum(heat_release_rate(burner)) > 0
    @test stats.jacobians > 0 && stats.newton_iterations >= stats.jacobians
    @test norm(flame_residual!(similar(burner.state),burner),Inf) < 1e-8
    @test maximum(abs.(sum(mass_fractions(burner);dims=1) .- 1)) < 1e-10
    @test minimum(mass_fractions(burner)) > -1e-12
end
