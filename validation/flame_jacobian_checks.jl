using Arrhenius,ForwardDiff,LinearAlgebra,NPZ,Test
const AF=Arrhenius
BLAS.set_num_threads(1)
parameters,output=ARGS[1:2];mkpath(output)

function dense_band(band,B,N)
    J=zeros(B*N,B*N)
    for column in 1:B*N,row in max(1,column-2B+1):min(B*N,column+2B-1)
        J[row,column]=band[4B-1+row-column,column]
    end
    return J
end
function dual_array_fields(value,::Type{D}) where D
    names=propertynames(value)
    return NamedTuple{names}(map(names) do name
        field=getproperty(value,name)
        field isa Array{Float64} ? D.(field) : field
    end)
end
function dual_workspace(w,gas,::Type{D}) where D
    names=propertynames(w)
    return NamedTuple{names}(map(names) do name
        value=getproperty(w,name)
        name==:kinetics ? KineticsWorkspace(gas.reaction,D) :
            name==:conservative ? dual_array_fields(value,D) :
            value isa Array{Float64} ? D.(value) : value
    end)
end
function automatic_jacobian(f,w;previous=nothing,dt=Inf)
    shape=size(f.state)
    rhs=inputs->begin
        state=reshape(inputs,shape);storage=dual_workspace(w,f.gas,eltype(inputs))
        residual=similar(state)
        AF.flame_residual!(residual,f,state,storage;previous,dt,update_transport=false)
        vec(residual)
    end
    return ForwardDiff.jacobian(rhs,vec(f.state))
end
function check_case(gas,data,mode,branch,pseudo;kind=:free)
    grid=collect(range(0.,.03;length=5))
    f=kind==:free ? FreeFlame(gas;T=300.,P=one_atm,X="H2:1.1,O2:1,AR:5",grid) :
        BurnerFlame(gas;T=300.,P=one_atm,X="H2:1.1,O2:1,AR:5",grid,mdot=.07)
    multi=startswith(mode,"multi");soret=endswith(mode,"soret")
    set_transport!(f,multi ? :multicomponent : :mixture_averaged;data,soret,
        flux_gradient_basis=startswith(mode,"mass") ? :mass : :mole)
    n,N=gas.n_species,length(f.grid)
    for j in 1:N
        f.state[1,j]=.53+.21*j
        for k in 1:n
            f.state[k+1,j]=.02+.003*k+.002*sin(k+j)
        end
        f.state[2:n+1,j]./=sum(@view(f.state[2:n+1,j]))
        f.state[end,j]=.07+.004*j
    end
    if branch!=:positive
        k=findfirst(==("H2O2"),gas.species_names)
        f.state[k+1,:].=branch==:zero ? 0.0 : -2e-8
    end
    if kind==:fixed
        f.state[1,1]=f.inlet_temperature/1000
        set_temperature_profile!(f,f.grid,temperature(f);relative=false)
    end
    previous=pseudo ? copy(f.state).*.99 : nothing;dt=pseudo ? .003 : Inf
    w=AF.FlameWorkspace(f);r=similar(f.state)
    AF.flame_residual!(r,f,f.state,w;previous,dt)
    ok=AF._conservative_flame_jacobian!(f,f.state,w,r;previous,dt)
    @test ok
    J=dense_band(w.band,n+2,N)
    original=copy(f.state)
    if branch==:zero
        # max(y,0) has two directional derivatives at its kink. Preserve the
        # positive-side convention of the production forward difference;
        # ForwardDiff's exact-tie convention instead selects the clipped side.
        f.state[findfirst(==("H2O2"),gas.species_names)+1,:].=1e-20
    end
    Jad=automatic_jacobian(f,w;previous,dt)
    f.state.=original
    rowscale=maximum(abs,Jad;dims=2)
    error=maximum(abs.(J-Jad)./max.(abs.(Jad),rowscale.*1e-8,1e-12))
    @test error<3e-6
    wf=AF.FlameWorkspace(f);AF.flame_residual!(r,f,f.state,wf;previous,dt)
    Jfd=dense_band(AF._flame_jacobian(f,f.state,wf,r;previous,dt,analytic=false),n+2,N)
    fd_error=maximum(abs.(J-Jfd)./max.(abs.(J),rowscale.*1e-4,1e-10))
    @test all(isfinite,Jfd)
    # Independently differentiate the unchanged residual with a second-order
    # one-sided Richardson stencil. The old per-species 1e-12 step can lose
    # flux-derivative digits on rows with large chemical residuals.
    direction=reshape([.03*sin(k) for k in eachindex(f.state)],size(f.state))
    if branch!=:positive
        index=findfirst(==("H2O2"),gas.species_names)+1
        direction[index,:].=branch==:zero ? abs.(direction[index,:]) : 0.0
    end
    reference=AF.FlameWorkspace(f);base=similar(f.state)
    AF.flame_residual!(base,f,f.state,reference;previous,dt)
    deltas=Vector{Float64}[]
    for step in [1e-5,5e-6]
        trial=f.state.+step.*direction;changed=similar(trial)
        AF.flame_residual!(changed,f,trial,reference;previous,dt,update_transport=false)
        push!(deltas,vec(changed.-base)./step)
    end
    richardson=2deltas[2]-deltas[1];prediction=J*vec(direction)
    directional_error=maximum(abs.(richardson.-prediction)./max.(abs.(prediction),vec(rowscale).*1e-6,1e-10))
    @test directional_error<1e-4
    println((;species=gas.n_species,kind,mode,branch,pseudo,ad_error=error,fd_error,directional_error));flush(stdout)
    npzwrite(joinpath(output,"$(gas.n_species)-$kind-$mode-$branch-$pseudo.npz"),Dict("analytic"=>J,"automatic"=>Jad,
        "finite_difference"=>Jfd,"direction"=>direction,"richardson"=>richardson,"predicted_direction"=>prediction))
end
mechanism=joinpath(parameters,"h2o2.yaml")
gas=CreateSolution(mechanism);data=MultiTransportData(mechanism*".multicomponent.npz",gas;mechanism)
@testset "Conservative premixed local Jacobian" begin
    for mode in ["mole","mass","mole-soret","mass-soret","multi","multi-soret"],
            branch in [:positive,:zero,:negative],pseudo in [false,true]
        check_case(gas,data,mode,branch,pseudo)
    end
    for kind in [:burner,:fixed],mode in ["mole","multi-soret"],branch in [:zero,:negative],pseudo in [false,true]
        check_case(gas,data,mode,branch,pseudo;kind)
    end
    mechanism=joinpath(parameters,"gri30.yaml")
    methane=CreateSolution(mechanism);methane_data=MultiTransportData(mechanism*".multicomponent.npz",methane;mechanism)
    for mode in ["mole","multi-soret"],branch in [:positive,:zero,:negative]
        check_case(methane,methane_data,mode,branch,true;kind=:fixed)
    end
end
