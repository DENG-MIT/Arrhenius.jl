using Test, Arrhenius, NPZ, LinearAlgebra

@testset "counterflow physical growth residual" begin
    ref=npzread(joinpath(@__DIR__,"data","cantera_counterflow_growth_norm.npz"))
    gas=CreateSolution(joinpath(@__DIR__,"..","mechanism","h2o2.yaml"))
    @test gas.species_names==split(String(vec(ref["species_names_utf8"])),"\n")
    f=CounterflowDiffusionFlame(gas;fuel="H2:1",oxidizer="O2:1",
        T_fuel=800.,T_oxidizer=800.,P=ref["pressure"][1],
        mdot_fuel=ref["mdot_fuel"][1],mdot_oxidizer=ref["mdot_oxidizer"][1],grid=ref["grid"])
    n=gas.n_species
    function native_state(values)
        T=values[1,:]; Y=values[2:n+1,:]
        rho=f.pressure ./ (Arrhenius.R .* T .* vec(sum(Y./gas.MW;dims=1)))
        u=copy(values)
        u[1,:]./=1000
        u[end-2,:]./=Arrhenius._counterflow_Vscale
        u[end-1,:]./=Arrhenius._counterflow_Lscale
        u[end,:].*=rho
        return u
    end
    before=native_state(ref["before_values"]); after=native_state(ref["after_values"])
    f.state=copy(before); w=CounterflowWorkspace(f); r=similar(before)
    physical!(r,x,p,dt)=Arrhenius._counterflow_growth_residual!(r,f,x,w;previous=p,dt)
    raw!(r,x,p,dt)=counterflow_residual!(r,f,x,w;previous=p,dt)
    for (state,expected) in ((before,ref["before_residual"]),(after,ref["after_residual"]))
        physical!(r,state,before,1e-5)
        # Only Lambda changes between states, so transient terms vanish.
        for k in axes(r,1)
            @test maximum(abs.(r[k,:].-expected[k,:]))<=1e-9*max(1,maximum(abs.(expected[k,:])))
        end
    end
    g=Arrhenius._timestep_growth(physical!,r,:residual_ratio,1.5,before,after,1e-5,1)
    @test g≈(maximum(abs.(ref["before_residual"]))/maximum(abs.(ref["after_residual"])))^.2 rtol=1e-10
    # The old scaled norm incorrectly selected the maximum allowed growth.
    @test Arrhenius._timestep_growth(raw!,r,:residual_ratio,1.5,before,after,1e-5,1)==1.5
    @test Arrhenius._timestep_growth(physical!,r,:residual_ratio,1.5,after,before,1e-5,1)==1
end


@testset "growth residual dispatch by counterflow kind" begin
    gas=CreateSolution(joinpath(@__DIR__,"..","mechanism","h2o2.yaml"))
    grid=collect(range(0.,.02;length=7))
    kinds=(
        CounterflowDiffusionFlame(gas;fuel="H2:1",oxidizer="O2:1",mdot_fuel=.24,mdot_oxidizer=.72,grid),
        CounterflowPremixedFlame(gas;reactants="H2:2,O2:1,AR:7",mdot_reactants=.12,mdot_products=.06,width=.02,grid),
        CounterflowTwinPremixedFlame(gas;reactants="H2:2,O2:1,AR:7",mdot=.12,width=.02,grid),
        ImpingingJet(gas;reactants="H2:2,O2:1,AR:7",mdot=.06,width=.02,grid))
    expected=(CounterflowDiffusionFlame,Arrhenius.AbstractPremixedCounterflow,
        CounterflowTwinPremixedFlame,ImpingingJet)
    for (f,T) in zip(kinds,expected)
        w=CounterflowWorkspace(f); r=similar(f.state); x=copy(f.state)
        method=which(Arrhenius._counterflow_growth_residual!,
            (typeof(r),typeof(f),typeof(x),typeof(w)))
        @test Base.unwrap_unionall(method.sig).parameters[3] === T
        if f isa CounterflowDiffusionFlame
            @test Arrhenius._counterflow_growth_residual!(r,f,x,w) === r
            @test Arrhenius._check_timestep_controls(1.5,:residual_ratio,0,:auto) === nothing
        else
            f isa ImpingingJet && @test_throws ArgumentError Arrhenius._counterflow_growth_residual!(r,f,x,w)
            err=try
                solve!(f;time_step_growth_strategy=:residual_ratio)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin(string(nameof(typeof(f))),sprint(showerror,err))
        end
    end
end

@testset "counterflow dt trajectory decisions" begin
    # Cantera's own accepted-step trajectory on the 15 MPa h2o2 diffusion
    # flame under residual-ratio growth: for every step whose attempt
    # evaluated no Jacobian, the native growth function evaluated on
    # Cantera's recorded before/after states and attempt size must return
    # the factor Cantera applied (Cantera gates growth on the Jacobian
    # evaluation count, so gated steps keep dt unchanged by construction).
    ref=npzread(joinpath(@__DIR__,"data","cantera_counterflow_dt_trajectory.npz"))
    gas=CreateSolution(joinpath(@__DIR__,"..","mechanism","h2o2.yaml"))
    f=CounterflowDiffusionFlame(gas;fuel="H2:1",oxidizer="O2:1",
        T_fuel=800.,T_oxidizer=711.,P=ref["pressure"][1],
        mdot_fuel=ref["mdot_fuel"][1],mdot_oxidizer=ref["mdot_oxidizer"][1],grid=ref["grid"])
    n=gas.n_species
    function native_state(values)
        T=values[1,:]; Y=values[2:n+1,:]
        rho=f.pressure ./ (Arrhenius.R .* T .* vec(sum(Y./gas.MW;dims=1)))
        u=copy(values)
        u[1,:]./=1000
        u[end-2,:]./=Arrhenius._counterflow_Vscale
        u[end-1,:]./=Arrhenius._counterflow_Lscale
        u[end,:].*=rho
        return u
    end
    values=ref["values"]; attempt_dt=ref["attempt_dt"]
    factor=ref["factor"]; jac_evals=ref["jac_evals"]
    w=CounterflowWorkspace(f); r=similar(f.state)
    physical!(r,x,p,dt)=Arrhenius._counterflow_growth_residual!(r,f,x,w;previous=p,dt)
    # The generator reconstructs every accepted attempt size at full
    # precision, cross-checked against Cantera's printed sizes (4 significant
    # digits) with zero mismatches at fixture generation. A fixture that had
    # fallen back to the printed sizes would carry print rounding into every
    # applied factor; assert that invariant here instead of loosening the
    # gate for it.
    @test ref["reconstruction_ok"][1]
    # Fixed tolerances justified by measured reconstruction error, not chosen
    # from the fixture: with full-precision sizes the worst deviation of the
    # reconstructed growth factor from Cantera's applied factor over all
    # non-gated steps is 1.2e-7 relative (transport-property rounding at the
    # trajectory's ~3970 K states), so 1e-6 bounds it with an 8x margin. The
    # first attempt size is the sole print-seeded quantity (its full-precision
    # value is unknowable from the log), bounded by 5e-4 relative rounding of
    # the 4-significant-digit print, so its factor gets 2e-3.
    step_rtol(k)=k==1 ? 2e-3 : 1e-6
    grown=gated=flat=0
    for k in eachindex(attempt_dt)
        g=Arrhenius._timestep_growth(physical!,r,:residual_ratio,1.5,
            native_state(values[k,:,:]),native_state(values[k+1,:,:]),attempt_dt[k],0)
        if jac_evals[k]>0
            gated+=1
            @test factor[k]≈1 atol=1e-9
            @test g>=1.0
        elseif factor[k]>1+1e-6
            grown+=1
            @test g≈factor[k] rtol=step_rtol(k)
        else
            flat+=1
            @test g≈factor[k] rtol=step_rtol(k)
        end
    end
    # the trajectory must exercise grown, gated and flat decisions alike
    @test grown>50 && gated>50 && flat>50
end

@testset "counterflow Cantera bound fraction and trace species" begin
    # Mirror MultiNewton::boundStep: compute the largest admissible fraction,
    # landing exactly on the bound. This makes a subsequent outward Newton
    # direction a failed damped step eligible for Cantera's fresh-Jacobian retry.
    bound=Arrhenius._counterflow_bound_fraction
    lower,upper=-1e-7,1e5
    for y in (0.0,1e-25,-1e-25,-1e-7)
        alpha=bound(y,-1e-6,lower,upper)
        reference_trial=y-alpha*1e-6
        @test reference_trial>=lower-1e-30
        @test reference_trial≈max(y-1e-6,lower) atol=1e-30
    end
    @test bound(0.0,-1e-6,lower,upper)≈0.1 rtol=1e-14
    @test bound(1e-25,-1e-6,lower,upper)≈0.1 rtol=1e-14
    @test bound(-1e-25,-1e-6,lower,upper)≈0.1 rtol=1e-14
    @test bound(-1e-7,-1e-6,lower,upper)==0.0
    # A first step may land on the lower bound exactly; the next outward
    # Newton direction must return zero, matching Cantera's refresh trigger.
    @test bound(0.0,-1e-6,lower,upper)==0.1
    @test bound(lower,-1e-12,lower,upper)==0.0
    # Cantera's lower branch intentionally has no zero clamp: if the current
    # value is already below lower and the step moves farther down, the
    # returned fraction is negative. The upper branch does clamp this case.
    @test bound(lower-1e-7,-1e-6,lower,upper)≈-0.1 rtol=1e-14
    @test bound(upper+1e-7,1e-6,lower,upper)==0.0
end

@testset "counterflow transient Jacobian differential mask" begin
    gas=CreateSolution(joinpath(@__DIR__,"..","mechanism","h2o2.yaml"))
    f=CounterflowDiffusionFlame(gas;fuel="H2:1,AR:1",oxidizer="O2:.2,AR:.8",
        mdot_fuel=.24,mdot_oxidizer=.72,grid=collect(range(0.,.02;length=7)))
    w=CounterflowWorkspace(f); r=similar(f.state); previous=copy(f.state)
    counterflow_residual!(r,f,f.state,w;previous,dt=1e-5)
    old=copy(Arrhenius._counterflow_jacobian!(f,f.state,w,r;previous,dt=1e-5,analytic=false))
    updated=copy(old)
    Arrhenius._counterflow_transient_update!(updated,f,1e-5,Inf)
    B,N=size(f.state); diagonal=4B-1
    # Every interior species equation is differential; closure rows exist only
    # at boundaries and therefore receive no pseudo-time shift.
    for j in 1:N,k in 1:B
        shift=1<j<N && k<=gas.n_species+2 ? 10. : 0.
        @test updated[diagonal,(j-1)*B+k]-old[diagonal,(j-1)*B+k]≈shift atol=1e-10
    end
    counterflow_residual!(r,f,f.state,w)
    fresh=copy(Arrhenius._counterflow_jacobian!(f,f.state,w,r;analytic=false))
    @test maximum(abs.(updated.-fresh))<1e-4
    for dt in (1e-4,1e-6,1e-5)
        Arrhenius._counterflow_transient_update!(updated,f,Inf,dt)
        Arrhenius._counterflow_transient_update!(updated,f,dt,Inf)
    end
    @test maximum(abs.(updated.-fresh))<1e-4
end
# A deterministic Newton oracle isolates time-step policy from nonlinear solves.
mutable struct CounterflowBudgetFixture
    grid::Vector{Float64}
    state::Matrix{Float64}
    transient_calls::Int
    steady_calls::Int
end
mutable struct CounterflowBudgetWorkspace
    jac_evals::Int
    source_gradient::Matrix{Float64}
end
Arrhenius.CounterflowWorkspace(::CounterflowBudgetFixture)=CounterflowBudgetWorkspace(0,zeros(0,1))
function Arrhenius._counterflow_newton!(f::CounterflowBudgetFixture,w;previous=nothing,kwargs...)
    if previous===nothing
        f.steady_calls+=1
        return false
    end
    f.transient_calls+=1
    w.jac_evals+=1 # no growth in this oracle
    f.transient_calls>2 || return false
    f.state[1,1]+=1
    return true
end
mutable struct BoundRefreshWorkspace
    band::Matrix{Float64}
    jacobian::Matrix{Float64}
    pivots::Vector{LinearAlgebra.BlasInt}
    jac_age::Int
    jac_dt::Float64
    jac_valid::Bool
    jac_evals::Int
end
mutable struct BoundRefreshNewtonFixture
    state::Matrix{Float64}
    gas::NamedTuple
    control_points::Nothing
    dependent_species::Int
    residual_calls::Int
end
function Arrhenius.counterflow_residual!(r,f::BoundRefreshNewtonFixture,u,w;kwargs...)
    f.residual_calls+=1
    fill!(r,0); r[2,1]=1.0
    return r
end
function Arrhenius._counterflow_jacobian!(f::BoundRefreshNewtonFixture,u,w,r;kwargs...)
    fill!(w.band,0)
    for k in 1:length(u)
        w.band[15,k]=1
    end
    return w.band
end
@testset "counterflow stale Jacobian retries a negative bound failure" begin
    # The identical production failure retries only when the cached Jacobian
    # is stale (age > 1): the second attempt evaluates one fresh Jacobian,
    # then the fresh-Jacobian failure ends the solve without another retry.
    lower=-1e-7
    state=ones(4,1); state[2]=lower-1e-7
    f=BoundRefreshNewtonFixture(state,(n_species=1,),nothing,1,0)
    band=zeros(22,4)
    for k in 1:4
        band[15,k]=1
    end
    _,pivots=LinearAlgebra.LAPACK.gbtrf!(7,7,4,band)
    w=BoundRefreshWorkspace(copy(band),copy(band),pivots,2,Inf,true,0)
    @test !Arrhenius._counterflow_newton!(f,w;maxiters=4)
    @test w.jac_evals==1
    @test w.jac_age==1 # retry ended on the fresh Jacobian
    @test f.residual_calls==2 # initial residual on each damped attempt
end

@testset "counterflow negative bound fraction fails before trial evaluation" begin
    # Exercise the production Newton path, not just boundStep: the native
    # species value is already below its lower bound and Newton points farther
    # out. Cantera's asymmetric boundStep returns a negative fraction, which
    # dampStep rejects (<1e-10) before evaluating any trial residual.
    lower=-1e-7
    state=ones(4,1); state[2]=lower-1e-7
    f=BoundRefreshNewtonFixture(state,(n_species=1,),nothing,1,0)
    band=zeros(22,4)
    for k in 1:4
        band[15,k]=1
    end
    _,pivots=LinearAlgebra.LAPACK.gbtrf!(7,7,4,band)
    w=BoundRefreshWorkspace(copy(band),copy(band),pivots,0,Inf,false,0)
    @test Arrhenius._counterflow_bound_fraction(state[2,1],-1.0,lower,1.00001)<0
    @test !Arrhenius._counterflow_newton!(f,w;maxiters=1)
    @test f.residual_calls==1 # initial residual only; no damped trial was evaluated
    @test w.jac_evals==1
end

@testset "counterflow time-step budget and steady cadence" begin
    # The transient budget caps successful steps (the first two attempts fail
    # in this oracle), and an accepted final step returns before another steady
    # retry, matching Cantera's timeStep(nsteps) check in SteadyStateSystem.cpp.
    # Steady retries use the Cantera cadence: every ten successful steps.
    for (budget,steadies,accepted) in ((10,1,10),(21,3,21))
        f=CounterflowBudgetFixture([0.],zeros(1,1),0,0)
        dts=Float64[]
        steps0=Arrhenius.COUNTERFLOW_STATS.time_steps
        failures0=Arrhenius.COUNTERFLOW_STATS.time_step_failures
        @test !Arrhenius._counterflow_steady_impl!(f;max_time_steps=budget,step_callback=dt->push!(dts,dt))
        @test length(dts)==accepted
        @test Arrhenius.COUNTERFLOW_STATS.time_steps-steps0==budget
        @test Arrhenius.COUNTERFLOW_STATS.time_step_failures-failures0==2
        @test f.transient_calls==budget+2
        @test f.steady_calls==steadies
    end
end


# A constant algebraic residual can have a small weighted correction without
# satisfying the mass-flux constraint. The fallback must reject that state.
struct AlgebraicNewtonFixture
    state::Matrix{Float64}
    gas::NamedTuple
    control_points::Nothing
    dependent_species::Int
end
mutable struct AlgebraicNewtonWorkspace
    band::Matrix{Float64}
    jacobian::Matrix{Float64}
    pivots::Vector{LinearAlgebra.BlasInt}
    jac_age::Int
    jac_dt::Float64
    jac_valid::Bool
    jac_evals::Int
end
function Arrhenius.counterflow_residual!(r,f::AlgebraicNewtonFixture,u,w;kwargs...)
    fill!(r,0); r[end,1]=1e-6
    return r
end
function Arrhenius._counterflow_jacobian!(f::AlgebraicNewtonFixture,u,w,r;kwargs...)
    fill!(w.band,0)
    for k in 1:4
        w.band[15,k]=1
    end
    return w.band
end
@testset "weighted counterflow fallback preserves algebraic constraints" begin
    state=ones(4,1); state[1]=3
    f=AlgebraicNewtonFixture(state,(n_species=0,),nothing,1)
    w=AlgebraicNewtonWorkspace(zeros(22,4),zeros(22,4),LinearAlgebra.BlasInt[],21,Inf,false,0)
    # One accepted damped step has correction norm <1; its algebraic mass
    # residual remains1e-6 >1e-8, so maxiters fallback must return failure.
    @test !Arrhenius._counterflow_newton!(f,w;maxiters=1)
    @test f.state[end,1]<1 # verify that an accepted step preceded the fallback
end
