using Test
using LinearAlgebra
using SparseArrays
using ForwardDiff
using Arrhenius

# Load through the public optional entry point, not a test-only module path.
include(joinpath(@__DIR__, "..", "example", "reactors", "ode_solver.jl"))
@assert @isdefined(StructuredReactorSolver)
const SRS = StructuredReactorSolver
const KS = SRS.KLUCoreStructuredSolver
include(joinpath(@__DIR__, "..", "validation", "numerical_threads.jl"))
const STRUCTURED_THREAD_SETTINGS = benchmark_julia_thread_settings()
const H2O2 = joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml")

bitwise_equal(a::Array{Float64}, b::Array{Float64}) =
    reinterpret(UInt8, vec(a)) == reinterpret(UInt8, vec(b))

struct StructuredTestTag end
function signed_ad_column(reactor, state, column, time=0.0)
    D = ForwardDiff.Dual{StructuredTestTag,Float64,1}
    rhs = SRS.signed_trial_rhs(reactor; scalar_type=D)
    dual_state = [D(state[k], ForwardDiff.Partials((k == column ? 1.0 : 0.0,))) for k in eachindex(state)]
    dual_derivative = similar(dual_state)
    rhs(dual_derivative, dual_state, nothing, time)
    [ForwardDiff.partials(value)[1] for value in dual_derivative]
end

function one_reaction_gas(base, reactant, product; orders=copy(reactant),
        reversible=false, arrhenius=reshape([1.0e6, 0.0, 0.0], 1, 3),
        third=Int64[], falloff=Int64[], falloff_troe=Int64[],
        efficiencies=spzeros(Float64, base.n_species, 1),
        plog=Arrhenius.PlogData(Int64[], Int64[1], Float64[], Int64[1], zeros(0,3)),
        blowers_masel=Arrhenius.BlowersMaselData(Int64[], zeros(0,4)))
    ns = base.n_species
    vk = sparse(product - reactant)
    reaction = Arrhenius.Reaction(sparse(product), sparse(reactant), sparse(orders),
        [reversible], Matrix{Float64}(arrhenius), zeros(0,3), zeros(0,4),
        Vector{Int64}(third), Vector{Int64}(falloff), Vector{Int64}(falloff_troe),
        sparse(efficiencies), [Vector{Int64}(findnz(reactant[:,1])[1])],
        [Vector{Int64}(findnz(product[:,1])[1])], 1, vk,
        Vector{Float64}(vec(sum(vk,dims=1))), plog, blowers_masel)
    Arrhenius.Solution(ns, 1, copy(base.MW), copy(base.species_names),
        copy(base.elements), copy(base.ele_matrix), base.thermo, base.trans, reaction)
end

function anyn_gas(base, mode::Symbol)
    index(name) = something(findfirst(==(name), base.species_names))
    four = [index("H2"), index("O2"), index("H"), index("OH")]
    compact = [(index("H2O"), 2.0), (index("O"), 1.0)]
    reactant = spzeros(Float64, base.n_species, 1)
    product = spzeros(Float64, base.n_species, 1)
    if mode === :forward
        reactant[four,1] .= 1.0
        for (k,v) in compact; product[k,1] = v; end
    elseif mode in (:reversible_reverse, :irreversible_reverse)
        for (k,v) in compact; reactant[k,1] = v; end
        product[four,1] .= 1.0
    else
        throw(ArgumentError("unknown AnyN fixture mode"))
    end
    one_reaction_gas(base, reactant, product;
        reversible=mode === :reversible_reverse), four
end

function compare_guard!(guard, state; fallback)
    expected = zeros(length(state), length(state)); actual = similar(expected)
    saved = copy(state)
    if fallback
        SRS.signed_trial_jacobian!(expected, state, guard.analytic.float_rhs, 0.0)
    else
        guard.analytic(expected, state, nothing, 0.0)
    end
    guard(actual, state, nothing, 0.0)
    @test bitwise_equal(actual, expected)
    @test isequal(state, saved)
end

function klu_fixture(B, q, u, t, e, a)
    ns = size(B,1); core = sparse(Float64.(B))
    J = zeros(ns+1,ns+1)
    J[1:ns,1:ns] .= B .+ u*q'; J[1:ns,end] .= t
    J[end,1:ns] .= e; J[end,end] = a
    Dict{String,Any}("J"=>J, "core_colptr"=>copy(core.colptr),
        "core_rowval"=>copy(core.rowval), "core_nzval"=>copy(core.nzval),
        "q"=>Float64.(q), "rank_vector"=>Float64.(u),
        "temperature_column"=>vcat(Float64.(t),Float64(a)),
        "energy_row"=>Float64.(e))
end

@testset "signed product values and gradients" begin
    rows = [1,1,2,1,2,3,1,2,3,4]
    cols = [1,2,2,3,3,3,4,4,4,4]
    stoich = sparse(rows, cols, ones(length(rows)), 4, 4)
    plans = SRS.signed_product_plans(stoich, stoich)
    C = [2.0,3.0,5.0,7.0]; rate = 4.0
    for plan in plans
        gradient = zeros(length(plan.indices))
        SRS.signed_product_gradient!(gradient, plan, C, rate)
        oracle = ForwardDiff.gradient(x -> SRS.signed_multiply(plan,x,rate), C)
        @test gradient == oracle[plan.indices]
    end
    branch_cases = ((plans[1],[0.0,3.0,5.0,7.0],[4.0]),
        (plans[1],[-2.0,3.0,5.0,7.0],[4.0]),
        (plans[2],[0.0,3.0,5.0,7.0],[12.0,0.0]),
        (plans[2],[-2.0,3.0,5.0,7.0],[12.0,-8.0]),
        (plans[3],[0.0,3.0,5.0,7.0],[60.0,0.0,0.0]),
        (plans[3],[-2.0,3.0,5.0,7.0],[60.0,-40.0,-24.0]))
    for (plan, concentration, expected) in branch_cases
        actual = zeros(length(plan.indices))
        SRS.signed_product_gradient!(actual, plan, concentration, rate)
        @test actual == expected
    end
    for plan in plans[2:3]
        g = fill(NaN,length(plan.indices))
        SRS.signed_product_gradient!(g,plan,[-2.0,-3.0,5.0,7.0],rate)
        @test all(iszero,g)
        @test iszero(SRS.signed_multiply(plan,[-2.0,-3.0,5.0,7.0],rate))
    end
    g = fill(NaN,4)
    for C0 in ([0.0,3.0,5.0,7.0],[-2.0,3.0,5.0,7.0])
        SRS.signed_product_gradient!(g,plans[4],C0,rate)
        @test all(iszero,g)
    end
    invalid = SRS.SignedProductPlan([99],[1.0],[99],1)
    @test_throws BoundsError SRS.signed_multiply(invalid,C,rate)
    @test_throws DimensionMismatch SRS.signed_product_gradient!(zeros(1),plans[2],C,rate)
end

@testset "signed RHS and structured columns" begin
    gas = CreateSolution(H2O2)
    Y = collect(range(1.0,2.0;length=gas.n_species)); Y ./= sum(Y)
    reactor = IdealGasReactor(gas;temperature=1400.0,pressure=3one_atm,
        mass_fractions=Y,constraint=:constant_pressure,energy=:adiabatic)
    state = reactor_state(reactor); saved = copy(state)
    public_value = similar(state); reactor_rhs(reactor)(public_value,state,nothing,0.0)
    signed = SRS.signed_trial_rhs(reactor)
    signed_value = similar(state); signed(signed_value,state,nothing,0.0)
    # Product association differs from the public clipped-rate loop.
    @test all(isapprox.(signed_value,public_value;rtol=32eps(Float64),atol=0))
    negative = copy(state)
    negative[something(findfirst(==("H"),gas.species_names))] = -1e-7
    negative_saved = copy(negative)
    clipped = similar(state); reactor_rhs(reactor)(clipped,negative,nothing,0.0)
    signed_negative = similar(state); signed(signed_negative,negative,nothing,0.0)
    @test clipped != signed_negative
    @test all(isfinite,clipped) && all(isfinite,signed_negative)
    @test isequal(negative,negative_saved)

    jacobian = SRS.structured_trial_jacobian(reactor)
    J = zeros(length(state),length(state)); jacobian(J,state,nothing,0.0)
    for column in eachindex(state)
        oracle = signed_ad_column(reactor,state,column)
        @test norm(@view(J[:,column])-oracle,Inf) <= 1e-9*max(norm(oracle,Inf),1.0)
    end
    @test state == saved == reactor_state(reactor)

    # Purpose-built fixed-pressure PLOG coverage.
    h2 = something(findfirst(==("H2"),gas.species_names))
    h = something(findfirst(==("H"),gas.species_names))
    reactant = spzeros(Float64,gas.n_species,1); reactant[h2,1] = 1
    product = spzeros(Float64,gas.n_species,1); product[h,1] = 2
    plog = Arrhenius.PlogData(Int64[1],Int64[0],Int64[1,3],
        Float64[0.5one_atm,2one_atm],Int64[1,2,3],Float64[1e6 0 0;2e6 0 0])
    plog_gas = one_reaction_gas(gas,reactant,product;plog)
    plog_reactor = IdealGasReactor(plog_gas;temperature=1200.0,pressure=one_atm,
        mass_fractions=Y,constraint=:constant_pressure,energy=:adiabatic)
    plog_state = reactor_state(plog_reactor)
    plog_jacobian = SRS.structured_trial_jacobian(plog_reactor)
    plog_J = zeros(length(plog_state),length(plog_state))
    plog_jacobian(plog_J,plog_state,nothing,0.0)
    for column in eachindex(plog_state)
        oracle = signed_ad_column(plog_reactor,plog_state,column)
        @test norm(@view(plog_J[:,column])-oracle,Inf) <= 1e-9*max(norm(oracle,Inf),1.0)
    end
end

@testset "AnyN zero guard and mutation propagation" begin
    base = CreateSolution(H2O2)
    Y = collect(range(1.0,2.0;length=base.n_species)); Y ./= sum(Y)
    gas,_ = anyn_gas(base,:forward)
    reactor = IdealGasReactor(gas;temperature=1000.0,pressure=one_atm,
        mass_fractions=Y,constraint=:constant_pressure,energy=:adiabatic)
    guard = SRS.guarded_structured_jacobian(reactor); state = reactor_state(reactor)
    compare_guard!(guard,state;fallback=false)
    plan = only(guard.analytic.float_rhs.forward_plans)
    probe_C = [0.0,2.0,1.5,0.8]
    synthetic = SRS.SignedProductPlan(collect(1:4),ones(4),Int[],4)
    @test SRS._ambiguous_anyn_zero(synthetic,probe_C)
    @test !SRS._ambiguous_anyn_zero(synthetic,[0.0,0.0,1.5,0.8])
    @test !SRS._ambiguous_anyn_zero(
        SRS.SignedProductPlan(collect(1:4),[2.0,1.0,1.0,1.0],Int[],4),probe_C)
    @test !SRS._ambiguous_anyn_zero(synthetic,[-1e-3,2.0,1.5,0.8])
    zero_species = plan.indices[something(findfirst(==(1.0),plan.orders))]
    for value in (0.0,nextfloat(0.0))
        trial = copy(state); trial[zero_species] = value
        guard.analytic.float_rhs(guard.analytic.float_rhs.jac_base,trial,nothing,0.0)
        C = guard.analytic.float_rhs.workspace.C
        @test iszero(C[zero_species])
        @test SRS._ambiguous_anyn_zero(plan,C)
        compare_guard!(guard,trial;fallback=true)
    end
    @test (guard.analytic_calls,guard.fallback_calls,guard.generation) == (1,2,3)

    reverse_gas,reverse_species = anyn_gas(base,:reversible_reverse)
    rr = IdealGasReactor(reverse_gas;temperature=1000.0,pressure=one_atm,
        mass_fractions=Y,constraint=:constant_pressure,energy=:adiabatic)
    rg = SRS.guarded_structured_jacobian(rr); rs = reactor_state(rr)
    rs[first(reverse_species)] = 0
    compare_guard!(rg,rs;fallback=true)
    @test rg.fallback_calls == 1
    irreversible_gas,irreversible_species = anyn_gas(base,:irreversible_reverse)
    ir = IdealGasReactor(irreversible_gas;temperature=1000.0,pressure=one_atm,
        mass_fractions=Y,constraint=:constant_pressure,energy=:adiabatic)
    ig = SRS.guarded_structured_jacobian(ir); is = reactor_state(ir)
    is[first(irreversible_species)] = 0
    compare_guard!(ig,is;fallback=false)
    @test ig.fallback_calls == 0

    before = (guard.analytic_calls,guard.fallback_calls,guard.generation)
    @test_throws DimensionMismatch guard(zeros(length(state),length(state)),state[1:end-1],nothing,0.0)
    nonfinite = copy(state); nonfinite[end] = NaN
    @test_throws DomainError guard(zeros(length(state),length(state)),nonfinite,nothing,0.0)
    original_mw = gas.MW[1]
    try
        gas.MW[1] = nextfloat(original_mw)
        @test_throws ArgumentError guard(zeros(length(state),length(state)),state,nothing,0.0)
    finally
        gas.MW[1] = original_mw
    end
    @test (guard.analytic_calls,guard.fallback_calls,guard.generation) == before
end

@testset "scope, corrupt data, and wrapper fallback" begin
    gas = CreateSolution(H2O2); Y = fill(1/gas.n_species,gas.n_species)
    cp = IdealGasReactor(gas;temperature=1200.0,pressure=one_atm,mass_fractions=Y,
        constraint=:constant_pressure,energy=:adiabatic)
    cv = IdealGasReactor(gas;temperature=1200.0,pressure=one_atm,mass_fractions=Y,
        constraint=:constant_volume,energy=:adiabatic)
    iso = IdealGasReactor(gas;temperature=1200.0,pressure=one_atm,mass_fractions=Y,
        constraint=:constant_pressure,energy=:isothermal)
    @test SRS.try_structured_adapter(reactor_problem(cp,(0.0,1e-8))) !== nothing
    @test SRS.try_guarded_jacobian(reactor_problem(cp,(0.0,1e-8))) !== nothing
    @test SRS.try_structured_adapter(reactor_problem(cv,(0.0,1e-8))) === nothing
    @test SRS.try_guarded_jacobian(reactor_problem(cv,(0.0,1e-8))) === nothing
    @test SRS.try_structured_adapter(reactor_problem(iso,(0.0,1e-8))) === nothing
    @test SRS.try_structured_adapter(merge(reactor_problem(cp,(0.0,1e-8)),
        (u0=Float32.(reactor_state(cp)),))) === nothing
    @test SRS.try_structured_adapter(merge(reactor_problem(cp,(0.0,1e-8)),
        (f=(du,u,p,t)->fill!(du,0),))) === nothing

    @test_throws SRS.UnsupportedStructuredReactor SRS._modal_efficiency([0.0;1.0;;],1)
    @test_throws SRS.UnsupportedStructuredReactor SRS._modal_efficiency([2.0;2.0;1.0;;],1)
    @test_throws ArgumentError SRS._modal_efficiency([-1.0;0.0;0.0;;],1)
    @test_throws DomainError SRS._modal_efficiency([NaN;0.0;0.0;;],1)

    reactant = spzeros(Float64,gas.n_species,1); reactant[1,1] = 1
    product = spzeros(Float64,gas.n_species,1); product[2,1] = 1
    custom = copy(reactant); custom[1,1] = 2
    custom_gas = one_reaction_gas(gas,reactant,product;orders=custom)
    @test_throws SRS.UnsupportedStructuredReactor SRS.structured_trial_jacobian(
        IdealGasReactor(custom_gas;temperature=1200.0,pressure=one_atm,mass_fractions=Y))
    fractional = copy(reactant); fractional[1,1] = 0.5
    fractional_gas = one_reaction_gas(gas,fractional,product;orders=copy(fractional))
    @test_throws SRS.UnsupportedStructuredReactor SRS.structured_trial_jacobian(
        IdealGasReactor(fractional_gas;temperature=1200.0,pressure=one_atm,mass_fractions=Y))
    bm = Arrhenius.BlowersMaselData(Int64[1],reshape([1e6,0.0,0.0,0.0],1,4))
    bm_gas = one_reaction_gas(gas,reactant,product;blowers_masel=bm)
    @test_throws SRS.UnsupportedStructuredReactor SRS.structured_trial_jacobian(
        IdealGasReactor(bm_gas;temperature=1200.0,pressure=one_atm,mass_fractions=Y))
    collider_plog = Arrhenius.PlogData(Int64[1],Int64[1],Int64[1,2],
        Float64[one_atm],Int64[1,2],Float64[1e6 0 0])
    plog_gas = one_reaction_gas(gas,reactant,product;plog=collider_plog)
    @test_throws SRS.UnsupportedStructuredReactor SRS.structured_trial_jacobian(
        IdealGasReactor(plog_gas;temperature=1200.0,pressure=one_atm,mass_fractions=Y))
    efficiencies = sparse(fill(1.0,gas.n_species,1))
    overlap_gas = one_reaction_gas(gas,reactant,product;third=Int64[1],falloff=Int64[1],
        falloff_troe=Int64[-1],efficiencies)
    @test_throws SRS.UnsupportedStructuredReactor SRS.structured_trial_jacobian(
        IdealGasReactor(overlap_gas;temperature=1200.0,pressure=one_atm,mass_fractions=Y))

    invalid = copy(reactant); invalid[1,1] = -1
    invalid_gas = one_reaction_gas(gas,invalid,product;orders=copy(invalid))
    @test_throws ArgumentError SRS.structured_trial_jacobian(
        IdealGasReactor(invalid_gas;temperature=1200.0,pressure=one_atm,mass_fractions=Y))
    nonfinite = copy(reactant); nonfinite[1,1] = NaN
    nonfinite_gas = one_reaction_gas(gas,nonfinite,product;orders=copy(nonfinite))
    @test_throws DomainError SRS.structured_trial_jacobian(
        IdealGasReactor(nonfinite_gas;temperature=1200.0,pressure=one_atm,mass_fractions=Y))
    out_of_bounds = one_reaction_gas(gas,reactant,product;third=Int64[2])
    @test_throws BoundsError SRS.structured_trial_jacobian(
        IdealGasReactor(out_of_bounds;temperature=1200.0,pressure=one_atm,mass_fractions=Y))
    malformed = one_reaction_gas(gas,reactant,product;
        orders=spzeros(Float64,gas.n_species,0))
    malformed_reactor = IdealGasReactor(
        malformed;temperature=1200.0,pressure=one_atm,mass_fractions=Y)
    @test_throws DimensionMismatch SRS.structured_trial_jacobian(malformed_reactor)
    malformed_problem = reactor_problem(malformed_reactor,(0.0,1e-8))
    @test_throws DimensionMismatch SRS.try_guarded_jacobian(malformed_problem)
    @test_throws DimensionMismatch SRS.try_structured_adapter(malformed_problem)
end

@testset "KLU fixed-pattern and singular-core recovery" begin
    data = klu_fixture([4.0 1.0;2.0 3.0],[0.1,0.2],[0.3,-0.1],
        [0.2,0.4],[-0.2,0.1],2.5)
    solver = KS.new_structured_solver(data); x = zeros(3)
    for (sigma,b) in ((1.0,[1.0,-2.0,0.5]),(1.5,[-0.5,1.2,2.0]))
        KS.factor_shift!(solver,sigma); KS.solve_factored!(x,solver,b)
        @test x ≈ (data["J"]-sigma*I)\b rtol=2e-13 atol=2e-13
    end
    KS.solve_factored!(x,solver,[0.2,0.3,-0.1])
    report = KS.solver_report(solver)
    @test report["fresh_factors"] == 1
    @test report["successful_refactors"] == 1
    @test report["solves"] == 3
    saved_row = solver.shifted.rowval[1]; solver.shifted.rowval[1] = 2
    @test_throws ErrorException KS.factor_shift!(solver,2.0)
    solver.shifted.rowval[1] = saved_row

    singular = klu_fixture(reshape([1.0],1,1),[1.0],[1.0],[0.0],[0.0],2.0)
    recovery = KS.new_structured_solver(singular)
    recovery_x = zeros(2); recovery_b = [1.0,2.0]
    for sigma in (0.5,1.0,3.0)
        KS.factor_shift!(recovery,sigma)
        KS.solve_factored!(recovery_x,recovery,recovery_b)
        @test recovery_x ≈ (singular["J"]-sigma*I)\recovery_b rtol=2e-13 atol=2e-13
    end
    rr = KS.solver_report(recovery)
    @test rr["fresh_factors"] == 2
    @test rr["discarded_refactors"] == 1
    @test rr["fallback_shift_calls"] == 1
    @test rr["fallback_solves"] == 1
    @test rr["recovery_rebuilds"] == 1
end

@testset "native_bdf selection, lifecycle, and balances" begin
    gas = CreateSolution(H2O2); mixture = Dict("H2"=>2.0,"O2"=>1.0,"AR"=>4.0)
    cp = IdealGasReactor(gas;temperature=1400.0,pressure=3one_atm,
        mole_fractions=mixture,constraint=:constant_pressure,energy=:adiabatic)
    problem = reactor_problem(cp,(0.0,1e-8))
    initial_problem_state = copy(problem.u0)
    @test_throws ArgumentError native_bdf(problem;linear_solver=:invalid)
    public_solution = native_bdf(problem;reltol=1e-9,abstol=1e-15,
        saveat=[0.0,1e-8],maxiters=100_000)
    @test SciMLBase.successful_retcode(public_solution)
    @test public_solution.prob.f.jac isa SRS.QNDFStructuredAdapter

    adapter = SRS.try_structured_adapter(problem); integrator_ref = Ref{Any}(nothing)
    precs = SRS.adapter_precs(adapter,integrator_ref)
    linsolve = SRS.adapter_linsolve(adapter,integrator_ref)
    f = SciMLBase.ODEFunction(adapter.guard.analytic.float_rhs;jac=adapter,tgrad=problem.tgrad)
    ode = SciMLBase.ODEProblem(f,problem.u0,problem.tspan,problem.p)
    integ = SciMLBase.init(ode,OrdinaryDiffEqBDF.QNDF(precs=precs,linsolve=linsolve);
        reltol=1e-9,abstol=1e-15,saveat=[0.0,1e-8],maxiters=100_000,
        isoutofdomain=(u,p,t)->!all(isfinite,u)||u[end]<=0)
    integrator_ref[] = integ
    structured_solution = SciMLBase.solve!(integ)
    @test SciMLBase.successful_retcode(structured_solution)
    @test structured_solution.destats.nnonlinconvfail == 0
    report = SRS.adapter_report(adapter)
    @test SRS._lifecycle_valid(adapter) && report["lifecycle_passed"]
    @test report["real_hook_calls"] == report["solves"] == structured_solution.destats.nsolve
    @test report["new_w_calls"] == report["factorizations"] == structured_solution.destats.nw
    @test report["jacobian_generations"] == adapter.guard.generation
    @test public_solution.u == structured_solution.u

    dense_W_observed = Ref(false)
    dense_algorithm_type = Ref("")
    dense_algorithm_enum = Ref("")
    dense_qr_fallback = Ref(false)
    dense_callback = SciMLBase.DiscreteCallback(
        (u,t,integrator)->true,
        integrator->begin
            cache = integrator.cache.nlsolver.cache.linsolve
            dense_W_observed[] = cache.A isa Matrix{Float64} &&
                size(cache.A) == (length(problem.u0),length(problem.u0))
            dense_algorithm_type[] = string(typeof(cache.alg))
            hasproperty(cache.alg,:alg) &&
                (dense_algorithm_enum[] = string(getproperty(cache.alg,:alg)))
            hasproperty(cache.cacheval,:fell_back_to_qr) &&
                (dense_qr_fallback[] |= cache.cacheval.fell_back_to_qr)
            SciMLBase.u_modified!(integrator,false)
        end;save_positions=(false,false))
    dense_solution = native_bdf(problem;linear_solver=:dense,reltol=1e-9,abstol=1e-15,
        saveat=[0.0,1e-8],maxiters=100_000,callback=dense_callback)
    @test SciMLBase.successful_retcode(dense_solution)
    dense_guard = dense_solution.prob.f.jac
    @test dense_guard isa SRS.GuardedStructuredJacobian
    @test dense_guard.generation == dense_solution.destats.njacs > 0
    @test dense_guard.analytic_calls + dense_guard.fallback_calls == dense_guard.generation
    @test dense_solution.destats.nnonlinconvfail == 0
    @test dense_W_observed[]
    @test !occursin("KLU",dense_algorithm_type[])
    @test problem.u0 == initial_problem_state == reactor_state(cp)
    @info "native_bdf dense backend" algorithm_type=dense_algorithm_type[] algorithm_enum=dense_algorithm_enum[] qr_fallback_observed=dense_qr_fallback[] dense_W=dense_W_observed[] jacobians=dense_guard.generation

    initial = reactor_properties(cp,problem.u0)
    final = reactor_properties(cp,structured_solution.u[end])
    @test abs(final.mass_fraction_sum-1) < 5e-10
    initial_elements = gas.ele_matrix*(problem.u0[1:end-1]./gas.MW)
    final_elements = gas.ele_matrix*(structured_solution.u[end][1:end-1]./gas.MW)
    @test norm(final_elements-initial_elements,Inf) < 5e-11
    @test abs(final.enthalpy/initial.enthalpy-1) < 1e-6
    dense_final = reactor_properties(cp,dense_solution.u[end])
    @test abs(dense_final.mass_fraction_sum-1) < 5e-10
    dense_elements = gas.ele_matrix*(dense_solution.u[end][1:end-1]./gas.MW)
    @test norm(dense_elements-initial_elements,Inf) < 5e-11
    @test abs(dense_final.enthalpy/initial.enthalpy-1) < 1e-6
    @test minimum(dense_solution.u[end][1:end-1]) > -1e-12

    for (constraint,energy) in ((:constant_volume,:adiabatic),(:constant_pressure,:isothermal))
        reactor = IdealGasReactor(gas;temperature=1200.0,pressure=one_atm,
            mole_fractions=mixture,constraint,energy)
        fallback_problem = reactor_problem(reactor,(0.0,1e-8))
        @test SRS.try_structured_adapter(fallback_problem) === nothing
        solution = native_bdf(fallback_problem;reltol=1e-8,abstol=1e-14,
            saveat=[0.0,1e-8],maxiters=100_000)
        @test SciMLBase.successful_retcode(solution)
        @test all(isfinite,reduce(vcat,solution.u))
        dense_fallback = native_bdf(fallback_problem;linear_solver=:dense,
            reltol=1e-8,abstol=1e-14,saveat=[0.0,1e-8],maxiters=100_000)
        @test SciMLBase.successful_retcode(dense_fallback)
        @test all(isfinite,reduce(vcat,dense_fallback.u))
        @test !(dense_fallback.prob.f.jac isa SRS.GuardedStructuredJacobian)
    end
    @test benchmark_julia_thread_settings(enforce=false) == STRUCTURED_THREAD_SETTINGS
end
