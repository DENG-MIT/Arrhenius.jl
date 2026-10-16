using Arrhenius, LinearAlgebra, Test

@testset "failed premixed factorization invalidates persistent Jacobian" begin
    gas = CreateSolution(joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml"))
    f = FreeFlame(gas; X=mole_fractions(gas, "H2:1.1,O2:1,AR:5"), discretization=:finite_difference)
    w = Arrhenius.FlameWorkspace(f)
    r = flame_residual!(similar(f.state), f, f.state, w)
    copyto!(w.steady_jacobian, Arrhenius._flame_jacobian(f, f.state, w, r))

    # Establish a valid cached factorization, then make a cached transient
    # refactorization singular inside Newton. The failed attempt must return
    # false, invalidate all cache fields, and leave recovery to the next solve.
    Arrhenius._flame_factorize!(f, w, 0.0)
    @test w.jacobian_valid[]
    old_pivots = copy(w.jacobian_pivots[])
    @test !isempty(old_pivots)
    @test w.jacobian_rdt[] == 0.0
    fill!(w.steady_jacobian, 0.0) # boundary rows remain singular at any rdt
    @test !Arrhenius._flame_newton!(f, w; previous=f.state, dt=1.0)
    @test !w.jacobian_valid[]
    @test isempty(w.jacobian_pivots[])
    @test isnan(w.jacobian_rdt[])

    # The next Newton call must recompute the steady Jacobian and factor it.
    evaluations = w.jacobian_evals[]
    Arrhenius._flame_newton!(f, w; maxiters=1)
    @test w.jacobian_evals[] == evaluations + 1
    @test w.jacobian_valid[]
end

@testset "premixed factorization error classification" begin
    classify = Arrhenius._flame_recoverable_factorization_error
    @test classify(SingularException(1))
    @test classify(LAPACKException(1))
    @test !classify(LAPACKException(0))
    @test !classify(LAPACKException(-1))
    @test !classify(ErrorException("unrelated"))
end
