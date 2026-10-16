using Arrhenius, LinearAlgebra, Test

@testset "failed premixed factorization invalidates persistent Jacobian" begin
    gas = CreateSolution(joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml"))
    f = FreeFlame(gas; X=mole_fractions(gas, "H2:1.1,O2:1,AR:5"), discretization=:finite_difference)
    w = Arrhenius.FlameWorkspace(f)
    r = flame_residual!(similar(f.state), f, f.state, w)
    copyto!(w.steady_jacobian, Arrhenius._flame_jacobian(f, f.state, w, r))

    # Establish a valid cached factorization, then make the next transient
    # factorization singular while an old cache exists.
    Arrhenius._flame_factorize!(f, w, 0.0)
    w.jacobian_valid[] = true
    old_pivots = copy(w.jacobian_pivots[])
    @test !isempty(old_pivots)
    @test w.jacobian_rdt[] == 0.0
    fill!(w.steady_jacobian, 0.0) # boundary rows remain singular at any rdt
    @test_throws LAPACKException Arrhenius._flame_factorize!(f, w, 1.0)
    @test !w.jacobian_valid[]
    @test isempty(w.jacobian_pivots[])
    @test isnan(w.jacobian_rdt[])

    # A subsequent Newton call must evaluate a fresh Jacobian instead of
    # consuming the failed factors and stale pivots.
    copyto!(w.steady_jacobian, Arrhenius._flame_jacobian(f, f.state, w, r))
    evaluations = w.jacobian_evals[]
    Arrhenius._flame_newton!(f, w; maxiters=1)
    @test w.jacobian_evals[] == evaluations + 1
    @test w.jacobian_valid[]
end
