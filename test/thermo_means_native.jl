module NativeThermoMeanTests
using Arrhenius, Test, LinearAlgebra, ForwardDiff

const tm_gas = CreateSolution(joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml"))

# Fixture style mirrors test/species_thermo_native.jl; duplicated here so this
# file does not pull in that file's testsets.
function tm_fixture(models)
    return Dict{String,Any}("phases"=>[Dict("species"=>["species_$i" for i in eachindex(models)])],
        "species"=>[Dict("name"=>"species_$i","thermo"=>model) for (i,model) in enumerate(models)])
end

function tm_wrap_solution(thermo, MW)
    n = size(thermo.nasa_low, 1)
    return Arrhenius.Solution(n, 0, MW, ["species_$i" for i in 1:n],
        tm_gas.elements, zeros(length(tm_gas.elements), n), thermo,
        tm_gas.trans, tm_gas.reaction)
end

const tm_mixed_gas = tm_wrap_solution(IdealGasThermo(tm_fixture([
    Dict("model"=>"NASA7","temperature-ranges"=>[200.,1000.,6000.],
        "data"=>[[2.8,1.2e-3,-2.5e-7,4.0e-11,-3.0e-15,950.,3.5],
                 [3.4,-5.0e-4,1.0e-7,-1.0e-11,5.0e-16,1150.,5.2]]),
    Dict("model"=>"NASA9","temperature-ranges"=>[200.,1000.,6000.,20000.],
        "data"=>[[120.,-1.5e-3,2.6+0.2r,1.0e-4*r,-2.0e-8*r,4.0e-12,-3.0e-16,
                  8200.0+300r,4.0+0.5r] for r in 0:2]),
    Dict("model"=>"Shomate","temperature-ranges"=>[200.,1000.,6000.],
        "data"=>[[28.0+2r,4.5-0.5r,-1.2,0.06,0.004,850.0+40r,190.0+10r] for r in 0:1]),
    Dict("model"=>"constant-cp","T0"=>"1000 K","h0"=>"9.22 kcal/mol",
        "s0"=>"-3.02 cal/mol/K","cp0"=>"5.95 cal/(mol*K)"),
])), [2.016, 18.015, 28.014, 44.010])

# Normalized and nonnormalized compositions, each with an exact zero and a
# trace fraction above the 1e-30 entropy cutoff.
function tm_compositions(n)
    base = collect(range(0.05, 0.30, n))
    base[1] = 0.0
    base[2] = 1.0e-15
    normalized = base / sum(base)
    return normalized, 2.7 * normalized
end

const tm_pairs = ((cal_cp, cal_cp_mean, cal_cpmass_mean),
            (cal_cv, cal_cv_mean, cal_cvmass_mean),
            (cal_h, cal_h_mean, cal_hmass_mean),
            (cal_u, cal_u_mean, cal_umass_mean),
            (cal_s0, cal_s0_mean, cal_s0mass_mean),
            (cal_s, cal_s_mean, cal_smass_mean),
            (cal_g, cal_g_mean, cal_gmass_mean),
            (cal_a, cal_a_mean, cal_amass_mean))

@testset "mixture means match species-vector contraction" begin
    for (gas, Ts) in ((tm_gas, (800., 2500.)), (tm_mixed_gas, (800., 2500., 7000.)))
        Xn, Xnn = tm_compositions(gas.n_species)
        for T in Ts, P in (one_atm, 3.0e5), X in (Xn, Xnn)
            for (species_fn, mean_fn, mass_fn) in tm_pairs
                reference = dot(X, species_fn(gas, T, P, X))
                @test mean_fn(gas, T, P, X) ≈ reference
                @test mass_fn(gas, T, P, X) ≈ reference / dot(X, gas.MW)
            end
        end
    end
end

@testset "mixture thermodynamic identities" begin
    for (gas, T, P) in ((tm_gas, 1500., 2.0e5), (tm_mixed_gas, 3200., 0.7e5))
        _, X = tm_compositions(gas.n_species)
        sumX = sum(X)
        cp_mean = cal_cp_mean(gas, T, P, X)
        h_mean = cal_h_mean(gas, T, P, X)
        s_mean = cal_s_mean(gas, T, P, X)
        u_mean = cal_u_mean(gas, T, P, X)
        g_mean = cal_g_mean(gas, T, P, X)
        @test cal_cv_mean(gas, T, P, X) ≈ cp_mean - R * sumX
        @test u_mean ≈ h_mean - R * T * sumX
        @test g_mean ≈ h_mean - T * s_mean
        @test cal_a_mean(gas, T, P, X) ≈ u_mean - T * s_mean
        @test cal_a_mean(gas, T, P, X) ≈ g_mean - R * T * sumX
        @test cal_s0_mean(gas, T, P, X) - s_mean ≈
            R * (dot(X, log.(max.(X, 1.0e-30))) + sumX * log(P / one_atm))
    end
end

@testset "mixture mean derivatives" begin
    for (gas, T, P) in ((tm_gas, 1500., 2.0e5), (tm_mixed_gas, 3200., 0.7e5))
        X, _ = tm_compositions(gas.n_species)
        @test ForwardDiff.derivative(t -> cal_h_mean(gas, t, P, X), T) ≈
            cal_cp_mean(gas, T, P, X) rtol=1e-8
        @test ForwardDiff.derivative(t -> cal_u_mean(gas, t, P, X), T) ≈
            cal_cv_mean(gas, T, P, X) rtol=1e-8
        @test ForwardDiff.derivative(p -> cal_s_mean(gas, T, p, X), P) ≈
            -R * sum(X) / P rtol=1e-8
        @test ForwardDiff.derivative(p -> cal_g_mean(gas, T, p, X), P) ≈
            R * T * sum(X) / P rtol=1e-8
        @test ForwardDiff.derivative(p -> cal_cp_mean(gas, T, p, X), P) ≈ 0 atol=1e-10
        @test ForwardDiff.gradient(x -> cal_cp_mean(gas, T, P, x), X) ≈ cal_cp(gas, T, P, X)
        @test ForwardDiff.gradient(x -> cal_h_mean(gas, T, P, x), X) ≈ cal_h(gas, T, P, X)
        grad_s = ForwardDiff.gradient(x -> cal_s_mean(gas, T, P, x), X)
        s_species = cal_s(gas, T, P, X)
        # Above the cutoff the mixing term contributes an extra -R; at and
        # below it the log is pinned and the gradient is the species value.
        expected = [x > 1.0e-30 ? s - R : s for (x, s) in zip(X, s_species)]
        @test grad_s ≈ expected
    end
end

@testset "Float32 thermo promotion" begin
    thermo32 = IdealGasThermo{Float32}(tm_gas.thermo.nasa_low, tm_gas.thermo.nasa_high,
        tm_gas.thermo.Trange, tm_gas.thermo.isTcommon)
    gas32 = tm_wrap_solution(thermo32, tm_gas.MW)
    X, _ = tm_compositions(gas32.n_species)
    X32 = Float32.(X)
    T32, P32 = 1500f0, Float32(2.0e5)
    for (_, mean_fn, _) in tm_pairs
        value32 = mean_fn(gas32, T32, P32, X32)
        @test value32 isa Float32
        @test value32 ≈ mean_fn(gas32, Float64(T32), Float64(P32), Float64.(X32)) rtol=1e-5
    end
end

@testset "allocation-free mixture means" begin
    X, _ = tm_compositions(tm_gas.n_species)
    Xm, _ = tm_compositions(tm_mixed_gas.n_species)
    T, P = 1500., 2.0e5
    cal_cp_mean(tm_gas, T, P, X)
    cal_h_mean(tm_gas, T, P, X)
    cal_s_mean(tm_gas, T, P, X)
    cal_g_mean(tm_gas, T, P, X)
    cal_cpmass_mean(tm_gas, T, P, X)
    cal_s_mean(tm_mixed_gas, T, P, Xm)
    @test @allocated(cal_cp_mean(tm_gas, T, P, X)) == 0
    @test @allocated(cal_h_mean(tm_gas, T, P, X)) == 0
    @test @allocated(cal_s_mean(tm_gas, T, P, X)) == 0
    @test @allocated(cal_g_mean(tm_gas, T, P, X)) == 0
    @test @allocated(cal_cpmass_mean(tm_gas, T, P, X)) == 0
    @test @allocated(cal_s_mean(tm_mixed_gas, T, P, Xm)) == 0
    @test_throws DimensionMismatch cal_cp_mean(tm_gas,T,P,X[2:end])
    @test cal_cp_mean(tm_gas,T,0.,X) == cal_cp_mean(tm_gas,T,P,X)
    @test cal_h_mean(tm_gas,T,0.,X) == cal_h_mean(tm_gas,T,P,X)
    @test cal_s0_mean(tm_gas,T,0.,X) == cal_s0_mean(tm_gas,T,P,X)
end
end
