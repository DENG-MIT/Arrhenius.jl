using Test, Arrhenius

@testset "linear-Burke effective pressure" begin
    # One PLOG reaction, two pressure groups (k = 1e10 at 1e4 Pa, 1e12 at 1e6 Pa),
    # and one listed collider (species 2) with temperature-independent efficiency 2.
    coefficients = [1e10 0.0 0.0; 1e12 0.0 0.0]
    plog = Arrhenius.PlogData([1], [0], [1, 3], [1e4, 1e6], [1, 2, 3], coefficients,
        [1, 2], [2], [2.0 0.0 0.0])
    plain = Arrhenius.PlogData([1], [0], [1, 3], [1e4, 1e6], [1, 2, 3], coefficients)
    T, P, C = 1000.0, 1e5, [1.0, 3.0, 0.0]
    # eps_mix = (1*1 + 3*2 + 0*1) / 4
    @test Arrhenius._plog_effective_pressure(plog, 1, T, P, C) ≈ 1.75P
    @test Arrhenius._plog_effective_pressure(plain, 1, T, P, C) == P
    @test !Arrhenius._plog_linear_burke(plain, 1) && Arrhenius._plog_linear_burke(plog, 1)
    @test Arrhenius._plog_log_slope(plog, 1, T, 1e5, log(T)) ≈ 1.0
    @test Arrhenius._plog_log_slope(plog, 1, T, 1e7, log(T)) == 0.0
    @test Arrhenius._plog_rate(plog, 1, T, 1.75P, log(T)) ≈ 1e11 * 1.75
end
