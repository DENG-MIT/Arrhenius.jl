using Test
include("jacobian_trace.jl")

@testset "Cantera per-attempt Jacobian attribution" begin
    fixture = joinpath(@__DIR__, "fixtures", "nakamura_diffusion_jacobian_trace.log")
    attempts = cantera_attempt_trace(fixture)
    per_attempt = getindex.(attempts, "jacobians")
    cumulative = getindex.(attempts, "jacobians_cumulative")
    entry_age = getindex.(attempts, "jacobian_age_start")

    # The recorded trace evaluates once on each retry, while Cantera's displayed
    # N_jac column is cumulative. These assertions fail if N_jac is attributed
    # directly to each attempt, or if the first stale-retry refresh is omitted.
    @test per_attempt == [1, 1, 1]
    @test cumulative == [1, 2, 3]
    @test entry_age == [1, 1, 1]
    @test cumulative != per_attempt
end

@testset "Cantera attempt closes before the next ten-step batch" begin
    fixture = joinpath(@__DIR__, "fixtures", "steady_batch_jacobian_trace.log")
    attempts = cantera_attempt_trace(fixture)
    @test getindex.(attempts, "dt") == [1e-5, 1e-5, 1.5e-5]
    @test getindex.(attempts, "jacobians") == [2, 0, 1]
    @test getindex.(attempts, "jacobians_cumulative") == [2, 2, 3]
    @test getindex.(attempts, "accepted") == [true, true, false]
end

@testset "Cantera solver stats distinguish cumulative Jacobians from steps" begin
    stats = Dict("jacobian_evals" => [3, 4, 4], "steps" => [1, 1, 1])
    @test cantera_stat_total(stats, "jacobian_evals") == 4
    @test cantera_stat_total(stats, "steps") == 3
end
