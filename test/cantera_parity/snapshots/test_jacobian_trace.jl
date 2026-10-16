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
