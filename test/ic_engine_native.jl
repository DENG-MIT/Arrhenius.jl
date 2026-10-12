module NativeEngineExampleTests
using Arrhenius
using Test
include(joinpath(@__DIR__,"..","example","reactors","ic_engine_setup.jl"))

@testset "engine geometry and valve/injection schedule" begin
    @test engine_volume(0.0) ≈ ENGINE_CLEARANCE
    @test engine_volume(0.01)/engine_volume(0.0) ≈ 20
    @test engine_volume(0.02) ≈ ENGINE_CLEARANCE
    @test engine_volume(0.03)/engine_volume(0.02) ≈ 20
    for t in (0.003,0.011,0.027,0.057,0.097,0.139)
        h = 1e-8
        derivative = (engine_volume(t+h)-engine_volume(t-h))/(2h)
        @test derivative ≈ -ENGINE_PISTON_AREA*engine_piston_speed(t) rtol=1e-8
    end
    @test engine_gate(0.0,-18,198)
    @test engine_gate(0.0,522,18)
    @test !engine_gate(0.0,350,365)
    @test engine_gate(0.020,350,365)
    @test !engine_gate(0.020,-18,198)
    @test !engine_gate(0.020,522,18)
    stops = engine_switching_times(0.16)
    @test first(stops) == 0.0
    @test last(stops) == 0.16
    @test all(diff(stops).>0)
    fuel_mass = sum(ENGINE_INJECTION_RATE*(stops[i+1]-stops[i])*
        engine_gate((stops[i]+stops[i+1])/2,350,365) for i in 1:length(stops)-1)
    @test fuel_mass ≈ 4*3.2e-5 rtol=1e-14
    for t in (0.003,0.014,0.0196,0.032)
        @test engine_gate(t,-18,198) == engine_gate(t+0.04,-18,198)
        @test engine_gate(t,350,365) == engine_gate(t+0.04,350,365)
        @test engine_gate(t,522,18) == engine_gate(t+0.04,522,18)
    end
    raw_times = collect(0.0:0.00001:0.0001)
    temperatures = [300.0,305,311,319,325,329,324,316,308,303,300]
    chosen = engine_output_indices(raw_times,temperatures)
    @test maximum(diff(raw_times[chosen])) <= 1/(360*ENGINE_FREQUENCY)
    @test maximum(abs,diff(temperatures[chosen])) <= 20
    @test maximum(temperatures[chosen]) == maximum(temperatures)
    @test first(chosen)==1 && last(chosen)==length(raw_times)
end
@testset "engine integral sampling" begin
    gas = (species_names=["co"],)
    output = Dict("time"=>collect(0.0:4.0),"heat_release_rate"=>collect(2.0:2.0:10.0),
        "work_rate"=>ones(5),"mean_molecular_weight"=>fill(20.0,5),
        "mdot_out"=>fill(0.2,5),"X"=>fill(1e-6,5,1))
    full = engine_integral_values(engine_quadrature_terms(output,gas;omit_initial=false))
    coarse = engine_integral_values(engine_quadrature_terms(output,gas;indices=[1,3,5],omit_initial=false))
    @test full.heat_J ≈ 24.0
    @test full.work_J ≈ 4.0
    @test full.CO_ppm ≈ 1.0
    @test full == coarse
    @test ic_engine_integrals(output,gas).heat_J ≈ 21.0
    t = [0.0,0.1,0.8,1.4,4.0]
    @test engine_quadratic_integral(2 .* t.^2 .+ 3 .* t .+ 5,t) ≈ 2*4^3/3+3*4^2/2+5*4
    @test engine_quadratic_integral([1.0,3.0],[0.0,1.0]) ≈ 2.0
    @test engine_quadratic_integral([1.0,2.0,3.0,4.0],[0.0,1.0,2.0,3.0]) ≈ 7.5
end
end
