using Arrhenius, NPZ, LinearAlgebra, Test
BLAS.set_num_threads(1)
gas = CreateSolution(ARGS[1])
data = MultiTransportData(ARGS[1]*".multicomponent.npz",gas;mechanism=ARGS[1])
f = FreeFlame(gas;X="H2:1.1,O2:1,AR:5")
slope = length(ARGS)>2 ? parse(Float64,ARGS[3]) : .06
curve = length(ARGS)>3 ? parse(Float64,ARGS[4]) : .12
@testset "premixed transport models" begin
for (model,basis,soret) in ((:mixture_averaged,:mole,false),
        (:mixture_averaged,:mass,false),(:mixture_averaged,:mass,true),(:multicomponent,:mass,false),
        (:multicomponent,:mass,true))
    set_transport!(f,model;data,flux_gradient_basis=basis,soret)
    elapsed = @elapsed solve!(f; slope,curve,loglevel=1)
    label = replace(String(model),"_"=>"-")*"-"*String(basis)*"-"*string(Int(soret))
    reference = npzread(joinpath(ARGS[2],label*".npz"))
    npzwrite(joinpath(ARGS[2],"julia-"*label*".npz"),
        Dict("grid"=>f.grid,"T"=>temperature(f),"Y"=>mass_fractions(f),"speed"=>[flame_speed(f)]))
    println("case=",label," speed=",flame_speed(f)," Tmax=",maximum(temperature(f)),
        " points=",length(f.grid)," seconds=",elapsed)
    @test f.converged
    @test flame_speed(f) ≈ only(reference["speed"]) rtol=.01
    @test maximum(temperature(f)) ≈ maximum(reference["T"]) rtol=.01
    @test norm(flame_residual!(similar(f.state),f),Inf) < 1e-8
    @test maximum(abs.(sum(mass_fractions(f);dims=1) .- 1)) < 1e-10
    npzwrite(joinpath(ARGS[2],"julia-"*label*".npz"),
        Dict("grid"=>f.grid,"T"=>temperature(f),"Y"=>mass_fractions(f),"speed"=>[flame_speed(f)]))
end
end
