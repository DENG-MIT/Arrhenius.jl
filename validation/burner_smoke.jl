using Arrhenius, NPZ, LinearAlgebra, Test
BLAS.set_num_threads(1)
gas = CreateSolution(ARGS[1])
for fixed in (false,true)
    f = BurnerFlame(gas; T=373.,P=.05*Arrhenius.one_atm,mdot=.06,
        X=Dict("H2"=>1.5,"O2"=>1.,"AR"=>7.),width=.5)
    filename = fixed ? "burner-fixed.npz" : "burner.npz"
    oracle = npzread(joinpath(ARGS[2],filename))
    if fixed
        set_temperature_profile!(f,oracle["positions"],oracle["temperatures"])
    end
    elapsed = @elapsed solve!(f; slope=.05,curve=.1,loglevel=1)
    residual = flame_residual!(similar(f.state),f)
    println("BURNER fixed=",fixed," seconds=",elapsed," Tmax=",maximum(temperature(f)),
        " points=",length(f.grid)," residual=",norm(residual,Inf))
    @test f.converged
    @test norm(residual,Inf) < 1e-8
    @test maximum(temperature(f)) ≈ maximum(oracle["T"]) rtol=.01
    @test maximum(abs.(sum(mass_fractions(f),dims=1).-1)) < 1e-10
    @test all(abs.(f.state[end,:] .- .06) .< 1e-10)
    npzwrite("julia-"*filename,Dict("grid"=>f.grid,"T"=>temperature(f),"Y"=>mass_fractions(f)))
end
