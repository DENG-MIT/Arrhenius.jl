using Arrhenius, NPZ, LinearAlgebra
BLAS.set_num_threads(1)
gas = CreateSolution(ARGS[1])
f = FreeFlame(gas; X=Dict("H2"=>1.1,"O2"=>1.,"AR"=>5.))
println("Equilibrium initialization T=",maximum(temperature(f)))
if length(ARGS) > 1
    # Diagnostic only: validate a discretized residual near an independent
    # solution. This is never the initialization used by the public solver.
    ref = npzread(ARGS[2])
    f.grid = ref["grid"]
    f.state = vcat(transpose(ref["T"]/1000),ref["Y"],transpose(ref["mdot"]))
    f.anchor = argmin(abs.(ref["T"] .- ref["fixed_temperature"][1]))
    f.fixed_temperature = ref["T"][f.anchor]
    f.inlet_Y .= ref["inlet_Y"]
    println("Reference residual=",norm(flame_residual!(similar(f.state),f),Inf))
end
elapsed = @elapsed solve!(f; loglevel=1)
r = flame_residual!(similar(f.state),f)
println("RESULT seconds=",elapsed," speed=",flame_speed(f)," Tmax=",maximum(temperature(f)),
    " points=",length(f.grid)," residual=",norm(r,Inf)," minY=",minimum(mass_fractions(f)))
npzwrite("flame-smoke.npz",Dict("grid"=>f.grid,"T"=>temperature(f),"Y"=>mass_fractions(f),
    "velocity"=>velocity(f),"residual"=>r))
