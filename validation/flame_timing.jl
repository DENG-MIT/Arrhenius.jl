using Arrhenius, LinearAlgebra, NPZ
BLAS.set_num_threads(1)
gas = CreateSolution(ARGS[1])
for trial in 1:5
    local f
    elapsed = @elapsed begin
        f = FreeFlame(gas; X=Dict("H2"=>1.1,"O2"=>1.,"AR"=>5.))
        solve!(f)
    end
    println("trial=",trial," seconds=",elapsed," speed=",flame_speed(f)," points=",length(f.grid))
end
