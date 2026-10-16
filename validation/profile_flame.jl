using Arrhenius, LinearAlgebra, Profile
BLAS.set_num_threads(1)
gas = CreateSolution(ARGS[1])
function run(gas)
    if length(ARGS)>1 && ARGS[2] == "burner"
        f = BurnerFlame(gas;T=373.,P=.05one_atm,mdot=.06,width=.5,X="H2:1.5,O2:1,AR:7")
        solve!(f;slope=.05,curve=.1)
    else
        f = FreeFlame(gas; X=Dict("H2"=>1.1,"O2"=>1.,"AR"=>5.))
        solve!(f)
    end
end
run(gas)
Profile.clear()
@profile for i in 1:10
    run(gas)
end
Profile.print(format=:flat, sortedby=:count, mincount=50)
