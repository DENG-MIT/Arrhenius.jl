using Arrhenius

gas = CreateSolution(joinpath(@__DIR__,"..","..","mechanism","gri30.yaml"))
X = set_equivalence_ratio(gas,1.0;fuel="CH4",oxidizer="O2:1,N2:3.76")

# Each result conserves the pair selected by mode, starting from the same
# composition. SP/SV use a hotter inlet to keep the final temperature inside
# the thermodynamic polynomial ranges.
for mode in (:TP,:TV,:HP,:UV,:SP,:SV)
    T = mode in (:SP,:SV) ? 1500.0 : 300.0
    result = equilibrate(gas;T,P=one_atm,X,mode)
    println(mode,": T = ",round(result.T;digits=3)," K, P = ",
            round(result.P/one_atm;digits=6)," atm")
end
