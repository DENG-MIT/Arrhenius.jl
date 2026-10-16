# Methane/air burner with the original flame_fixed_T.py temperature profile.
# julia --project=. example/flames/flame_fixed_T.jl /path/to/stock/gri30.yaml [output-directory]
# Prepare its chemistry and multicomponent sidecars before running.
using Arrhenius
include("source_flame_sequence.jl")

isempty(ARGS) && throw(ArgumentError("supply prepared stock Cantera gri30.yaml (53 species, 325 reactions)"))
mechanism=ARGS[1]
output=length(ARGS)>1 ? ARGS[2] : mktempdir()
mkpath(output)
gas=CreateSolution(mechanism)
data=MultiTransportData(mechanism*".multicomponent.npz",gas;mechanism)
profile=source_flame_temperature_profile()
function write_fixed_stage(flame,mode)
    label=mode=="mole" ? "mixture" : "multicomponent"
    save_flame(joinpath(output,label*".npz"),flame)
    if mode=="multi"
        save_flame(joinpath(output,"flame-fixed-T.csv"),flame;basis=:mole)
        println("Solved ",length(flame.grid)," points; inlet velocity ",velocity(flame)[1]," m/s")
    end
end
run_source_sequence(gas,data,"fixed",profile;after_stage=write_fixed_stage)
println("Results: ",output)
