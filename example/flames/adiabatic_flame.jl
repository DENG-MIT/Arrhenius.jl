# Hydrogen/oxygen/argon flame: mixture, mixture/Soret, multi, and multi/Soret.
# julia --project=. example/flames/adiabatic_flame.jl /path/to/stock/h2o2.yaml [output-directory]
# Prepare its chemistry and multicomponent sidecars before running.
using Arrhenius
include("source_flame_sequence.jl")

isempty(ARGS) && throw(ArgumentError("supply prepared stock Cantera h2o2.yaml (10 species, 29 reactions)"))
mechanism=ARGS[1]
output=length(ARGS)>1 ? ARGS[2] : mktempdir()
mkpath(output)
gas=CreateSolution(mechanism)
data=MultiTransportData(mechanism*".multicomponent.npz",gas;mechanism)
labels=Dict("mass"=>"mixture","mass-soret"=>"mixture-soret",
    "multi"=>"multicomponent","multi-soret"=>"multicomponent-soret")
function write_adiabatic_stage(flame,mode)
    label=labels[mode]
    println(label," speed: ",flame_speed(flame)," m/s")
    save_flame(joinpath(output,label*".npz"),flame)
    mode=="multi-soret" && save_flame(joinpath(output,"adiabatic-flame.csv"),flame;basis=:mole)
end
run_source_sequence(gas,data,"free";after_stage=write_adiabatic_stage)
println("Results: ",output)
