# Low-pressure hydrogen/oxygen/argon burner, mixture then multicomponent transport.
# julia --project=. example/flames/burner_flame.jl /path/to/stock/h2o2.yaml [output-directory]
# Prepare its chemistry and multicomponent sidecars before running.
using Arrhenius
include("source_flame_sequence.jl")

isempty(ARGS) && throw(ArgumentError("supply prepared stock Cantera h2o2.yaml (10 species, 29 reactions)"))
mechanism=ARGS[1]
output=length(ARGS)>1 ? ARGS[2] : mktempdir()
mkpath(output)
gas=CreateSolution(mechanism)
data=MultiTransportData(mechanism*".multicomponent.npz",gas;mechanism)
function write_burner_stage(flame,mode)
    label=mode=="mole" ? "mixture" : "multicomponent"
    println(label," peak temperature: ",maximum(temperature(flame))," K")
    save_flame(joinpath(output,label*".npz"),flame)
    if mode=="multi"
        save_flame(joinpath(output,"burner.csv"),flame)
        save_flame(joinpath(output,"burner.npz"),flame)
    end
end
run_source_sequence(gas,data,"burner";after_stage=write_burner_stage)
println("Results: ",output)
