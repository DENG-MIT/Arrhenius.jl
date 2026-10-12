# Diagnostic evaluation at reference profiles; not a solver benchmark or initialization.
using Arrhenius, NPZ, LinearAlgebra
gas = CreateSolution(ARGS[1])
data = MultiTransportData(ARGS[1]*".multicomponent.npz",gas)
f = FreeFlame(gas;X="H2:1.1,O2:1,AR:5",transport_model=:multicomponent,
    multicomponent_data=data,soret=true)
ref = npzread(ARGS[2])
f.grid = ref["grid"]
f.state = vcat(transpose(ref["T"]/1000),ref["Y"],transpose(ref["mdot"]))
f.anchor = argmin(abs.(ref["T"] .- only(ref["fixed_temperature"])))
f.fixed_temperature = ref["T"][f.anchor]
w = Arrhenius.FlameWorkspace(f)
r = flame_residual!(similar(f.state),f,f.state,w)
println("reference residual by row ",maximum(abs.(r);dims=2))
println("DT abs ",maximum(abs.(w.thermal_diffusion-ref["thermal_diffusion"])),
    " normalized ",maximum(abs.(w.thermal_diffusion-ref["thermal_diffusion"])) / maximum(abs.(ref["thermal_diffusion"])))
println("DT max native=",maximum(w.thermal_diffusion)," ref=",maximum(ref["thermal_diffusion"]))
println("flux sum ",maximum(abs.(sum(w.flux;dims=1))))
solve!(f;refine_grid=false)
println("Same grid solved speed=",flame_speed(f)," reference speed=",only(ref["speed"]))
