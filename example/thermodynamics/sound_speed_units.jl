include("sound_speed.jl")

"The source's 25 temperatures in °F and three sound speeds in ft/s."
function sound_speed_units_calculation(gas)
    fahrenheit = collect(range(80.,4880.;length=25))
    kelvin = (fahrenheit .- 32).*(5/9) .+ 273.15
    result = sound_speed_calculation(gas;temperatures=kelvin)
    data = copy(result.data)
    data[:,1] .= fahrenheit
    data[:,2:4] ./= .3048
    return (;data,final_states=result.final_states)
end

if abspath(PROGRAM_FILE) == @__FILE__
    gas = CreateSolution(joinpath(@__DIR__,"..","..","mechanism","gri30.yaml"))
    println("T [°F], equilibrium [ft/s], frozen [ft/s], frozen at perturbed equilibrium [ft/s]")
    display(sound_speed_units_calculation(gas).data)
end
