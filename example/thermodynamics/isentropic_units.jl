include("isentropic.jl")

"The ten source-example nozzle states, with temperature in K and area in m²."
function isentropic_units_calculation(gas)
    result = isentropic_calculation(gas;points=10)
    data = copy(result.data)
    data[:,3] .*= 1200. # 2160 °R = 1200 K
    return (;data,states=result.states,throat_area=result.throat_area)
end

if abspath(PROGRAM_FILE) == @__FILE__
    gas = CreateSolution(joinpath(@__DIR__,"..","..","mechanism","gri30.yaml"))
    println("area ratio, Mach number, temperature [K], pressure ratio")
    display(isentropic_units_calculation(gas).data)
end
