# Critical parameters and compressibility of eight pure-fluid models.
# Equivalent calculation: Cantera thermo/critical_properties.py.
using Arrhenius

const CRITICAL_FLUID_NAMES = ["water", "nitrogen", "methane", "hydrogen",
    "oxygen", "carbon dioxide", "heptane", "HFC-134a"]

function critical_properties_calculation(fluids=CRITICAL_FLUID_NAMES)
    output = Matrix{Float64}(undef, 5, length(fluids))
    for (column, fluid) in enumerate(fluids)
        point = critical_properties(fluid)
        output[1, column] = point.T
        output[2, column] = point.P
        output[3, column] = point.density
        output[4, column] = point.molecular_weight
        output[5, column] = point.Z
    end
    return output
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = critical_properties_calculation()
    println("fluid: Tc [K], Pc [Pa], density [kg/m³], molecular weight [kg/kmol], Zc")
    for (column, fluid) in enumerate(CRITICAL_FLUID_NAMES)
        println(fluid, ": ", join(result[:, column], ", "))
    end
end
