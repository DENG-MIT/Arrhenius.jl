# Critical parameters of Cantera's TPX pure-fluid models:
# https://github.com/Cantera/cantera/tree/726522be4e2a13454d8415b7ef799d621f665cf3/src/tpx
# These are model parameters, including each model's molecular weight, rather
# than a replacement table of recommended experimental critical properties.
const _PURE_FLUID_CRITICAL_PARAMETERS = Dict(
    "water" => (647.286, 22.089e6, 317.0, 18.016),
    "nitrogen" => (126.200, 3.4e6, 314.03, 28.01348),
    "methane" => (190.555, 4.5988e6, 160.43, 16.04996),
    "hydrogen" => (32.938, 1.2838e6, 31.36, 2.0159),
    "oxygen" => (154.581, 5.0429e6, 436.15, 31.9994),
    "carbon dioxide" => (304.21, 7.38350e6, 464.00, 44.01),
    "heptane" => (537.68, 2.6199e6, 197.60, 100.20),
    # HFC134a::Tcrit/Pcrit/Vcrit expose the critical point, which differs from
    # the reducing constants Tc/Pc/Roc used inside that model's EOS.
    "hfc-134a" => (374.21, 4059280.0, 511.95, 102.032),
)

"""
    critical_properties(fluid::AbstractString)

Return the TPX model's critical temperature `T` [K], pressure `P` [Pa], density
`density` [kg/m³], molecular weight `molecular_weight` [kg/kmol], and critical
compressibility `Z = P*molecular_weight/(density*R*T)`.

Supported names (case insensitive) are `water`, `nitrogen`, `methane`,
`hydrogen`, `oxygen`, `carbon dioxide`, `heptane`, and `HFC-134a`. This function
provides critical properties; it does not construct an equation-of-state model.
For liquid/vapor water states, use [`PureWater`](@ref).
"""
function critical_properties(fluid::AbstractString)
    parameters = get(_PURE_FLUID_CRITICAL_PARAMETERS, lowercase(fluid), nothing)
    isnothing(parameters) && throw(ArgumentError("unsupported pure fluid: $fluid"))
    T, P, density, molecular_weight = parameters
    Z = P*molecular_weight/(density*R*T)
    return (; T, P, density, molecular_weight, Z)
end

export critical_properties
