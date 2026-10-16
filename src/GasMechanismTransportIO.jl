# Transport is independent of gas kinetics. YAML transport metadata is retained
# by species definitions; legacy polynomial fits remain an optional sidecar.
function _read_gas_transport(arrays, sidecar, reader)
    reader===:native && sidecar===nothing && return arrays
    reader===:sidecar && return arrays
    sidecar===nothing && return arrays
    merged=Dict{String,Any}(String(k)=>v for (k,v) in arrays)
    for key in ("transport_model_utf8","species_viscosities_poly","thermal_conductivity_poly","binary_diff_coeffs_poly")
        haskey(sidecar,key) && (merged[key]=sidecar[key])
    end
    merged
end

