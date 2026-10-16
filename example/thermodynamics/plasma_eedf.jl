# Standalone air EEDF at 200 Td using the Phelps collision dataset.
# julia --project=. example/thermodynamics/plasma_eedf.jl air-plasma-Phelps.yaml [eedf.csv]
using Arrhenius

function air_eedf(path)
    model = read_eedf_model(path)
    Set(model.target_names) == Set(["N2", "O2"]) ||
        throw(ArgumentError("this example requires the N2/O2 Phelps air model"))
    state = EEDFState(model; T=300., P=one_atm,
        mole_fractions=Dict("N2"=>.79, "O2"=>.21, "N2+"=>1e-10, "Electron"=>1e-10),
        molecular_weights=Dict("N2"=>28.014, "O2"=>31.998), reduced_field=200e-21)
    return solve_eedf(model, state)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (1,2) || error("provide air-plasma-Phelps.yaml and optional output CSV")
    result = air_eedf(ARGS[1])
    destination = length(ARGS)==2 ? ARGS[2] : "eedf.csv"
    open(destination, "w") do io
        println(io, "energy_eV,eedf_eV_to_minus_3_over_2")
        for (energy, value) in zip(result.edges, result.edge_eedf)
            println(io, energy, ',', value)
        end
    end
    println("Converged in ", result.iterations, " iterations; mobility = ",
        result.mobility, " m²/(V·s)")
end
