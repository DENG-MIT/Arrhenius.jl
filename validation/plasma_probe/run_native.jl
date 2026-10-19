using Arrhenius

input = ARGS[1]
csv = ARGS[2]
model = read_eedf_model(input)
state = EEDFState(model; T=300.0, P=133.32,
    mole_fractions=Dict("Ar" => 1.0), molecular_weights=Dict("Ar" => 39.948),
    reduced_field=100e-21, frequency=0.0)
result = solve_eedf(model, state)
result.converged || error("native EEDF did not converge")
open(csv, "w") do io
    println(io, "energy_eV,eedf")
    for (energy, value) in zip(result.edges, result.edge_eedf)
        println(io, energy, ",", value)
    end
end
println("NATIVE_OK cells=$(length(model.energy_edges)-1) iterations=$(result.iterations) converged=$(result.converged) mobility=$(result.mobility)")
