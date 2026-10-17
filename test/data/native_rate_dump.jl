include(joinpath(@__DIR__, "..", "..", "src", "Arrhenius.jl"))
using .Arrhenius

const STATES = ((500.0, 0.1 * Arrhenius.one_atm),
                (1200.0, Arrhenius.one_atm),
                (2500.0, 10.0 * Arrhenius.one_atm))

println("# mechanism|T|P|reaction|forward|reverse")
for path in ARGS
    gas = Arrhenius.CreateSolution(path)
    X = fill(1.0 / gas.n_species, gas.n_species)
    for (T, P) in STATES
        rates = Arrhenius.reaction_rate_constants(gas; T, P, X)
        for i in eachindex(rates.forward)
            println(join((basename(path), repr(T), repr(P), i,
                          repr(rates.forward[i]), repr(rates.reverse[i])), Char(124)))
        end
    end
end
