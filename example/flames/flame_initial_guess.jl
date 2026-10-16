# Freely-propagating premixed methane/air flame: saving, loading, and
# restarting 1D calculations with different initial guesses.
# Native port of https://cantera.org/dev/examples/python/onedim/flame_initial_guess.html
# julia --project=. example/flames/flame_initial_guess.jl /path/to/gri30.yaml [output-directory]
using Arrhenius

isempty(ARGS) && throw(ArgumentError("supply stock Cantera gri30.yaml"))
mechanism = ARGS[1]
output = length(ARGS) > 1 ? ARGS[2] : mktempdir()
mkpath(output)

gas = CreateSolution(mechanism)

# Simulation parameters
p = one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
reactants = "CH4:0.45, O2:1.0, N2:3.76"
width = 0.03  # m

function describe(flame)
    oh_index = findfirst(==("OH"), flame.gas.species_names)
    oh_peak = flame.grid[argmax(mole_fractions(flame)[oh_index, :])]
    println("\nFlame speed                    = ",
        round(flame_speed(flame) * 100, digits=2), " cm/s")
    println("Maximum temperature            = ",
        round(maximum(temperature(flame)), digits=0), " K")
    println("Peak OH mole fraction location = ", round(oh_peak * 100, digits=2), " cm")
    println("Solved with ", length(flame.grid), " grid points\n")
end

# %%
# Initialization
f = FreeFlame(gas; T=Tin, P=p, X=reactants, width, discretization=:conservative)
# These refine criteria are for a relatively quick convergence for the example
# and may not be suitable for a publication-quality flame speed calculation.
# The explicit :conservative discretization converges every restart below;
# discretization=:finite_difference reproduces Cantera's own converged
# methane flame speed even more closely (34.45 cm/s and 117 points here).
solve!(f; ratio=3, slope=0.1, curve=0.2)
describe(f)

# %%
# Save the flame as a native restart snapshot and as a CSV table
# (Cantera's example saves YAML/HDF containers plus a mole-basis CSV; the
# native snapshot plays the role of the containers)
println("Save snapshot (npz)")
save_flame(joinpath(output, "flame.npz"), f)
println("Save CSV")
csv = joinpath(output, "flame.csv")
save_flame(csv, f; basis=:mole, overwrite=true)

# %%
# Restore the flame from the snapshot (Cantera: restore from YAML/HDF)
restored = FreeFlame(gas; T=Tin, P=p, X=reactants, width)
restore_flame!(restored, joinpath(output, "flame.npz"))
describe(restored)

# %%
# Load the initial guess from the CSV file directly
table = read_flame_csv(csv)
guess = FreeFlame(gas; T=Tin, P=p, X=reactants, width)
set_initial_guess!(guess, table)
describe(guess)

# %%
# Restart flame simulations with modified initial guesses
# Load the CSV table, removing half of the grid points
pruned = (; table.columns, data=table.data[1:2:end, :])
f2 = FreeFlame(gas; T=Tin, P=p, X=reactants, width)
set_initial_guess!(f2, pruned)
solve!(f2; ratio=3, slope=0.1, curve=0.2)

# %%
# We wouldn't expect the flame solutions to be exactly the same
describe(f2)
save_flame(joinpath(output, "flame_pruned.csv"), f2; basis=:mole, overwrite=true)

# %%
# Load the CSV table, removing half of the grid points and all but the first
# 20 species columns
pruned24 = (; columns=table.columns[1:24], data=table.data[1:2:end, 1:24])
f3 = FreeFlame(gas; T=Tin, P=p, X=reactants, width)
set_initial_guess!(f3, pruned24)
solve!(f3; ratio=3, slope=0.1, curve=0.2)

# %%
# We wouldn't expect the flame solutions to be exactly the same
describe(f3)
save_flame(joinpath(output, "flame_pruned_species.csv"), f3; basis=:mole, overwrite=true)

# %%
# Load the CSV table, removing half of the grid points, and raise the inlet
# temperature by 50 K
hot = FreeFlame(gas; T=Tin + 50, P=p, X=reactants, width)
set_initial_guess!(hot, pruned)
solve!(hot; ratio=3, slope=0.1, curve=0.2)

# %%
# We expect these flames to be different because we raised the temperature
describe(hot)
save_flame(joinpath(output, "flame_hot.csv"), hot; basis=:mole, overwrite=true)
println("Results: ", output)
