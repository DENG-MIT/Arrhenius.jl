module Arrhenius
    include("Constants.jl")
    include("Solution.jl")

    # Initially, we don't have to do it here, all of those infomation can be interpreted from Cantera and ReacTorch

    # YAML reads in the species and thermo data
    # List: species (file)
    # Array: NASA polynomials (file)
    # Array: Arrhenius parameters (file)
    # List: three body reaction (file)
    # List: falloff reaction (file)
    # Array: reactants orders (file)
    # Array: stoichiometric cofficients (file)
    # Array: molecular weight
end
