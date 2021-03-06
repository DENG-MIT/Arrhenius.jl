module Arrhenius
    using YAML
    using NPZ
    using LinearAlgebra
    using SparseArrays

    include("Constants.jl")
    include("DataStructure.jl")
    include("Solution.jl")
    include("Magic.jl")
    include("Thermo.jl")
    include("Kinetics.jl")
    include("Transport.jl")
end