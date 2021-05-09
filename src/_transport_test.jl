if pwd()[end - 2:end] != "src"
    cd("src")
end

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


```
Following codes are used during development phase only.
```
using LinearAlgebra
using DifferentialEquations
using Sundials
using Profile
using ForwardDiff
using SparseArrays
using PyCall
using Test

R = Arrhenius.R
one_atm = Arrhenius.one_atm

mech = "../mechanism/h2o2.yaml"
gas = Arrhenius.CreateSolution(mech);
ns = gas.n_species

Y0 = zeros(ns);
Y0[Arrhenius.species_index(gas, "H2")] = 0.1;
Y0[Arrhenius.species_index(gas, "O2")] = 0.2;
Y0[Arrhenius.species_index(gas, "AR")] = 0.7;

P = Arrhenius.one_atm
T = 1200.0
Y = Y0;
mean_MW = 1. / dot(Y, 1 ./ gas.MW)
ρ_mass = P / R / T * mean_MW
X =  Arrhenius.Y2X(gas, Y, mean_MW);
C =  Arrhenius.Y2C(gas, Y, ρ_mass);

η_mix, λ_mix, Dkm = Arrhenius.mix_trans(gas, P, T, X, mean_MW)

## Call Cantera
ct = pyimport("cantera")
ct_gas = ct.Solution(mech)

ct_gas.TPY = T, P, "H2:0.1,O2:0.2,AR:0.7"

@test η_mix ≈ ct_gas.viscosity rtol = 1.e-4

@test λ_mix ≈ ct_gas.thermal_conductivity rtol = 1.e-4

@test Dkm ≈ ct_gas.mix_diff_coeffs rtol = 1.e-3

println("Dkm from Arrhenius.jl | Cantera")
hcat(Dkm, ct_gas.mix_diff_coeffs)