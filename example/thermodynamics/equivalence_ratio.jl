using Arrhenius
using Printf

"""
    equivalence_ratio_calculation(gas)

Reproduce every calculation of Cantera's `thermo/equivalenceRatio.py` example
(commit 726522be4e2a13454d8415b7ef799d621f665cf3) in source order with the
native stateless API and return all results as a named tuple. `X_`/`Y_` arrays
hold mole/mass fractions in `gas.species_names` order; scalars are
dimensionless except `_K` (kelvin) and `_Pa` (pascal) fields.

The native `Solution` carries no state, so the source's implicit gas state is
made explicit: mixtures created before the constant-HP equilibration keep the
default 300 K / 1 atm, and every mixture created afterwards retains the burnt
temperature (`T_post_burnt_K`).
"""
function equivalence_ratio_calculation(gas)
    air = "O2:0.21,N2:0.79"
    mole_fraction(name,X) = mole_fractions(gas,X)[findfirst(==(name),gas.species_names)]
    mass_fraction(name,X) = mass_fractions(gas,X)[findfirst(==(name),gas.species_names)]

    # Stoichiometric CH4/air mixture, fuel and oxidizer on a mole basis.
    X_stoich_mole = set_equivalence_ratio(gas,1.0;fuel="CH4:1",oxidizer=air)

    # Same equivalence ratio with compositions given on a mass basis.
    X_stoich_mass = set_equivalence_ratio(gas,1.0;fuel="CH4:1",
        oxidizer="O2:0.233,N2:0.767",basis=:mass)

    # Equivalence ratio of the current mixture; amounts need not sum to one.
    phi_mass_basis = equivalence_ratio(gas,X_stoich_mass;fuel="CH4:1",
        oxidizer="O2:233,N2:767",basis=:mass)

    # Without fuel/oxidizer: assume all C/H/S come from fuel, all O from oxidizer.
    phi_assumed_origin = equivalence_ratio(gas,X_stoich_mass)

    # Mixture fraction is always kg fuel / (kg fuel + kg oxidizer).
    Z_default = mixture_fraction(gas,X_stoich_mass;fuel="CH4:1",oxidizer=air)
    Z_bilger = mixture_fraction(gas,X_stoich_mass;fuel="CH4:1",oxidizer=air,element="Bilger")
    Z_carbon = mixture_fraction(gas,X_stoich_mass;fuel="CH4:1",oxidizer=air,element="C")

    # Pure-methane fuel with air: Z equals the CH4 mass fraction.
    Y_CH4_stoich = mass_fraction("CH4",X_stoich_mass)

    # Set a mixture holding 5.5 mass-% fuel.
    X_Z055 = set_mixture_fraction(gas,0.055;fuel="CH4:1",oxidizer=air)
    Y_CH4_Z055 = mass_fraction("CH4",X_Z055)

    # phi and Z are invariant to reaction progress; compositions may be dicts.
    fuel_dict = Dict("CH4"=>1)
    X_fresh_burnt_case = set_equivalence_ratio(gas,1.0;fuel=fuel_dict,oxidizer=air)
    burnt = equilibrate(gas;T=300.0,P=one_atm,X=X_fresh_burnt_case,mode=:HP)
    T_burnt_K, P_burnt_Pa = burnt.T, burnt.P
    X_burnt, Y_burnt = burnt.X, burnt.Y
    phi_burnt = equivalence_ratio(gas,X_burnt;fuel=fuel_dict,oxidizer=air)
    Z_burnt = mixture_fraction(gas,X_burnt;fuel=fuel_dict,oxidizer=air)

    # Consistent arbitrary fuel/oxidizer; the source retains the burnt
    # temperature for this and all following mixtures.
    fuel_arbitrary = "CH4:1,O2:0.01,CO:0.05,N2:0.1"
    oxidizer_arbitrary = "O2:0.2,N2:0.8,CO2:0.05,CH4:0.01"
    X_arbitrary = set_equivalence_ratio(gas,2.5;fuel=fuel_arbitrary,oxidizer=oxidizer_arbitrary)
    phi_arbitrary = equivalence_ratio(gas,X_arbitrary;fuel=fuel_arbitrary,oxidizer=oxidizer_arbitrary)
    phi_arbitrary_assumed = equivalence_ratio(gas,X_arbitrary)

    # Dilute a phi=2 H2/O2 mixture so the product holds 30 mol-% H2O.
    X_diluted_H2O = set_equivalence_ratio(gas,2.0;fuel="H2:1",oxidizer="O2:1",
        diluent="H2O",fraction=(diluent=0.3,))
    X_H2O_diluted = mole_fraction("H2O",X_diluted_H2O)
    H2_O2_mole_ratio = mole_fraction("H2",X_diluted_H2O)/mole_fraction("O2",X_diluted_H2O)

    # Same mixture diluted on a mass basis so the fuel mass fraction is 0.1.
    X_diluted_mass = set_equivalence_ratio(gas,2.0;fuel="H2",oxidizer="O2",
        diluent="CO2:0.5,H2O:0.5",fraction=(fuel=0.1,),basis=:mass)
    Y_H2_diluted = mass_fraction("H2",X_diluted_mass)

    # Equivalence ratio of the diluted mixture, ignoring diluent species.
    phi_include_species = equivalence_ratio(gas,X_diluted_mass;fuel="H2",oxidizer="O2",
        include_species=["H2","O2"])

    # Treat the diluent as part of the fuel stream instead.
    fuel_diluted = "H2:0.5,H2O:0.5"
    X_diluted_fuel = set_equivalence_ratio(gas,2.0;fuel=fuel_diluted,oxidizer=air)
    phi_diluted_fuel = equivalence_ratio(gas,X_diluted_fuel;fuel=fuel_diluted,oxidizer=air)

    return (;species_names=gas.species_names,
        T_fresh_K=300.0,P_fresh_Pa=one_atm,
        X_stoich_mole,X_stoich_mass,
        phi_mass_basis,phi_assumed_origin,
        Z_default,Z_bilger,Z_carbon,Y_CH4_stoich,
        X_Z055,Y_CH4_Z055,
        X_fresh_burnt_case,T_burnt_K,P_burnt_Pa,X_burnt,Y_burnt,phi_burnt,Z_burnt,
        T_post_burnt_K=T_burnt_K,
        X_arbitrary,phi_arbitrary,phi_arbitrary_assumed,
        X_diluted_H2O,X_H2O_diluted,H2_O2_mole_ratio,
        X_diluted_mass,Y_H2_diluted,phi_include_species,
        X_diluted_fuel,phi_diluted_fuel)
end

function main()
    gas = CreateSolution(joinpath(@__DIR__,"..","..","mechanism","gri30.yaml"))
    r = equivalence_ratio_calculation(gas)
    f3(x) = @sprintf("%.3f",x)
    println("phi = ",f3(r.phi_mass_basis))
    println("phi = ",f3(r.phi_assumed_origin))
    println("Z = ",f3(r.Z_default))
    println("Z(Bilger mixture fraction) = ",f3(r.Z_bilger))
    println("Z(mixture fraction based on C) = ",f3(r.Z_carbon))
    println("mass fraction of CH4 = ",f3(r.Y_CH4_stoich))
    println("mass fraction of CH4 = ",f3(r.Y_CH4_Z055))
    println("adiabatic equilibrium T = ",round(r.T_burnt_K;digits=3)," K")
    println("phi(burnt) = ",f3(r.phi_burnt))
    println("Z(burnt) = ",f3(r.Z_burnt))
    println("phi = ",f3(r.phi_arbitrary))
    println("phi = ",f3(r.phi_arbitrary_assumed))
    println("mole fraction of H2O = ",f3(r.X_H2O_diluted))
    println("ratio of H2/O2: ",f3(r.H2_O2_mole_ratio))
    println("mass fraction of H2 = ",f3(r.Y_H2_diluted))
    println("phi = ",f3(r.phi_include_species))
    println("phi = ",f3(r.phi_diluted_fuel))
    return r
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
