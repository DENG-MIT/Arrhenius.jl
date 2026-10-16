# Boltzmann plasma thermochemistry

`PlasmaMechanism` reads the gas and charged-species NASA7 data, imported reaction sections, and collision tables of a `Boltzmann-two-term` YAML phase. Supported rates include reversible Arrhenius reactions, simple third bodies, Lindemann/Troe falloff, and irreversible two-temperature, Chebyshev, and electron-collision rates.

The gas temperature, stored electron temperature, electric field, and electron distribution are distinct state variables. Setting E/N stores the absolute field using the current density. `update_eedf!` updates the distribution and mobility while retaining the stored electron temperature. Gas-state and enthalpy setters retain the cached field, distribution, and mobility.

```julia
using Arrhenius
m = PlasmaMechanism("methane-plasma-pavan-2023.yaml";
    data_paths=["/path/to/thermodynamic-data"])
s = PlasmaState(m; temperature=300.0, pressure=101325.0,
    mole_fractions=Dict("CH4"=>0.095, "O2"=>0.19, "N2"=>0.715, "e"=>1e-11))
set_reduced_electric_field!(s, 190e-21 * exp(-32))
update_eedf!(s)
set_reduced_electric_field!(s, 190e-21)
update_eedf!(s)
rates = plasma_rates(s)
h = plasma_thermodynamics(s).h_mass
set_plasma_state!(s; temperature=330.0)
set_plasma_enthalpy!(s, h; pressure=101325.0)
```

The [state validation results](results/cantera4_wsl_pulse_thermochemistry.json) compare the methane pulse mechanism's 71 species and 440 reactions at five states against Cantera 4.0.0a2. They cover gas and electron thermodynamics, all rate families, density, field/EEDF caching, and enthalpy inversion. Tiny net production rates are also checked against high-precision sums of independently parsed net stoichiometry: separate production/destruction accumulation can leave floating-point cancellation residuals.

`PlasmaEnergyReactor` supplies a closed, constant-pressure energy equation for Boltzmann phases. Its ODE state contains mass, total enthalpy, and all species mass fractions. Gas temperature is recovered from enthalpy while electron temperature remains fixed. Joule heating changes total enthalpy; the RHS retains the field, mobility, and collision-rate cache between explicit EEDF updates.

```julia
r = PlasmaEnergyReactor(s; volume=1.0)
rhs = reactor_rhs(r)
u = reactor_state(r)
du = similar(u)
rhs(du, u)
properties = reactor_properties(rhs, u)
update_eedf!(rhs, u; reduced_field=190e-21)
```

Each `reactor_rhs(r)` or `reactor_problem(r, tspan)` starts from the reactor snapshot. For piecewise field histories, retain the prepared RHS and accepted terminal state between intervals, update that RHS's EEDF, and restart the integrator's history and caches.

The [energy-reactor validation results](results/cantera4_wsl_pulse_energy.json) compare the 73-state RHS, thermodynamics, and cached EEDF at six prescribed methane-pulse states against a converged Cantera reference. Species sources use independently accumulated net stoichiometry. The optional [integration test](../test/plasma_energy_integration.jl) checks two inert charged-mixture heating intervals against a closed-form solution:

```sh
julia --project=example/reactors test/plasma_energy_integration.jl
```

The [complete native methane pulse](results/cantera4_wsl_nanosecond_pulse.json) passes all 664 trajectory and snapshot checks, including its 901 output times and 91 field-update intervals. The [repeated-run benchmark](results/cantera4_wsl_nanosecond_pulse_repeated.json) passes 7,550 accuracy checks and measures 0.633 times Cantera's speed in WSL (0.507 s native versus 0.321 s Cantera), below the 0.95 target. `PlasmaReactor` remains an isothermal model and does not accept Boltzmann phases requiring energy coupling.
