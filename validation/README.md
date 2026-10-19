# Validation entry points

This directory contains opt-in validation and benchmark drivers. Package tests run with `Pkg.test()`; the drivers below prepare or compare focused cases and may require Cantera 4.0 for reference data. Mechanism files used by native examples are kept in `mechanism/`.

| Family | Entry points |
| --- | --- |
| Flames | [`flame_benchmarks.py`](flame_benchmarks.py), [`flame_benchmarks.jl`](flame_benchmarks.jl), [`burner_oracle.py`](burner_oracle.py), [`burner_smoke.jl`](burner_smoke.jl), [`counterflow_native.jl`](counterflow_native.jl) |
| Thermochemistry | [`blowers_masel_cases.py`](blowers_masel_cases.py), [`blowers_masel_cases.jl`](blowers_masel_cases.jl), [`equilibrium_cases.py`](equilibrium_cases.py), [`equilibrium_cases.jl`](equilibrium_cases.jl), [`species_thermo_cases.py`](species_thermo_cases.py), [`species_thermo_cases.jl`](species_thermo_cases.jl) |
| Reactors | [`reactor_cases.py`](reactor_cases.py), [`reactor_cases.jl`](reactor_cases.jl), [`reactor_network_cases.py`](reactor_network_cases.py), [`reactor_network_cases.jl`](reactor_network_cases.jl), [`ic_engine_case.py`](ic_engine_case.py), [`ic_engine_case.jl`](ic_engine_case.jl) |
| Mixing and walls | [`mixing_case.jl`](mixing_case.jl), [`moving_wall_cases.py`](moving_wall_cases.py), [`moving_wall_cases.jl`](moving_wall_cases.jl) |
| Transport | [`parallel_transport_cases.py`](parallel_transport_cases.py), [`parallel_transport_cases.jl`](parallel_transport_cases.jl) |
| Real gas and equations of state | [`real_gas_cases.py`](real_gas_cases.py), [`real_gas_cases.jl`](real_gas_cases.jl), [`equations_of_state_cases.py`](equations_of_state_cases.py), [`equations_of_state_cases.jl`](equations_of_state_cases.jl) |
| Water | [`pure_water_cases.py`](pure_water_cases.py), [`pure_water_cases.jl`](pure_water_cases.jl) |
| Surfaces | [`surface_cases.py`](surface_cases.py), [`surface_flow_cases.py`](surface_flow_cases.py), [`surface_reactor_cases.py`](surface_reactor_cases.py) |
| Other focused cases | [`custom2_case.jl`](custom2_case.jl), [`equivalence_ratio_case.py`](equivalence_ratio_case.py), [`isentropic_cases.py`](isentropic_cases.py), [`multicomponent_reference.py`](multicomponent_reference.py) |

Most Python case scripts create reference inputs or compare results; their paired Julia scripts exercise the native implementation. See each script's usage text for its arguments. Focused regression checks in this directory are separate from the package test suite.
