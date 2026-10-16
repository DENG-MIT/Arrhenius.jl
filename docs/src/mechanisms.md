# Mechanisms and preprocessing

CreateSolution("mechanism.yaml") loads species thermodynamics, reaction data, and transport fits from the YAML file and its .yaml.npz sidecar. The mechanism exporter uses Cantera's YAML parser during preprocessing; Cantera is not required when running Julia calculations.

The sidecar exporter supports elementary and three-body Arrhenius reactions, Lindemann and Troe falloff, pressure-dependent Arrhenius (PLOG) rates, Blowers–Masel rates, explicit reaction orders, and supported NASA7, NASA9, Shomate, and constant-cp thermodynamic models. Unsupported rate models are rejected during export.

YAML phase definitions may select local or imported species and reaction sections. Imported files are resolved relative to the mechanism file; optional data_paths can provide additional search directories:

~~~julia
gas = CreateSolution("mechanism.yaml"; data_paths=["shared-data"])
~~~

Regenerate the sidecar after changing the mechanism or any imported YAML file. For multicomponent diffusion or Soret diffusion, export matching collision-integral data with mechanism/export_multicomponent.py; see [flame transport](flames.md#transport).
