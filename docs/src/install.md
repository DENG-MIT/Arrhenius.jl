# Installation

Arrhenius.jl requires Julia. Install the package with Julia's package manager:

~~~julia
]
add https://github.com/DENG-MIT/Arrhenius.jl
~~~

A local checkout can be developed with Pkg.develop(path="."). The optional reactor example environment contains SciML solver dependencies used by the zero-dimensional examples.

## Mechanism files

CreateSolution reads a Cantera YAML mechanism and a same-name .yaml.npz sidecar. The repository includes prepared sidecars for its example mechanisms. For another mechanism, install Cantera and NumPy, then run:

~~~bash
python mechanism/export_sidecar.py path/to/mechanism.yaml
~~~

For multicomponent or Soret transport, also create a .multicomponent.npz file:

~~~bash
python mechanism/export_multicomponent.py path/to/mechanism.yaml path/to/mechanism.yaml.multicomponent.npz
~~~

Regenerate sidecars whenever the YAML or an imported file changes. See [mechanism support](mechanisms.md) for supported rate and thermodynamic models.
