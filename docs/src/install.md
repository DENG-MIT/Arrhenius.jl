# Installation

Arrhenius.jl requires Julia. Install the package with Julia's package manager:

~~~julia
]
add https://github.com/DENG-MIT/Arrhenius.jl
~~~

A local checkout can be developed with Pkg.develop(path="."). The optional reactor example environment contains SciML solver dependencies used by the zero-dimensional examples.

## Mechanism files

`CreateSolution` reads gas species, thermochemistry, reactions, and supported rate models directly from YAML. Python, Cantera, and a kinetics sidecar are not needed to load gas chemistry. A same-name `.yaml.npz` sidecar is optional and supplies fitted transport polynomials when available. For a mechanism that needs those transport fits, install Cantera and NumPy, then run:

~~~bash
python mechanism/export_sidecar.py path/to/mechanism.yaml
~~~

For multicomponent or Soret transport, also create a `.multicomponent.npz` file:

~~~bash
python mechanism/export_multicomponent.py path/to/mechanism.yaml path/to/mechanism.yaml.multicomponent.npz
~~~

Regenerate transport sidecars whenever the YAML or an imported file changes. See [mechanism support](mechanisms.md) for supported rate and thermodynamic models.
