# API

The API reference is generated from exported Arrhenius.jl docstrings. It includes constructors, solver functions, thermodynamic and transport methods, and diagnostics.

~~~@autodocs
Modules = [Arrhenius]
Private = false
Order = [:constant, :type, :function]
~~~

For usage paths, see [reactors](reactors.md), [flames](flames.md), and the [examples index](examples.md).


## Plasma module

Plasma functionality is organized under `Arrhenius.Plasma`. The EEDF solver,
plasma mechanism and state readers, plasma rate and thermochemistry functions,
and plasma reactor types are available there. Existing top-level names remain
available for compatibility, for example `Arrhenius.solve_eedf` and
`Arrhenius.Plasma.solve_eedf` refer to the same function. Plasma composes the
core ideal-gas thermochemistry, kinetics and reactor APIs. Surface kinetics and
coverage thermochemistry are available through `Arrhenius.Plasma` as well as
`Arrhenius`; their implementations remain owned by the core module.
