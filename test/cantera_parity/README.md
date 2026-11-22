# Cross-residual tool

`xresid.py` compares native and Cantera residuals at two states and checks the state-transfer round trips. The native subprocess is `xresid_native.jl`; both are kept in this directory. `snapshots/xresid_adapter.py` imports the local `xresid.py`, so no separate checkout or path variable is needed. Python requires NumPy and Cantera 4; Julia resolves Arrhenius from the repository root. `XRESID_JULIA` can optionally select a Julia executable.

Run the compact h2o2 regression from the repository root:

```sh
python test/cantera_parity/xresid_h2o2.py --kinds free
python test/cantera_parity/xresid_h2o2.py
```

The first command is a single-kind smoke check; the second covers free, burner, fixed-temperature, and stagnation flames. For recorded states, run `xresid.py SPEC STATE_A.npz STATE_B.npz` or use the adapter with the same arguments.

## Free flame diffusion gradient basis

`flux_gradient_mass.jl` ports [`TestFreeFlame.test_flux_gradient_mass`](https://github.com/Cantera/cantera/blob/726522be4e2a13454d8415b7ef799d621f665cf3/test/python/test_onedim.py) from Cantera 4.0.0a2. It runs in `Pkg.test()` with frozen [mechanism and transport data](../data/cantera_flux_gradient_mass/) and needs no Cantera runtime. It solves on a molar gradient basis, then continues the same flame on a mass gradient basis. The speeds must differ by no more than 10% of the mass-basis speed and must be unequal.

Run the single case from the repository root:

```sh
julia --project=. test/cantera_parity/flux_gradient_mass.jl
```

## Burner flame blowoff

`burner_blowoff.jl` ports [`TestBurnerFlame.test_blowoff`](https://github.com/Cantera/cantera/blob/726522be4e2a13454d8415b7ef799d621f665cf3/test/python/test_onedim.py#L2001) from Cantera 4.0.0a2. It reuses the frozen [10-species mechanism and transport data](../data/cantera_flux_gradient_mass/) and records the Cantera inlet and outlet values. `Pkg.test()` checks the original within-flame endpoint criteria for temperature, velocity, and each species mass fraction for both the recorded reference and the native solution. The test needs no Cantera runtime.

```sh
julia --project=. test/cantera_parity/burner_blowoff.jl
```
