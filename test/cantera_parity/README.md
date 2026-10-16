# Cross-residual tool

`xresid.py` compares native and Cantera residuals at two states and checks the state-transfer round trips. The native subprocess is `xresid_native.jl`; both are kept in this directory. `snapshots/xresid_adapter.py` imports the local `xresid.py`, so no separate checkout or path variable is needed. Python requires NumPy and Cantera 4; Julia resolves Arrhenius from the repository root. `XRESID_JULIA` can optionally select a Julia executable.

Run the compact h2o2 regression from the repository root:

```sh
python test/cantera_parity/xresid_h2o2.py --kinds free
python test/cantera_parity/xresid_h2o2.py
```

The first command is a single-kind smoke check; the second covers free, burner, fixed-temperature, and stagnation flames. For recorded states, run `xresid.py SPEC STATE_A.npz STATE_B.npz` or use the adapter with the same arguments.
