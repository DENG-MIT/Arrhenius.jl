# Cantera auto Jacobian test fixture

The tracked `free_auto_stagni.npz` contains Cantera's 9-point FreeFlame initial state, residuals, and full Jacobian from Cantera 4.0.0a2. `free_auto_stagni_source.npz` contains the Cantera `net_production_rates_ddCi` source Jacobian transformed to fixed-pressure mass-fraction coordinates for the same interior nodes.

To regenerate from the committed spec using the existing dump harness and Cantera environment:

```sh
# activate a Python environment with Cantera at the pinned version
python test/cantera_parity/ct_dump.py test/data/cantera/free_auto_stagni.json test/data/cantera/generated state
mv test/data/cantera/generated/ct_state.npz test/data/cantera/free_auto_stagni.npz
mv test/data/cantera/generated/ct_state_meta.json test/data/cantera/free_auto_stagni_meta.json
python test/cantera_parity/generate_auto_species_fixture.py
```

`ct_dump.py` uses the repository mechanism path in the spec and evaluates the Cantera auto Jacobian at the initial 9-point state. The tests compare the native source-derivative routine against the independently generated Cantera fixture.
