import cantera as ct
import torch
import numpy as np

import reactorch as rt

import json

mech_yaml = 'JP10skeletal.yaml'
# mech_yaml = 'h2o2.yaml'

sol = rt.Solution(mech_yaml=mech_yaml, device=torch.device('cpu'),
                  vectorize=True,
                  is_clip=False, is_norm=False, is_wdot_vec=False)

gas = sol.gas

gas.TPY = 1200, ct.one_atm, np.ones(gas.n_species)/gas.n_species

molecular_weights = gas.molecular_weights.tolist()
reactant_stoich_coeffs = gas.reactant_stoich_coeffs()
product_stoich_coeffs = gas.product_stoich_coeffs()
Arrhenius_coeffs = sol.Arrhenius_coeffs
efficiencies_coeffs = sol.efficiencies_coeffs

np.savez(mech_yaml,
         molecular_weights=molecular_weights,
         reactant_stoich_coeffs=reactant_stoich_coeffs,
         product_stoich_coeffs=product_stoich_coeffs,
         reactant_orders=sol.reactant_orders,
         is_reversible=sol.is_reversible,
         Arrhenius_coeffs=Arrhenius_coeffs,
         efficiencies_coeffs=efficiencies_coeffs,
         Arrhenius_A0=sol.Arrhenius_A0,
         Arrhenius_b0=sol.Arrhenius_b0,
         Arrhenius_Ea0=sol.Arrhenius_Ea0,
         Troe_A=sol.Troe_A,
         Troe_T1=sol.Troe_T1,
         Troe_T2=sol.Troe_T2,
         Troe_T3=sol.Troe_T3
         )
