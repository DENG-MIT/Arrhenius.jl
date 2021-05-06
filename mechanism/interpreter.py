import argparse

import cantera as ct
import numpy as np
import reactorch as rt
import torch


def fit_transport_data(gas, T_min=300, T_max=3000, n_points=2700):

    p = ct.one_atm
    comp = gas.species_name(1) + ":1.0"

    T_list = np.linspace(start=T_min, stop=T_max, num=n_points, endpoint=True)

    gas.TPX = gas.T, p, comp

    species_viscosities_list, species_viscosities_poly = rt.species_viscosities(
        gas, T_list, p, comp)
    binary_diff_coeffs_list, binary_diff_coeffs_poly = rt.binary_diff_coeffs(
        gas, T_list, p, comp)
    thermal_conductivity_list, thermal_conductivity_poly = rt.thermal_conductivity(
        gas, T_list, p, comp)

    return species_viscosities_poly, binary_diff_coeffs_poly, thermal_conductivity_poly


parser = argparse.ArgumentParser()

parser.add_argument('-i', '--yaml', default='gri30.yaml', 
                    required=False, help="yaml mech file")

args = parser.parse_args()

print(f'Processing mechanism {args.yaml}')

mech_yaml = args.yaml

sol = rt.Solution(mech_yaml=mech_yaml, device=torch.device('cpu'),
                  vectorize=True,
                  is_clip=False, is_norm=False, is_wdot_vec=False)

gas = sol.gas

gas.TPY = 1200, ct.one_atm, np.ones(gas.n_species)/gas.n_species

(species_viscosities_poly, binary_diff_coeffs_poly,
 thermal_conductivity_poly) = fit_transport_data(gas)


if len(sol.list_reaction_type4_Troe) == 0:
    np.savez(mech_yaml,
             molecular_weights=gas.molecular_weights.tolist(),
             reactant_stoich_coeffs=gas.reactant_stoich_coeffs(),
             product_stoich_coeffs=gas.product_stoich_coeffs(),
             reactant_orders=sol.reactant_orders,
             is_reversible=sol.is_reversible,
             Arrhenius_coeffs=sol.Arrhenius_coeffs,
             efficiencies_coeffs=sol.efficiencies_coeffs,
             species_viscosities_poly=species_viscosities_poly,
             thermal_conductivity_poly=thermal_conductivity_poly,
             binary_diff_coeffs_poly=binary_diff_coeffs_poly,
             )
else:
    np.savez(mech_yaml,
             molecular_weights=gas.molecular_weights.tolist(),
             reactant_stoich_coeffs=gas.reactant_stoich_coeffs(),
             product_stoich_coeffs=gas.product_stoich_coeffs(),
             reactant_orders=sol.reactant_orders,
             is_reversible=sol.is_reversible,
             Arrhenius_coeffs=sol.Arrhenius_coeffs,
             efficiencies_coeffs=sol.efficiencies_coeffs,
             species_viscosities_poly=species_viscosities_poly,
             thermal_conductivity_poly=thermal_conductivity_poly,
             binary_diff_coeffs_poly=binary_diff_coeffs_poly,
             Arrhenius_A0=sol.Arrhenius_A0,
             Arrhenius_b0=sol.Arrhenius_b0,
             Arrhenius_Ea0=sol.Arrhenius_Ea0,
             Troe_A=sol.Troe_A,
             Troe_T1=sol.Troe_T1,
             Troe_T2=sol.Troe_T2,
             Troe_T3=sol.Troe_T3
             )
