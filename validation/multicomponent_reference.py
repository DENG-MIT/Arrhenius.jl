"""Generate independent Cantera coefficient and native C++ flux references."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import cantera as ct
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mechanism"))
from export_multicomponent import export_multicomponent
from export_sidecar import export as export_mechanism

CPP = r'''
#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/transport/Transport.h"
#include <iostream>
#include <iomanip>
#include <vector>
int main(int argc, char** argv) {
    auto sol = Cantera::newSolution(argv[1], "", "multicomponent");
    auto gas = sol->thermo();
    auto tr = sol->transport();
    size_t n = gas->nSpecies();
    std::vector<double> x(n), gx(n), f(n);
    double T, P, gT;
    std::cout << std::setprecision(17);
    while (std::cin >> T >> P >> gT) {
        for (auto& a : x) std::cin >> a;
        for (auto& a : gx) std::cin >> a;
        gas->setState_TPX(T, P, x);
        // Refresh composition before getSpeciesFluxes, as required by this API.
        tr->thermalConductivity();
        tr->getSpeciesFluxes(1, std::span<const double>(&gT,1),n,gx,n,f);
        for (auto a : f) std::cout << a << ' ';
        std::cout << '\n';
    }
}
'''


def generate(output: Path, prefix: Path):
    output.mkdir(parents=True, exist_ok=True)
    source = output / "multicomponent_flux_reference.cpp"
    binary = output / "multicomponent_flux_reference"
    source.write_text(CPP)
    command = ["/usr/bin/clang++", "-std=c++20", "-O2", str(source), "-o", str(binary),
               "-I"+str(prefix/"include"), "-L"+str(prefix/"lib"),
               "-lcantera_shared", "-Wl,-rpath,"+str(prefix/"lib")]
    if sys.platform == "darwin":
        command += ["-isystem", "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1"]
    subprocess.run(command, check=True)
    summary = {"cantera_version":ct.__version__, "flux_reference":"C++ Transport::getSpeciesFluxes",
               "source_commit":"726522be4e2a13454d8415b7ef799d621f665cf3", "cases":{}}
    for name, fuel in [("h2o2", "H2"), ("gri30", "CH4")]:
        mechanism = str(prefix/"share/cantera/data"/(name+".yaml"))
        export_multicomponent(mechanism, output/(name+"-transport.npz"))
        if name == "h2o2":
            native_mechanism=output/(name+".yaml")
            shutil.copyfile(mechanism,native_mechanism)
            export_mechanism(native_mechanism,output/(name+".yaml.npz"))
        gas = ct.Solution(mechanism, transport_model="multicomponent")
        mixture = ct.Solution(mechanism, transport_model="mixture-averaged")
        compositions = []
        gas.set_equivalence_ratio(1.,fuel,"O2:1,N2:3.76")
        compositions.append(gas.X.copy())
        gas.TP = 300, ct.one_atm
        gas.equilibrate("HP")
        compositions.append(gas.X.copy())
        gas.X = "H2:0.03, AR:0.97"
        compositions.append(gas.X.copy())
        gas.X = "H2:2, O2:1, N2:3.76, AR:5"
        compositions.append(gas.X.copy())
        gas.X = "AR:1"
        compositions.append(gas.X.copy())
        rng = np.random.default_rng(413)
        gas.X = rng.uniform(.1,1,gas.n_species)
        compositions.append(gas.X.copy())
        rows=[]
        for composition_index, X in enumerate(compositions):
            for T in [300.,700.,1500.,2500.]:
                for P in [.2*ct.one_atm,ct.one_atm,10*ct.one_atm]:
                    gas.TPX=T,P,X
                    mixture.TPX=T,P,X
                    # A physical local gradient proportional to X also covers exact zeros.
                    gx=X*rng.normal(size=gas.n_species)
                    gx-=X*sum(gx)
                    gT=321.
                    rows.append(dict(T=T,P=P,X=gas.X.copy(),cp_R=gas.standard_cp_R.copy(),
                        diffusion=gas.multi_diff_coeffs.copy(),thermal_diffusion=gas.thermal_diff_coeffs.copy(),
                        conductivity=gas.thermal_conductivity,grad_X=gx,grad_T=gT,
                        mixture_thermal_diffusion=mixture.thermal_diff_coeffs.copy(),
                        mixture_diffusion=mixture.mix_diff_coeffs.copy(),
                        composition_index=composition_index))
        stdin="\n".join(" ".join(map(str,[r['T'],r['P'],r['grad_T'],*r['X'],*r['grad_X']])) for r in rows)+"\n"
        env=os.environ.copy()
        env["CANTERA_DATA"]=str(prefix/"share/cantera/data")
        result=subprocess.run([str(binary),mechanism],input=stdin,text=True,capture_output=True,check=True,env=env)
        fluxes=np.array([np.fromstring(line,sep=" ") for line in result.stdout.splitlines()])
        assert fluxes.shape==(len(rows),gas.n_species), result.stdout
        data={key:np.array([row[key] for row in rows]) for key in rows[0]}
        data['flux']=fluxes
        np.savez_compressed(output/(name+"-reference.npz"),**data)
        summary['cases'][name]=len(rows)
    (output/'reference-summary.json').write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--cantera-prefix",type=Path,required=True)
    args=parser.parse_args()
    generate(args.output,args.cantera_prefix)
