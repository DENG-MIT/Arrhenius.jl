using Arrhenius, YAML

# Pass an explicit-species Cantera YAML file, for example upstream airNASA9.yaml.
# Standard-state species properties need no reaction or transport sidecar.
yaml = YAML.load_file(ARGS[1])
thermo = IdealGasThermo(yaml)
names = yaml["phases"][1]["species"]
for T in (300.0,1000.0,6000.0,10000.0)
    properties = species_thermo(thermo,T)
    println("T = ",T," K; P = 1 atm")
    for (i,name) in enumerate(names)
        println(name,": cp/R = ",properties.cp_R[i],", h/(RT) = ",properties.h_RT[i],
                ", s/R = ",properties.s_R[i])
    end
end
