using Arrhenius

function methane_surface_pfr(path)
    m = SurfaceMechanism(path)
    gas = CreateSolution(m.gas_file)
    composition = Dict("CH4"=>1.0,"O2"=>1.5,"AR"=>0.1)
    inlet = IdealGasReactor(gas;temperature=1073.15,pressure=one_atm,mole_fractions=composition)
    area,porosity,velocity = 1e-4,0.3,0.4/60
    mdot = velocity*inlet.density*area*porosity
    return SurfaceFlowReactor(m;gas,temperature=1073.15,pressure=one_atm,
        mole_fractions=composition,area,mass_flow_rate=mdot,
        surface_area_per_length=1e5*porosity*area)
end
