using Arrhenius

"Construct the mixing configuration from Cantera's mix1.py."
function mixing_network(gas, air=gas)
    air_composition = Dict("O2"=>0.21, "N2"=>0.78, "AR"=>0.01)
    inlet_air = Reservoir(air; temperature=300.0, mole_fractions=air_composition)
    inlet_fuel = Reservoir(gas; temperature=300.0, mole_fractions=Dict("CH4"=>1))
    outlet = Reservoir(air; temperature=300.0, mole_fractions=air_composition)
    mixer = WellStirredReactor(gas; temperature=300.0, mole_fractions=air_composition, volume=1.0)
    mfc_air = MassFlowController(:inlet_air, :mixer; mdot=inlet_air.initial.density*2.5/0.21)
    mfc_fuel = MassFlowController(:inlet_fuel, :mixer; mdot=inlet_fuel.initial.density)
    valve = Valve(:mixer, :outlet; K=10.0)
    return ReactorNetwork((inlet_air=inlet_air, inlet_fuel=inlet_fuel, mixer=mixer, outlet=outlet);
                          flows=(mfc_air,mfc_fuel,valve))
end

"Construct the combustor.py vessel with mass-dependent residence-time control."
function combustor_network(gas; residence_time=0.1, burned_state=nothing)
    composition = Dict("CH4"=>1.0, "O2"=>4.0, "N2"=>15.04)
    inlet = Reservoir(gas; temperature=300.0, mole_fractions=composition)
    if isnothing(burned_state)
        burned = equilibrate(gas; T=300.0, P=one_atm, X=composition, mode=:HP)
        temperature, Y = burned.T, burned.Y
    else
        temperature, Y = burned_state.temperature, burned_state.mass_fractions
    end
    combustor = WellStirredReactor(gas; temperature, mass_fractions=Y, volume=1.0)
    exhaust = Reservoir(gas; temperature, mass_fractions=Y)
    # A Ref permits continuation to another residence time without rebuilding
    # the mechanism or replacing the initial burned/exhaust thermodynamic state.
    residence = residence_time isa Ref ? residence_time : Ref(Float64(residence_time))
    mfc = MassFlowController(:inlet,:combustor; mdot=(states,t) -> states.combustor.mass/residence[])
    outlet = PressureController(:combustor,:exhaust; primary=mfc,K=0.01)
    return ReactorNetwork((inlet=inlet,combustor=combustor,exhaust=exhaust); flows=(mfc,outlet))
end

"Two closed vessels exchanging mass and wall heat; total mass and energy are conserved."
function closed_pair_network(gas)
    hot = WellStirredReactor(gas; temperature=800.0, pressure=2one_atm,
        mole_fractions=Dict("N2"=>1), volume=0.02, chemistry=false)
    cold = WellStirredReactor(gas; temperature=400.0,
        mole_fractions=Dict("N2"=>1), volume=0.03, chemistry=false)
    mfc = MassFlowController(:hot,:cold; mdot=t -> 0.0005*(1+0.2sin(3t)))
    valve = Valve(:cold,:hot;K=1e-8)
    wall = HeatTransferWall(:hot,:cold;area=0.1,U=10.0,heat_flux=t -> 20sin(t))
    return ReactorNetwork((hot=hot,cold=cold);flows=(mfc,valve),walls=(wall,))
end
