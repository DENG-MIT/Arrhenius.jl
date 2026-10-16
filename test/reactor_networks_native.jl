using Arrhenius
using LinearAlgebra
using Test

@testset "native reactor network balances and devices" begin
    gas = CreateSolution(joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml"))
    gri = CreateSolution(joinpath(@__DIR__, "..", "mechanism", "gri30.yaml"))
    mixture = Dict("H2" => 2.0, "O2" => 1.0, "AR" => 4.0)
    closed = IdealGasReactor(gas; temperature=1400.0, mole_fractions=mixture,
                             constraint=:constant_volume)
    vessel = WellStirredReactor(closed; volume=0.2)
    network = ReactorNetwork((vessel=vessel,))
    u = network_state(network)
    rhs = network_rhs(network)
    derivative = similar(u)
    rhs(derivative, u, nothing, 0.0)
    closed_derivative = similar(reactor_state(closed))
    reactor_rhs(closed)(closed_derivative, reactor_state(closed), nothing, 0.0)
    mass = closed.density * vessel.volume
    @test derivative[1:end-1] ≈ mass .* closed_derivative[1:end-1] rtol=1e-12
    @test derivative[end] ≈ closed_derivative[end] rtol=1e-12
    balances = network_diagnostics(rhs, u)
    @test abs(balances.total_mass_rate) < 1e-12 * norm(derivative[1:end-1], 1)
    @test abs(balances.total_energy_rate) < 1e-7
    @test balances.external_mass_rate == 0
    @test balances.external_energy_rate == 0

    Y = collect(1.0:gas.n_species)
    left = WellStirredReactor(gas; temperature=1200.0, mass_fractions=Y,
                              volume=0.03, chemistry=false)
    right = WellStirredReactor(gas; temperature=600.0, mass_fractions=reverse(Y),
                               volume=0.05, chemistry=false)
    forward = MassFlowController(:left, :right; mdot=0.005)
    backward = MassFlowController(:right, :left; mdot=(states,t) -> states.right.mass / 4)
    wall = HeatTransferWall(:left, :right; area=0.1, U=25.0, heat_flux=t -> 2t)
    connected = ReactorNetwork((left=left, right=right); flows=(forward, backward), walls=(wall,))
    connected_rhs = network_rhs(connected)
    connected_state = network_state(connected)
    balance = network_diagnostics(connected_rhs, connected_state, 0.3)
    @test abs(balance.total_mass_rate) < 1e-16
    @test abs(balance.total_energy_rate) < 1e-9
    @test maximum(abs, values(balance.element_rates)) < 1e-16
    @test balance.external_energy_rate ≈ 0 atol=1e-9
    @test balance.wall_heat_rates[1] ≈ 0.1 * (25 * 600 + 0.6)
    @test balance.nodes.left.mass ≈ left.initial.density * left.volume
    @test balance.nodes.right.pressure ≈ one_atm

    J = zeros(length(connected_state), length(connected_state))
    before = copy(connected_state)
    network_jacobian!(J, connected_state, connected_rhs, 0.3)
    direction = connected_state .* collect(range(-0.1, 0.2; length=length(connected_state)))
    plus, minus = similar(connected_state), similar(connected_state)
    step = 1e-5
    connected_rhs(plus, connected_state + step * direction, nothing, 0.3)
    connected_rhs(minus, connected_state - step * direction, nothing, 0.3)
    @test J * direction ≈ (plus - minus) / (2step) rtol=1e-6
    @test connected_state == before
    @test !network_isoutofdomain(connected, connected_state)
    invalid = copy(connected_state)
    invalid[1] = -1.0
    @test network_isoutofdomain(connected, invalid)

    source = Reservoir(gas; temperature=900.0, pressure=2one_atm, mole_fractions=Dict("AR"=>1))
    sink = Reservoir(gas; temperature=300.0, pressure=one_atm, mole_fractions=Dict("AR"=>1))
    chamber = WellStirredReactor(gas; temperature=600.0, pressure=1.5one_atm,
        mole_fractions=Dict("AR"=>1), volume=0.01, chemistry=false)
    inlet = MassFlowController(:source, :chamber; mdot=t -> 0.01 * (1+t))
    outlet = PressureController(:chamber, :sink; primary=inlet, K=1e-6)
    pressure_valve = Valve(:source, :chamber; K=1e-6, time_function=t -> 2.0)
    reverse_valve = Valve(:sink, :chamber; K=1e-6)
    negative = MassFlowController(:source, :chamber; mdot=-0.1)
    flowing = ReactorNetwork((source=source, chamber=chamber, sink=sink);
        flows=(inlet, outlet, pressure_valve, reverse_valve, negative),
        walls=(HeatTransferWall(:chamber, :sink; area=2, U=3),))
    flow_balance = network_diagnostics(flowing, network_state(flowing), 2.0)
    @test flow_balance.mass_flow_rates ≈ [0.03, 0.03+0.5one_atm*1e-6, one_atm*1e-6, 0, 0]
    expected_mass = 0.03 + one_atm*1e-6 - (0.03 + 0.5one_atm*1e-6)
    @test flow_balance.total_mass_rate ≈ expected_mass
    @test flow_balance.external_mass_rate ≈ expected_mass
    expected_energy = (0.03 + one_atm*1e-6) * flow_balance.nodes.source.enthalpy -
        flow_balance.mass_flow_rates[2] * flow_balance.nodes.chamber.enthalpy - 2*3*(600-300)
    @test flow_balance.total_energy_rate ≈ expected_energy rtol=1e-12
    @test flow_balance.external_energy_rate ≈ expected_energy rtol=1e-12
    problem = network_problem(flowing, (0.0, 1.0))
    gradient = similar(problem.u0)
    problem.tgrad(gradient, problem.u0, nothing, 2.0)
    @test abs(sum(gradient[1:end-1])) < 1e-8 # primary pressure controller tracks the inlet
    @test length(problem.u0) == gas.n_species + 1

    bath = WellStirredReactor(gas; temperature=600.0, mole_fractions=Dict("AR"=>1),
        chemistry=false, energy=:isothermal)
    thermostated = ReactorNetwork((bath=bath, cold=sink);
        walls=(HeatTransferWall(:bath, :cold; U=1.0),))
    thermostat = network_diagnostics(thermostated)
    @test thermostat.thermostat_heat_rates[1] ≈ 300.0
    @test thermostat.total_energy_rate == 0
    @test thermostat.external_energy_rate == 0

    mixed_receiver = WellStirredReactor(gri; temperature=300.0,
        mole_fractions=Dict("O2"=>0.21,"N2"=>0.78,"AR"=>0.01), chemistry=false)
    cross_mechanism = ReactorNetwork((source=source, target=mixed_receiver);
        flows=(MassFlowController(:source,:target; mdot=0.01),))
    @test network_diagnostics(cross_mechanism).total_mass_rate ≈ 0.01
    methane = Reservoir(gri; temperature=300.0, mole_fractions=Dict("CH4"=>1))
    @test_throws ArgumentError ReactorNetwork((methane=methane, chamber=chamber);
        flows=(MassFlowController(:methane,:chamber;mdot=0.01),))
    @test_throws ArgumentError ReactorNetwork((source=source,))
    @test_throws ArgumentError ReactorNetwork((chamber=chamber,);
        flows=(MassFlowController(:missing,:chamber;mdot=0.01),))
    @test_throws ArgumentError ReactorNetwork((source=source,chamber=chamber,sink=sink);
        flows=(outlet,inlet))
    @test_throws ArgumentError Valve(:source,:chamber;K=-1)
    @test_throws ArgumentError HeatTransferWall(:source,:chamber;U=-1)
    @test_throws ArgumentError WellStirredReactor(closed;volume=0)
    @test_throws ArgumentError WellStirredReactor(IdealGasReactor(gas;
        temperature=300.0,mole_fractions=Dict("AR"=>1)))
    @test_throws ArgumentError network_problem(network,(1.0,0.0))
end
