using Arrhenius
using LinearAlgebra
using Test

@testset "native variable-volume and pressure-work balances" begin
    gas = CreateSolution(joinpath(@__DIR__,"..","mechanism","h2o2.yaml"))
    left = WellStirredReactor(gas;temperature=900.0,pressure=2one_atm,
        mole_fractions=Dict("AR"=>1.0),volume=0.3,chemistry=false)
    right = WellStirredReactor(gas;temperature=600.0,pressure=one_atm,
        mole_fractions=Dict("AR"=>1.0),volume=0.2,chemistry=false)
    base = ReactorNetwork((left=left,right=right))
    ubase = network_state(base)
    fixed_rhs = network_rhs(base)
    ddefault,dvolumes = similar(ubase),similar(ubase)
    fixed_rhs(ddefault,ubase,nothing,0.0)
    fixed_rhs(dvolumes,ubase,nothing,0.0;volumes=[0.3,0.2])
    @test ddefault == dvolumes
    # The optional volume affects density and rate scaling; default resets cached volume.
    fixed_rhs(dvolumes,ubase,nothing,0.0;volumes=[0.6,0.4])
    @test fixed_rhs.states.left.pressure ≈ one_atm
    @test fixed_rhs.states.left.volume == 0.6
    fixed_rhs(ddefault,ubase,nothing,0.0)
    @test fixed_rhs.states.left.pressure ≈ 2one_atm
    @test fixed_rhs.states.left.volume == 0.3
    @test_throws DimensionMismatch fixed_rhs(dvolumes,ubase,nothing,0.0;volumes=[0.3])
    @test_throws DomainError fixed_rhs(dvolumes,ubase,nothing,0.0;volumes=[0.0,0.2])

    reacting = WellStirredReactor(gas;temperature=1400.0,pressure=2one_atm,volume=0.2,
        mole_fractions=Dict("H2"=>2.0,"O2"=>1.0,"AR"=>4.0))
    expanded = WellStirredReactor(gas;temperature=1400.0,pressure=one_atm,volume=0.4,
        mole_fractions=Dict("H2"=>2.0,"O2"=>1.0,"AR"=>4.0))
    reacting_network = ReactorNetwork((chamber=reacting,))
    expanded_network = ReactorNetwork((chamber=expanded,))
    chemical_state = network_state(reacting_network)
    chemical_actual,chemical_expected = similar(chemical_state),similar(chemical_state)
    network_rhs(reacting_network)(chemical_actual,chemical_state,nothing,0.0;volumes=[0.4])
    network_rhs(expanded_network)(chemical_expected,network_state(expanded_network),nothing,0.0)
    @test chemical_state ≈ network_state(expanded_network) rtol=1e-15
    @test chemical_actual ≈ chemical_expected rtol=1e-12

    stationary = MovingWallNetwork(base;walls=())
    du = similar(moving_wall_state(stationary))
    moving_wall_rhs(stationary)(du,moving_wall_state(stationary),nothing,0.0)
    @test du[1:length(ubase)] == ddefault
    @test all(iszero,du[length(ubase)+1:end])

    wall = MovingWall(:left,:right;area=0.2,K=1e-6,velocity=t->0.01t,U=20.0,heat_flux=3.0)
    model = MovingWallNetwork(base;walls=(wall,))
    u = moving_wall_state(model)
    rhs = moving_wall_rhs(model)
    d = moving_wall_diagnostics(rhs,u,0.5)
    v = 1e-6one_atm+0.005
    @test d.wall_velocities[1] ≈ v
    @test d.volume_rates ≈ [0.2v,-0.2v]
    @test d.wall_heat_rates[1] ≈ 0.2*(20*300+3)
    @test sum(d.volume_rates) == 0
    @test d.total_internal_energy_rate ≈ -one_atm*0.2v rtol=1e-12
    @test d.pressure_work_rate ≈ one_atm*0.2v
    @test d.energy_input_rate == 0
    @test d.external_mass_rate == 0
    @test abs(d.total_internal_energy_rate+d.pressure_work_rate-d.energy_input_rate) < 1e-9
    @test !moving_wall_isoutofdomain(model,u)
    invalid = copy(u)
    invalid[model.volume_indices[1]] = 0
    @test moving_wall_isoutofdomain(model,invalid)
    @test_throws ArgumentError moving_wall_problem(model,(0.0,1.0);initial_state=invalid)
    @test_throws ArgumentError MovingWall(:left,:right;K=-1)
    @test_throws ArgumentError MovingWallNetwork(base;walls=(MovingWall(:left,:left),))
    @test_throws ArgumentError InertialWall(:left,:right;mass=0)

    J = zeros(length(u),length(u))
    moving_wall_jacobian!(J,u,rhs,0.5)
    direction = u.*range(-0.1,0.1;length=length(u))
    plus,minus = similar(u),similar(u)
    rhs(plus,u+1e-5direction,nothing,0.5)
    rhs(minus,u-1e-5direction,nothing,0.5)
    @test J*direction ≈ (plus-minus)/(2e-5) rtol=1e-6
    @test u == moving_wall_state(model)
    @test all(iszero,@view J[:,model.ledger_offset:end])

    inertia = MovingWallNetwork(base;walls=(InertialWall(:left,:right;area=0.2,mass=2.0,initial_velocity=-0.3),))
    ui = moving_wall_state(inertia)
    di = similar(ui)
    irhs = moving_wall_rhs(inertia)
    irhs(di,ui,nothing,0.0)
    @test di[inertia.velocity_indices[1]] ≈ 0.1one_atm
    balance = moving_wall_diagnostics(irhs,ui)
    @test balance.wall_kinetic_energy ≈ 0.09
    @test balance.total_internal_energy_rate + 2*(-0.3)*di[inertia.velocity_indices[1]] ≈ 0 atol=1e-10

    iso = WellStirredReactor(gas;temperature=600.0,pressure=one_atm,energy=:isothermal,
        mole_fractions=Dict("AR"=>1.0),chemistry=false)
    reservoir = Reservoir(gas;temperature=300.0,pressure=one_atm,mole_fractions=Dict("AR"=>1.0))
    open_base = ReactorNetwork((chamber=iso,boundary=reservoir);
        flows=(MassFlowController(:boundary,:chamber;mdot=0.01),))
    iso_model = MovingWallNetwork(open_base;walls=(MovingWall(:chamber,:boundary;velocity=0.1,U=2),))
    iso_d = moving_wall_diagnostics(iso_model)
    @test iso_d.external_mass_rate == 0.01
    @test iso_d.volume_rates == [0.1,0.0]
    @test iso_d.total_internal_energy_rate ≈ iso_d.energy_input_rate-iso_d.pressure_work_rate rtol=1e-12
    @test iso_d.thermostat_heat_rates[1] > 0
end
