using Arrhenius, Test

function network_allocation_measurements(rhs, u, jacobian!)
    du, J = zero(u), zeros(length(u), length(u))
    rhs(du, u, nothing, 0.003)
    jacobian!(J, u, rhs, 0.003)
    rhs_bytes = @allocated rhs(du, u, nothing, 0.003)
    jacobian_bytes = @allocated jacobian!(J, u, rhs, 0.003)
    return (; rhs_bytes, jacobian_bytes)
end

@testset "allocation-free heterogeneous network callbacks" begin
    gas = CreateSolution(joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml"))
    boundary = Reservoir(gas; temperature=300.0, pressure=2one_atm,
        mole_fractions=Dict("AR"=>1))
    chamber = WellStirredReactor(gas; temperature=900.0, pressure=one_atm,
        mole_fractions=Dict("H2"=>2, "O2"=>1, "AR"=>4), volume=0.03)
    feed = MassFlowController(:inlet, :chamber; mdot=(states,t)->0.01states.chamber.mass)
    exhaust = PressureController(:chamber, :outlet; primary=feed, K=1e-6)
    network = ReactorNetwork((inlet=boundary, chamber=chamber, outlet=boundary);
        flows=(feed, exhaust), walls=(HeatTransferWall(:chamber,:outlet;U=2),))
    rhs = network_rhs(network)
    @test all(rhs.states[i] === rhs.state_vector[i] for i in eachindex(rhs.state_vector))
    fixed = network_allocation_measurements(rhs, network_state(network), network_jacobian!)
    @test fixed.rhs_bytes == 0
    @test fixed.jacobian_bytes == 0

    moving = MovingWallNetwork(network; walls=(
        MovingWall(:chamber,:outlet;velocity=t->0.01sin(t),U=3),
        InertialWall(:inlet,:chamber;mass=0.2,initial_velocity=0.001)))
    mrhs = moving_wall_rhs(moving)
    mobile = network_allocation_measurements(mrhs, moving_wall_state(moving), moving_wall_jacobian!)
    @test mobile.rhs_bytes == 0
    @test mobile.jacobian_bytes == 0
end
