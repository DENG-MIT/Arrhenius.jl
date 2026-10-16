using Arrhenius, Test
include(joinpath(@__DIR__,"..","example","reactors","mixing_solver.jl"))
function mixer_contracts(input)
    gas=CreateSolution(joinpath(input,"gri30.yaml"))
    air=CreateSolution(joinpath(input,"air.yaml"))
    expected=solve_mixing_network(gas,air;save_history=true)
    @testset "Mixer fresh-state and failure contracts" begin
        @test_throws ArgumentError solve_mixing_network(gas,air;save_history=1)
        @test_throws ArgumentError solve_mixing_network(gas,gas)
        @test_throws MethodError solve_mixing_network(gas)
        for value in (0.,-1.,Inf,NaN)
            @test_throws ArgumentError solve_mixing_network(gas,air;residual_tolerance=value)
        end
        @test_throws ArgumentError solve_mixing_network(gas,air;max_iterations=0)
        @test_throws ArgumentError solve_mixing_network(gas,air;max_iterations=1.5)
        @test_throws ArgumentError solve_mixing_network(gas,air;max_backtracks=-1)
        @test_throws ArgumentError solve_mixing_network(gas,air;max_backtracks=1.5)
        stopped=solve_mixing_network(gas,air;max_iterations=1,save_history=true)
        @test !stopped.converged
        @test stopped.iterations==1
        @test stopped.physical_residual>1e-9
        @test size(stopped.states)==(54,2)
        @test stopped.states[:,1]==network_state(stopped.network)
        early=solve_mixing_network(gas,air;residual_tolerance=10.,save_history=true)
        @test early.converged && early.iterations==0
        @test size(early.states)==(54,1)
        @test early.state==network_state(early.network)
        replay=solve_mixing_network(gas,air)
        @test replay.converged
        @test replay.history===nothing && replay.states===nothing
        @test reinterpret(UInt64,replay.state)==reinterpret(UInt64,expected.state)
        @test replay.network!==expected.network
        @test network_state(replay.network)==network_state(expected.network)
        replay.state[1]+=1.
        @test replay.state!=expected.state
        again=solve_mixing_network(gas,air)
        @test reinterpret(UInt64,again.state)==reinterpret(UInt64,expected.state)
    end
end
length(ARGS)==1 || error("supply prepared mixer input directory")
mixer_contracts(ARGS[1])
