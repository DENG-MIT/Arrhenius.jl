using Arrhenius, Test, NPZ, SHA

@testset "inert ideal-gas mechanism" begin
    mktempdir() do directory
        path = joinpath(directory,"argon.yaml")
        write(path,"""
        phases:
        - name: argon
          thermo: ideal-gas
          elements: [Ar]
          species: [AR]
        species:
        - name: AR
          composition: {Ar: 1}
          thermo:
            model: NASA7
            temperature-ranges: [200., 6000.]
            data: [[2.5, 0., 0., 0., 0., -745.375, 4.366]]
        """)
        # Minimal numerical sidecar emitted for a phase without reactions or transport.
        npzwrite(path*".npz",Dict("molecular_weights"=>[39.95],
            "sidecar_format_utf8"=>collect(codeunits("arrhenius-sidecar-v2")),
            "source_sha256_utf8"=>collect(codeunits(bytes2hex(sha256(read(path)))))))
        gas = CreateSolution(path)
        @test gas.n_reactions == 0
        @test size(gas.reaction.vk) == (1,0)
        @test set_states(gas,1000.,one_atm,[1.]) == [0.]
        reactor = IdealGasReactor(gas;temperature=1000.,mole_fractions=Dict("AR"=>1.))
        rhs = reactor_rhs(reactor)
        state = reactor_state(reactor)
        derivative = similar(state)
        rhs(derivative,state,nothing,0.)
        @test derivative == [0.,0.]
    end
end
