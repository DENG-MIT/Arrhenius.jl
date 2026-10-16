using Arrhenius, Test, LinearAlgebra

@testset "temperature cache preserves concentration-dependent rates" begin
    for name in ("h2o2","gri30")
        gas = CreateSolution(joinpath(@__DIR__,"..","mechanism",name*".yaml"))
        cache = Arrhenius._KineticsTemperatureCache(gas.reaction)
        w = KineticsWorkspace(gas.reaction)
        uncached = similar(w.rates_of_progress)
        cached = similar(uncached)
        production = zeros(gas.n_species)
        # Repeated temperatures with changing P/X exercise cached falloff and
        # third-body factors; revisiting temperatures exercises invalidation.
        for T in (700.,1300.,1300.,700.,2400.), P in (one_atm,10one_atm), shift in (0,1)
            X = Float64.(circshift(collect(1:gas.n_species),shift))
            X ./= sum(X)
            C = P/(R*T) .* X
            H = R*T .* cal_h_RT(gas,T,P,X)
            S = R .* cal_s0_R(gas,T,P,X)
            uncached .= wdot!(production,gas.reaction,T,C,S,H,w;get_qdot=true)
            cached .= wdot!(production,gas.reaction,T,C,S,H,w;get_qdot=true,temperature_cache=cache)
            @test cached ≈ uncached rtol=1e-13 atol=1e-20
        end
    end
end
