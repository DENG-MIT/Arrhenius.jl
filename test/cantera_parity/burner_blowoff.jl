using Arrhenius, SHA, Test

@testset "Cantera BurnerFlame blowoff" begin
    mechanism = joinpath(@__DIR__, "..", "data", "cantera_flux_gradient_mass", "h2o2.yaml")
    @test bytes2hex(sha256(read(mechanism))) ==
          "0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4"
    @test bytes2hex(sha256(read(mechanism * ".npz"))) ==
          "1cc087ea45a8451ebf575ee35d9edee90445c6d08260f0a92669ee4457f528b5"

    gas = CreateSolution(mechanism)
    # Cantera's equivalence ratio 0.4 for H2 with O2:1, AR:5.
    flame = BurnerFlame(gas; T=300.0, P=one_atm, X="H2:0.8, O2:1, AR:5",
        width=0.1, mdot=1.2)
    solve!(flame; auto=true, ratio=3.0, slope=0.3, curve=0.5, prune=0.0)
    @test flame.converged
    @test length(flame.grid) > 1 && all(>(0), diff(flame.grid))

    # Cantera 4.0.0a2 (726522be4e2a13454d8415b7ef799d621f665cf3),
    # TestBurnerFlame.test_blowoff on the same YAML. Species order:
    # H2, H, O, O2, OH, H2O, HO2, H2O2, AR, N2.
    reference = (
        inlet=(T=300.0, u=0.8607976668154603, Y=[
            0.006911186454623055, 2.846587401482386e-26, 4.562091760132423e-25,
            0.1371181449497945, -4.269074392437515e-24, -8.130158740962498e-26,
            7.416266840264795e-24, 1.929251266477797e-24, 0.8559706685955825,
            -3.571091511397917e-36,
        ]),
        outlet=(T=300.00007901413596, u=0.8607969642651884, Y=[
            0.006911191401055586, 1.867506442117825e-18, 1.5203125381116922e-10,
            0.13711812519029973, -2.0892713083551982e-13, 1.7736166765373135e-08,
            -6.411618793980614e-09, 2.8558833665499445e-10, 0.8559706716466868,
            -2.3419287293824742e-20,
        ]),
    )
    t, u, y = temperature(flame), velocity(flame), mass_fractions(flame)
    @test all(isfinite, t) && all(isfinite, u) && all(isfinite, y)
    native = (
        inlet=(T=t[1], u=u[1], Y=y[:, 1]),
        outlet=(T=t[end], u=u[end], Y=y[:, end]),
    )

    for (name, endpoints) in (("Cantera reference", reference), ("Arrhenius", native))
        @testset "$name" begin
            inlet, outlet = endpoints.inlet, endpoints.outlet
            # pytest.approx scales by the expected endpoint. T and u use the
            # inlet and the default absolute floor; Y uses the outlet per species.
            @test abs(outlet.T - inlet.T) <= max(1e-6 * abs(inlet.T), 1e-12)
            @test abs(outlet.u - inlet.u) <= max(1e-6 * abs(inlet.u), 1e-12)
            @test all(abs.(inlet.Y .- outlet.Y) .<= max.(1e-6 .* abs.(outlet.Y), 1e-6))
        end
    end
end
