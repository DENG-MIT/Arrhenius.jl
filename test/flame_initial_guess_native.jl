using Arrhenius, LinearAlgebra, Test

@testset "flame initial guess from tabulated data" begin
    gas = CreateSolution(joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml"))
    X = mole_fractions(gas, "H2:1.1,O2:1,AR:5")
    f = FreeFlame(gas; X, discretization=:finite_difference)
    solve!(f)
    @test f.converged

    mktempdir() do directory
        csv = joinpath(directory, "flame.csv")
        save_flame(csv, f; basis=:mole, overwrite=true)
        table = read_flame_csv(csv)
        @test table.columns[1:5] ==
              ["z_m", "velocity_m_s", "temperature_K", "density_kg_m3", "heat_release_W_m3"]
        @test size(table.data) == (length(f.grid), length(table.columns))
        @test table.data[:, 1] == f.grid
        @test table.data[:, 2] ≈ velocity(f)

        # A restart from the saved CSV reproduces the flame; only the
        # fixed-temperature point is reselected on the loaded profile. The
        # comparison normalizes each node, as the restart fractions are.
        normalize(Y) = Y ./ sum(Y; dims=1)
        restarted = FreeFlame(gas; X, discretization=:finite_difference)
        set_initial_guess!(restarted, csv)
        @test restarted.grid == f.grid
        @test !restarted.converged
        @test 2 <= restarted.anchor <= length(restarted.grid) - 1
        @test temperature(restarted) ≈ temperature(f)
        @test velocity(restarted) ≈ velocity(f)
        @test normalize(mass_fractions(restarted)) ≈ normalize(mass_fractions(f)) atol=1e-14
        solve!(restarted)
        @test restarted.converged
        @test flame_speed(restarted) ≈ flame_speed(f) rtol=.001

        # Pruned tables restart the same calculation (Cantera's example removes
        # half the grid points, then solves).
        pruned = (; table.columns, data=table.data[1:2:end, :])
        half = FreeFlame(gas; X, discretization=:finite_difference)
        set_initial_guess!(half, pruned)
        @test half.grid ≈ f.grid[1:2:end]
        @test temperature(half) ≈ temperature(f)[1:2:end]
        solve!(half)
        @test half.converged
        @test flame_speed(half) ≈ flame_speed(f) rtol=.01

        # Tables without species columns keep the inlet composition, as a fresh
        # Cantera flame keeps its phase composition on a composition-free
        # restore.
        speciesless = (; columns=table.columns[1:5], data=table.data[1:2:end, 1:5])
        empty = FreeFlame(gas; X, discretization=:finite_difference)
        set_initial_guess!(empty, speciesless)
        @test mass_fractions(empty) ≈ repeat(empty.inlet_Y, 1, length(empty.grid))
        @test velocity(empty) ≈ velocity(f)[1:2:end]
        solve!(empty)
        @test empty.converged
        @test flame_speed(empty) ≈ flame_speed(f) rtol=.01

        # Cantera's CSV layout (grid,velocity,T,D,X_<species>) reads back too.
        ctcsv = joinpath(directory, "cantera.csv")
        Xcols = "X_" .* gas.species_names
        open(ctcsv, "w") do io
            println(io, join(["grid", "velocity", "T", "D", Xcols...], ","))
            data = hcat(f.grid, velocity(f), temperature(f), density(f), transpose(mole_fractions(f)))
            for row in eachrow(data)
                println(io, join(row, ","))
            end
        end
        cantera_restart = FreeFlame(gas; X, discretization=:finite_difference)
        set_initial_guess!(cantera_restart, ctcsv)
        @test cantera_restart.grid == f.grid
        @test temperature(cantera_restart) ≈ temperature(f)
        @test velocity(cantera_restart) ≈ velocity(f)
        @test normalize(mass_fractions(cantera_restart)) ≈ normalize(mass_fractions(f)) atol=1e-14

        # A fresh inlet temperature shifts the loaded profile, matching
        # Cantera's boundary-condition adjustment on restart.
        hot = FreeFlame(gas; X, T=350, discretization=:finite_difference)
        set_initial_guess!(hot, csv)
        @test temperature(hot)[1] == 350
        @test temperature(hot) ≈ temperature(f) .+ 50
        @test velocity(hot) ≈ velocity(f)

        # Error paths.
        @test_throws ArgumentError read_flame_csv(joinpath(directory, "missing.csv"))
        @test_throws ArgumentError set_initial_guess!(restarted, "not a table")
        bad_grid = (; columns=table.columns, data=table.data[:, :])
        bad_grid.data[:, 1] .= f.grid[1]
        @test_throws ArgumentError set_initial_guess!(FreeFlame(gas; X), bad_grid)
        unknown_species = (;
            columns=["z_m", "velocity_m_s", "temperature_K", "X_Unobtainium"],
            data=table.data[:, [1, 2, 3, 6]],
        )
        @test_throws ArgumentError set_initial_guess!(FreeFlame(gas; X), unknown_species)
        no_velocity = (; columns=table.columns[[1, 3]], data=table.data[:, [1, 3]])
        @test_throws ArgumentError set_initial_guess!(FreeFlame(gas; X), no_velocity)
    end

    # Mole-fraction tables renormalize per point and convert through the
    # species molecular weights; missing species stay zero.
    mktempdir() do directory
        tiny = joinpath(directory, "tiny.csv")
        open(tiny, "w") do io
            println(io, "grid,velocity,T,D,X_H2,X_O2")
            for z in range(0, 0.03, 5)
                println(io, join([z, 1.0, 400.0, 1.0, 0.5, 0.25], ","))
            end
        end
        g = FreeFlame(gas; X, discretization=:finite_difference)
        set_initial_guess!(g, tiny)
        iH2 = findfirst(==("H2"), gas.species_names)
        iO2 = findfirst(==("O2"), gas.species_names)
        iAr = findfirst(==("AR"), gas.species_names)
        @test mass_fractions(g)[iAr, :] == zeros(5)
        h2 = 0.5 * gas.MW[iH2] / (0.5 * gas.MW[iH2] + 0.25 * gas.MW[iO2])
        @test mass_fractions(g)[iH2, 1] ≈ h2
        @test sum(mass_fractions(g)[:, 1]) ≈ 1
    end

    # Burner flames retain their prescribed mass flux through a restart.
    burner = BurnerFlame(gas; X="H2:1.5,O2:1,AR:7", T=373., P=.05one_atm,
        mdot=.06, width=.5, discretization=:finite_difference)
    solve!(burner; slope=.05, curve=.1)
    @test burner.converged
    mktempdir() do directory
        csv = joinpath(directory, "burner.csv")
        save_flame(csv, burner; basis=:mole, overwrite=true)
        restarted = BurnerFlame(gas; X="H2:1.5,O2:1,AR:7", T=373., P=.05one_atm,
            mdot=.06, width=.5, discretization=:finite_difference)
        set_initial_guess!(restarted, csv)
        @test restarted.grid == burner.grid
        @test all(==(0.06), restarted.state[end, :])
        @test temperature(restarted) ≈ temperature(burner)
        normalize = Y -> Y ./ sum(Y; dims=1)
        @test normalize(mass_fractions(restarted)) ≈ normalize(mass_fractions(burner)) atol=1e-14
        solve!(restarted; slope=.05, curve=.1)
        @test restarted.converged
        @test flame_speed(restarted) ≈ flame_speed(burner) rtol=.001
    end
end
