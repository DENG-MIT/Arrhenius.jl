using Test, Arrhenius

@testset "planar flame refinement interpolates velocity like Cantera" begin
    ref = Arrhenius.npzread(joinpath(@__DIR__, "data", "cantera_flame_refine_velocity.npz"))
    gaspath = joinpath(pkgdir(Arrhenius), "mechanism", "h2o2.yaml")
    cases = (
        ("free", FreeFlame(CreateSolution(gaspath); X="H2:1.6,O2:1,AR:7", width=0.03,
            discretization=:finite_difference)),
        ("burner", BurnerFlame(CreateSolution(gaspath); mdot=0.04, X="H2:1.6,O2:1,AR:7",
            width=0.03, discretization=:finite_difference)),
        ("fixedT", BurnerFlame(CreateSolution(gaspath); mdot=0.04, X="H2:1.6,O2:1,AR:7",
            width=0.03, discretization=:finite_difference)),
    )
    for (kind, f) in cases
        kind == "fixedT" && set_temperature_profile!(f, [0., .25, .6, 1.], [300., 420., 670., 900.])
        pre = "$(kind)/pre/"
        f.grid = Vector{Float64}(ref[pre*"grid"])
        Y = Matrix{Float64}(ref[pre*"Y"])
        T = Vector{Float64}(ref[pre*"T"])
        u = Vector{Float64}(ref[pre*"velocity"])
        f.state[1,:] .= T ./ 1000
        f.state[2:end-1,:] .= Y
        f.state[end,:] .= [Arrhenius._flame_density(f,T[j],@view(Y[:,j])).density*u[j] for j in eachindex(T)]
        @test temperature(f) == T
        @test mass_fractions(f) == Y
        @test velocity(f) ≈ u rtol=2e-14 atol=2e-14
        @test maximum(abs.(sum(Y;dims=1) .- 1)) < 1e-12
        if kind == "free"
            f.anchor = findfirst(==(Float64(ref[pre*"fixed_point"])), f.grid)
            f.fixed_temperature = Float64(ref[pre*"fixed_temperature"])
        end
        old_anchor_z = f.grid[f.anchor]
        @test Arrhenius._refine_flame!(f; ratio=3., slope=.03, curve=.06, prune=0.)
        post = "$(kind)/post/"
        @test f.grid == Vector{Float64}(ref[post*"grid"])
        @test temperature(f) ≈ Vector{Float64}(ref[post*"T"]) rtol=1e-12 atol=1e-12
        @test mass_fractions(f) ≈ Matrix{Float64}(ref[post*"Y"]) rtol=1e-12 atol=1e-14
        @test velocity(f) ≈ Vector{Float64}(ref[post*"velocity"]) rtol=1e-12 atol=1e-12
        @test f.grid[f.anchor] == old_anchor_z
        for name in ("spreadRate", "Lambda", "eField", "Uo")
            @test maximum(abs.(Vector{Float64}(ref[post*"component/"*name]))) == 0
        end
        kind == "burner" && @test f.mass_flux == 0.04
        @test f.state[end,1] ≈ Arrhenius._flame_density(f,temperature(f)[1],@view(mass_fractions(f)[:,1])).density * velocity(f)[1] rtol=2e-14
        if kind == "fixedT"
            expected_profile = [300., 420., 670., 900.]
            expected_positions = [0., .25, .6, 1.] .* (f.grid[end]-f.grid[1]) .+ f.grid[1]
            expected = [Arrhenius._interpolate_profile(expected_positions,expected_profile,z) for z in f.grid]
            @test f.imposed_temperature ≈ expected rtol=1e-14 atol=1e-12
        end

        nr = flame_residual!(similar(f.state), f)
        for name in ("velocity", "T", f.gas.species_names...)
            native_row = name == "velocity" ? size(nr,1) : name == "T" ? 1 : findfirst(==(name), f.gas.species_names)+1
            cr = Vector{Float64}(ref[post*"residual/"*name])
            selected = [i for i in 2:length(cr)-1 if abs(cr[i]) > 1e-3*maximum(abs.(cr))]
            if isempty(selected)
                println("RESIDUAL kind=$kind component=$name no-interior-signal")
            else
                ratios = nr[native_row,selected] ./ cr[selected]
                med = sort(ratios)[cld(length(ratios),2)]
                spread = maximum(abs.(ratios ./ med .- 1))
                println("RESIDUAL kind=$kind component=$name ctmax=$(maximum(abs.(cr))) natmax=$(maximum(abs.(nr[native_row,:]))) ratio=$med spread=$spread")
            end
        end
    end
end