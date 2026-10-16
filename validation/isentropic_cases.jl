using Arrhenius, NPZ, Test, TOML
include(joinpath(@__DIR__,"..","example","thermodynamics","isentropic_units.jl"))
include(joinpath(@__DIR__,"..","example","thermodynamics","sound_speed_units.jl"))
root = ARGS[1]
reference = npzread(joinpath(root,"reference.npz"))
gases = Dict(name=>CreateSolution(joinpath(root,name)) for name in ("h2o2.yaml","gri30.yaml","gri30_highT.yaml"))
native = Dict{String,Any}()
report = Dict{String,Any}("reference"=>"Cantera 4.0.0a2; frozen entropy refinement and equilibrium rtol=1e-13",
    "performance_qualified"=>false,"cases"=>Dict{String,Any}())
@testset "Cantera 4 nozzle and acoustic examples" begin
    for (name,calculation,mechanism) in (("isentropic",isentropic_calculation,"h2o2.yaml"),
            ("isentropic_units",isentropic_units_calculation,"gri30.yaml"),
            ("sound_speed",sound_speed_calculation,"gri30_highT.yaml"),
            ("sound_speed_units",sound_speed_units_calculation,"gri30.yaml"))
        result = calculation(gases[mechanism])
        metrics = Dict{String,Any}()
        report["cases"][name] = metrics
        @testset "$name" begin
            for key in keys(result)
                value = result[key] isa Number ? [result[key]] : result[key]
                native[name*"_"*string(key)] = value
                expected = reference[name*"_"*string(key)]
                if startswith(name,"sound_speed")
                    # Pressure finite differences amplify the equilibrium tolerance.
                    # Compare both the published and independently refined references.
                    refined = reference[name*"_refined2_"*string(key)]
                    previous = reference[name*"_refined_"*string(key)]
                    if key == :final_states
                        @test all(isapprox.(value[1:2,:],refined[1:2,:];rtol=2e-8,atol=1e-8))
                        @test all(isapprox.(value[3:end,:],refined[3:end,:];rtol=2e-7,atol=5e-9))
                    else
                        @test all(isapprox.(value,refined;rtol=1e-7,atol=1e-9))
                        @test all(isapprox.(previous,refined;rtol=1e-6,atol=1e-9))
                    end
                    metrics[string(key)] = Dict("max_refined_error"=>maximum(abs.(value-refined)),
                        "max_reference_refinement_change"=>maximum(abs.(previous-refined)),
                        "max_published_reference_error"=>maximum(abs.(expected-refined)),
                        "published_agreement_at_1e_4"=>all(isapprox.(value,expected;rtol=1e-4,atol=1e-9)))
                    println(name," ",key," max refined error ",maximum(abs.(value-refined)),
                            "; published/refined ",maximum(abs.(expected-refined)))
                else
                    @test all(isapprox.(value,expected;rtol=2e-8,atol=1e-8))
                    metrics[string(key)] = Dict("max_reference_error"=>maximum(abs.(value-expected)))
                    println(name," ",key," max error ",maximum(abs.(value-expected)))
                end
            end
        end
    end
    npzwrite(joinpath(root,"native.npz"),native)
    open(joinpath(root,"validation.toml"),"w") do io
        TOML.print(io,report)
    end
end
