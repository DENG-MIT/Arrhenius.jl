using Arrhenius, NPZ, Test, LinearAlgebra
using SHA: sha256
include(joinpath(@__DIR__,"..","example","thermodynamics","equivalence_ratio.jl"))

"""Compare every native equivalence_ratio_calculation output against the
pinned Cantera 4 thermo/equivalenceRatio.py reference (NPZ + JSON produced by
validation/equivalence_ratio_case.py from the same gri30.yaml)."""
function validate_equivalence_ratio_case(mechanisms, references)
    mechanism = joinpath(mechanisms,"gri30.yaml")
    data = npzread(joinpath(references,"equivalence_ratio-gri30.npz"))
    text(key) = String(vec(UInt8.(data[key])))
    println("Checking against Cantera ",text("cantera_version_utf8"),
            " (commit ",text("source_commit_utf8"),")")
    gas = CreateSolution(mechanism)
    result = equivalence_ratio_calculation(gas)
    mixtures = split(text("mixture_names_utf8"),'\n')
    Xref, Yref, Tref, Pref = data["X"], data["Y"], data["T_K"], data["P_Pa"]
    # Mixtures recorded before the HP equilibration stay at the 300 K default;
    # later mixtures retain the burnt temperature on both sides.
    fresh = ("X_stoich_mole","X_stoich_mass","X_Z055","X_fresh_burnt_case")
    max_composition_error = 0.0
    @testset "Cantera 4 thermo/equivalenceRatio.py" begin
        @test text("cantera_version_utf8")[1] == '4'
        @test bytes2hex(sha256(read(mechanism))) == text("mechanism_sha256_utf8")
        @test gas.species_names == split(text("species_names_utf8"),'\n')
        @test gas.n_species == 53
        for (i,name) in enumerate(mixtures)
            @testset "$name" begin
                Xnative = getproperty(result,Symbol(name))
                if name == "X_burnt"
                    # Native equilibrium tolerances; see validation/equilibrium_cases.jl.
                    @test result.T_burnt_K ≈ Tref[i] rtol=2e-7 atol=2e-5
                    @test result.P_burnt_Pa == one_atm
                    # Cantera's default HP solve drifts by 8.5e-6 Pa here.
                    # Check the reference within its 1e-9 solver tolerance,
                    # while requiring exact prescribed pressure natively.
                    @test result.P_burnt_Pa ≈ Pref[i] rtol=1e-9
                    @test maximum(abs,Xnative-Xref[:,i]) < 1e-7
                    @test maximum(abs,result.Y_burnt-Yref[:,i]) < 1e-7
                else
                    @test all(isapprox.(Xnative,Xref[:,i];rtol=1e-12,atol=1e-14))
                    Ynative = Xnative .* gas.MW ./ dot(Xnative,gas.MW)
                    @test all(isapprox.(Ynative,Yref[:,i];rtol=1e-12,atol=1e-14))
                    Tnative = name in fresh ? result.T_fresh_K : result.T_post_burnt_K
                    @test Tnative ≈ Tref[i] rtol=2e-7 atol=2e-5
                    @test one_atm ≈ Pref[i] rtol=(name in fresh ? 1e-13 : 1e-9)
                end
                max_composition_error = max(max_composition_error,maximum(abs,Xnative-Xref[:,i]))
            end
        end
        for name in split(text("scalar_names_utf8"),'\n')
            @testset "$name" begin
                native = getproperty(result,Symbol(name))
                tolerance = name in ("phi_burnt","Z_burnt") ? 1e-8 : 1e-9
                @test native ≈ data["scalar_"*name] rtol=tolerance atol=1e-12
            end
        end
    end
    println("Maximum composition deviation: ",max_composition_error)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 2 ||
        error("usage: julia validation/equivalence_ratio_case.jl <mechanisms_dir> <reference_dir>")
    validate_equivalence_ratio_case(ARGS[1],ARGS[2])
end
