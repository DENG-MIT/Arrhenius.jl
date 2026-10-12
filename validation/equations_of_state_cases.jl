# julia --project=EOS_ENV validation/equations_of_state_cases.jl INPUT_DIR OUTPUT.npz
# EOS_ENV needs Arrhenius and Clapeyron 0.6.28.
include(joinpath(@__DIR__,"..","example","thermodynamics","equations_of_state.jl"))
using Arrhenius.NPZ, Libdl, SHA, Test

function main(args)
    length(args) == 2 || error("supply INPUT_DIR OUTPUT.npz")
    input,output = args
    path = joinpath(input,"co2-thermo.yaml")
    ideal = CreateSolution(path)
    rk = RedlichKwongThermo(path;phase="CO2-RK")
    parameter_path = joinpath(input,"carbon-dioxide.json")
    helmholtz = SingleFluid(read(parameter_path,String);coolprop_userlocations=false)
    result = equations_of_state_calculation(ideal,rk,helmholtz)
    arrays = Dict{String,Any}("T"=>[result.T],"pressure"=>result.pressure,
        "saturation_pressure"=>[result.helmholtz.saturation_pressure],
        "helmholtz_phase"=>[p == :vapour ? 1 : 2 for p in result.helmholtz.phase])
    for model in (:ideal,:rk,:helmholtz)
        state = getproperty(result,model)
        arrays[String(model)] = state.properties
        arrays[String(model)*"_density"] = state.density
    end
    arrays["parameter_sha256_utf8"] = collect(codeunits(bytes2hex(sha256(read(parameter_path)))))
    arrays["example_sha256_utf8"] = collect(codeunits(bytes2hex(sha256(read(joinpath(@__DIR__,"..","example","thermodynamics","equations_of_state.jl"))))))
    # Save the full output before any diagnostic can reject it.
    npzwrite(output,arrays)
    @testset "CO2 native phase, identities and replay" begin
        @test size(result.helmholtz.properties) == (5,1000)
        @test Set(result.helmholtz.phase) == Set((:vapour,:liquid))
        for model in (:ideal,:rk,:helmholtz)
            state = getproperty(result,model)
            p = state.properties
            @test all(isfinite,p)
            @test all(>(0),state.density)
            @test all(>(0),p[4:5,:])
            pv = result.pressure ./ state.density ./ 1000
            @test all(isapprox.(p[1,:].-p[2,:],pv.-pv[1];rtol=1e-8,atol=1e-8))
            model == :ideal || @test all(isapprox.(state.pressure,result.pressure;rtol=1e-8,atol=0.))
        end
        repeated = equations_of_state_calculation(ideal,rk,helmholtz)
        @test isequal(result,repeated)
        @test !any(path -> occursin(r"(?i)coolprop|cantera|python",path),Libdl.dllist())
    end
    println("Saved all 15,000 properties and 3,000 densities; performance not measured.")
end
main(ARGS)
