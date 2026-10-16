using Arrhenius, NPZ, Test
isdefined(Arrhenius,:PureWater) || Base.include(Arrhenius,joinpath(@__DIR__,"..","src","PureWater.jl"))

water = PureWater()
references = ARGS[1]
keys = (:T,:P,:Q,:rho,:h,:u,:s,:cp,:cv)
function compare(state,reference;heat_capacity=true)
    for (i,key) in enumerate(keys)
        !heat_capacity && key in (:cp,:cv) && continue
        expected = reference[i]
        isnan(expected) && continue
        value = getproperty(state,key)
        tolerance = key in (:cp,:cv) ? 2e-4 : key == :Q ? 2e-7 : 2e-7
        @test isapprox(value,expected;rtol=tolerance,atol=key in (:h,:u) ? .05 : key == :s ? 1e-4 : 2e-6)
    end
end

@testset "Cantera water vapor dome" begin
    data = npzread(joinpath(references,"water-dome.npz"))
    for i in eachindex(data["T"]), (j,Q) in enumerate((0.,1.))
        @testset "T=$(data["T"][i]), Q=$Q" begin
            state = water_state(water;T=data["T"][i],Q)
            compare(state,data["states"][:,j,i])
        end
    end
end
@testset "Cantera water TP/TQ/PQ states" begin
    data = npzread(joinpath(references,"water-states.npz"))
    for i in axes(data["inputs"],2)
        mode,a,b = data["inputs"][:,i]
        @testset "case $i/$mode/$a/$b" begin
            state = mode == 1 ? water_state(water;T=a,P=b) : mode == 2 ? water_state(water;T=a,Q=b) : water_state(water;P=a,Q=b)
            compare(state,data["states"][:,i])
        end
    end
end
@testset "Cantera Rankine cycles" begin
    data = npzread(joinpath(references,"water-rankine.npz"))
    for (i,(initialT,maximumP)) in enumerate(((300.,8e5),((80.33-32)*5/9+273.15,116.03*6894.757293168364)))
        first = water_state(water;T=initialT,Q=0.)
        pumpideal = water_state(water;s=first.s,P=maximumP)
        pump_work = (pumpideal.h-first.h)/.6
        second = water_state(water;h=first.h+pump_work,P=maximumP)
        third = water_state(water;P=maximumP,Q=1.)
        expansionideal = water_state(water;s=third.s,P=first.P)
        turbine_work = (third.h-expansionideal.h)*.8
        fourth = water_state(water;h=third.h-turbine_work,P=first.P)
        for (j,state) in enumerate((first,pumpideal,second,third,expansionideal,fourth))
            compare(state,data["states"][:,j,i])
        end
        heat = third.h-second.h
        values = [pump_work,turbine_work,heat,(turbine_work-pump_work)/heat]
        @test values ≈ data["metrics"][:,i] rtol=2e-6 atol=1e-3
        println("Rankine $i: ",values)
    end
end
