using Arrhenius

# Exact temperature grid and triple-point reference convention of the Cantera
# vapor-dome example. Write the table in the customary °C/bar/kJ units.
temperature_C = [WATER_TMIN-273.15;4.;5.;6.;8.;collect(10.:36.);38.;
    collect(40.:5.:95.);collect(100.:10.:290.);collect(300.:20.:360.);
    collect(370.:373.);WATER_TC-273.15]
water = PureWater()
reference = water_state(water;T=WATER_TMIN,Q=0.)
columns = ("T","P","vf","vfg","vg","uf","ufg","ug","hf","hfg","hg","sf","sfg","sg")
output = isempty(ARGS) ? "saturated_steam_T.csv" : ARGS[1]
open(output,"w") do stream
    println(stream,join(columns,','))
    for degC in temperature_C
        T = degC+273.15
        liquid,vapor = water_state(water;T,Q=0.),water_state(water;T,Q=1.)
        uf,ug = (liquid.u-reference.u)/1000,(vapor.u-reference.u)/1000
        hf,hg = (liquid.h-reference.h+reference.P*reference.v)/1000,(vapor.h-reference.h+reference.P*reference.v)/1000
        sf,sg = (liquid.s-reference.s)/1000,(vapor.s-reference.s)/1000
        values = (degC,water_saturation(water,T).P/1e5,liquid.v,vapor.v-liquid.v,vapor.v,
                  uf,ug-uf,ug,hf,hg-hf,hg,sf,sg-sf,sg)
        println(stream,join(values,','))
    end
end
println("Wrote ",length(temperature_C)," saturated-water rows to ",abspath(output))
