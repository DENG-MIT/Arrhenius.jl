# Shared calculation for the three published premixed-flame examples.
# Prepare mechanism data before calling; requested numerical profiles are timed, file output is not.
using Arrhenius

"Temperature profile from the original flame_fixed_T.py example, in metres and kelvin."
function source_flame_temperature_profile()
    positions=[
        0.,.00015625,.00023437,.00039063,.00046875,.00050781,
        .00054688,.000625,.00066406,.00070312,.00074219,.00078125,
        .00082031,.00085938,.00089844,.0009375,.00101563,.00105469,
        .00109375,.00113281,.00117187,.00121094,.00125,.00128906,
        .00132813,.00136719,.00140625,.00144531,.00148438,.00152344,
        .0015625,.00160156,.00164062,.00171875,.00175781,.00179688,
        .00183594,.001875,.00191406,.00195312,.00199219,.00203125,
        .00207031,.00210938,.00214844,.0021875,.00222656,.00226562,
        .00230469,.00234375,.00238281,.00242187,.00246094,.0025,
        .00257813,.00265625,.00273437,.0028125,.00289062,.00296875,
        .00304688,.003125,.00328125,.0034375,.00359375,.00375,
        .00390625,.0087,.01]
    temperatures=[
        373.7,465.4070428,510.4311676,599.5552837,643.8342938,
        665.9335545,688.0122338,732.1284327,754.1744755,776.2170662,
        798.2588757,820.3020011,842.348001,864.3979228,886.4523159,
        908.5112198,952.6396629,974.7018199,996.7515831,1018.777651,
        1040.765863,1062.69948,1084.558639,1106.320078,1127.956918,
        1149.438472,1170.730129,1191.793309,1212.585506,1233.060477,
        1253.168589,1272.857384,1292.072391,1328.859767,1346.323998,
        1363.101361,1379.147594,1394.425274,1408.905834,1422.569115,
        1435.40408,1447.410648,1458.597668,1468.982722,1478.590978,
        1487.453914,1495.607879,1503.092709,1509.950449,1516.224147,
        1521.956853,1527.19079,1531.966722,1536.32348,1543.891739,
        1550.203579,1555.480771,1559.908135,1563.637879,1566.794144,
        1569.477867,1571.77099,1575.385829,1578.108169,1580.194856,
        1581.820666,1583.106578,1589.51315,1589.578955]
    return Dict{String,Any}("positions"=>positions,"temperatures"=>temperatures)
end

"Run every published transport stage; return stage seconds, grid sizes and optional snapshots."
function run_source_sequence(gas,data,case,profile=nothing;save_profiles=false,after_stage::F=nothing,clock_ns::C=time_ns) where {F,C}
    case in ("free","burner","fixed") || throw(ArgumentError("unknown source flame case"))
    free=case=="free";fixed=case=="fixed"
    expected=fixed ? (53,325) : (10,29)
    (gas.n_species,gas.n_reactions)==expected ||
        throw(ArgumentError("supply the stock source mechanism: $(expected[1]) species and $(expected[2]) reactions"))
    modes=free ? ["mass","mass-soret","multi","multi-soret"] : ["mole","multi"]
    seconds=zeros(length(modes));points=zeros(Int,length(modes));snapshots=Dict{String,Any}()
    f=nothing
    for (stage,mode) in enumerate(modes)
        multi=startswith(mode,"multi")
        stage_start_ns=clock_ns()
        begin
            if stage==1
                X=free ? "H2:1.1,O2:1,AR:5" : fixed ? "CH4:.65,O2:1,N2:3.76" : "H2:1.5,O2:1,AR:7"
                f=free ? FreeFlame(gas;T=300.,P=one_atm,X,width=.03) :
                    BurnerFlame(gas;T=fixed ? 373.7 : 373.,P=fixed ? one_atm : .05one_atm,X,
                        mdot=fixed ? .04 : .06,width=fixed ? .01 : .5)
                f.discretization==:conservative || error("source example requires the conservative discretization")
                if fixed
                    set_temperature_profile!(f,vec(profile["positions"]),vec(profile["temperatures"]);relative=false,grid_policy=:adaptive)
                end
            end
            set_transport!(f,multi ? :multicomponent : :mixture_averaged;data,soret=endswith(mode,"soret"),
                flux_gradient_basis=free ? :mass : :mole)
            slope=free ? .06 : fixed ? .1 : .05
            curve=free ? .12 : fixed ? .2 : .1
            if fixed && stage == 1
                # Continue the source's coarse mixture solve to the final
                # resolution before saving it; both solves remain timed.
                solve!(f;ratio=3.,slope=.3,curve=1.)
            end
            solve!(f;ratio=3.,slope,curve)
            if save_profiles
                snapshots[mode]=Dict("grid"=>copy(f.grid),"T"=>temperature(f),"Y"=>mass_fractions(f),
                    "velocity"=>velocity(f),"inlet_Y"=>copy(f.inlet_Y),"state"=>copy(f.state),"P"=>[f.pressure])
            end
        end
        seconds[stage]=(clock_ns()-stage_start_ns)/1e9
        f.converged || error("$case $mode failed")
        points[stage]=length(f.grid)
        after_stage === nothing || after_stage(f,mode)
    end
    return seconds,points,snapshots
end
