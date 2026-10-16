using Arrhenius, Random, LinearAlgebra, YAML

length(ARGS) in (4,5) || error("usage: julia --project scripts/gas_sidecar_parity.jl <repo-root> <nh3-root> <npz-dir> <output.md> [mechanism-filter]")
repo,nh3,npzdir,out=abspath.(ARGS[1:4])
filter_name=length(ARGS)==5 ? strip(ARGS[5]) : ""
inputs=[(joinpath(repo,"mechanism",n),"repo_"*n) for n in ("2S_CH4_CM2.yaml","JP10skeletal.yaml","LiDryer.yaml","gri30.yaml","h2o2.yaml","nDodecane_Reitz.yaml","nc7_ver3.1_mech.yaml")]
append!(inputs,[(joinpath(nh3,n),"nh3_"*n) for n in ("alzueta_2023.yaml","creck_nh3.yaml","glarborg_2018.yaml","h2o2.yaml","nakamura_2017.yaml","shrestha_2020.yaml","stagni_2023.yaml")])
!isempty(filter_name) && filter!(item->occursin(lowercase(filter_name),lowercase(item[2])),inputs)
isempty(inputs) && error("no mechanisms matched filter: $filter_name")
ordered(x::Float64)=(b=reinterpret(UInt64,x); (b>>63)==1 ? ~b : b|(UInt64(1)<<63))
ulp(a,b)=abs(Int128(ordered(Float64(a)))-Int128(ordered(Float64(b))))
lines=["# Native gas YAML parity against Cantera 4.0.0a2 sidecars","","Integer and logical arrays require exact equality. Floating arrays report maximum ULP difference and number of unequal entries. Wdot and thermo are evaluated at matching random mass fractions and 500, 1200, and 2600 K.","","## Array parity","","| Mechanism | Array | Result | Max ULP | Unequal |","|---|---|---|---:|---:|"]
for (path,refname) in inputs
    refpath=joinpath(npzdir,refname*".npz")
    root=YAML.load_file(path); phase=root["phases"][1]
    gas=CreateSolution(path;data_paths=[dirname(path)])
    parsed=Arrhenius._native_gas_arrays(root,phase,gas.species_names;mechanism=path,data_paths=[dirname(path)])
    raw=Arrhenius.npzread(refpath)
    if "molecular_weights" in keys(parsed) && "molecular_weights" in keys(raw)
        for i in eachindex(parsed["molecular_weights"])
            a,b=parsed["molecular_weights"][i],raw["molecular_weights"][i]
            a == b || println("MW mismatch mechanism=$(basename(path)) species=$(gas.species_names[i]) native=$(repr(a)) sidecar=$(repr(b))")
        end
    end
    for key in sort!(collect(intersect(Set(keys(parsed)),Set(keys(raw)))))
        a,b=parsed[key],raw[key]; size(a)==size(b) || (push!(lines,"| $(basename(path)) | $key | SHAPE $(size(a)) != $(size(b)) | — | — |"); continue)
        isempty(a) && continue
        if eltype(a)<:AbstractFloat && eltype(b)<:AbstractFloat
            n=0; maxulp=0
            for (x,y) in zip(a,b); x==y || (n+=1); maxulp=max(maxulp,ulp(x,y)); end
            push!(lines,"| $(basename(path)) | $key | $(n==0 ? "exact" : "ULP") | $maxulp | $n |")
        else
            push!(lines,"| $(basename(path)) | $key | $(a==b ? "exact" : "DIFF") | — | $(count(a .!= b)) |")
        end
    end
    state_result="unavailable"
    try
        mktempdir() do dir
            localpath=joinpath(dir,basename(path)); cp(path,localpath)
            cp(refpath,localpath*".npz")
            side=CreateSolution(localpath;data_paths=[dirname(path)],kinetics_reader=:sidecar)
            maxrel=0.0; exact=true; thermo_exact=true; thermo_rel=0.0
            Random.seed!(20261016)
            for T in (500.0,1200.0,2600.0)
                Y=rand(gas.n_species); Y./=sum(Y)
                Xn=Y2X(gas,Y,1/dot(Y,1 ./ gas.MW)); Xs=Y2X(side,Y,1/dot(Y,1 ./ side.MW))
                an=set_states(gas,T,one_atm,Y); bs=set_states(side,T,one_atm,Y)
                scale=max(maximum(abs,bs;init=0.0),floatmin(Float64))
                maxrel=max(maxrel,maximum(abs.(an.-bs);init=0.0)/scale); exact &= an==bs
                for f in (cal_cp_R,cal_h_RT,cal_s0_R)
                    x=f(gas,T,one_atm,Xn); y=f(side,T,one_atm,Xs)
                    thermo_exact &= x==y
                    thermo_rel=max(thermo_rel,maximum(abs.(x.-y);init=0.0)/max(maximum(abs,y;init=0.0),floatmin(Float64)))
                end
            end
            state_result="wdot_bit_identical=$(exact), wdot_relmax=$(maxrel), thermo_bit_identical=$(thermo_exact), thermo_relmax=$(thermo_rel)"
        end
    catch err
        state_result="unavailable: $(sprint(showerror,err))"
    end
    push!(lines,"| $(basename(path)) | wdot + thermo (3 states) | $state_result | — | — |")
    println(basename(path)," done")
    GC.gc()
end
write(out,join(lines,"\n")*"\n")
println("wrote ",out)
