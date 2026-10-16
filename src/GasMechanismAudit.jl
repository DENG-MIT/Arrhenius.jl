# Lightweight gas-mechanism checks and opt-in NASA7 continuity refit.

"""Audit a gas mechanism without constructing or solving a flame.

The report is a named tuple with `species`, `reactions`, `issues`, and `checks`.
Collision-limit and transport-range findings are explicitly heuristic screens.
"""
function audit_mechanism(path; phase=nothing, data_paths=String[])
    root = YAML.load_file(path)
    phaseidx = phase === nothing ? 1 : findfirst(p -> get(p,"name","")==String(phase),root["phases"])
    phaseidx === nothing && throw(ArgumentError("phase $phase not found"))
    ph = root["phases"][phaseidx]
    solroot = _solution_species_yaml(root,String(path),data_paths)
    names = String.(solroot["phases"][1]["species"])
    root = solroot
    # Resolve selectors and parser-level errors through the same input path as CreateSolution.
    arrays = _native_gas_arrays(root,ph,names;mechanism=String(path),data_paths)
    defs = Dict(String(s["name"])=>s for s in root["species"])
    used = Set{String}(); issues = NamedTuple[]
    nr = size(arrays["reactant_stoich_coeffs"],2)
    compositions = [get(defs[n],"composition",Dict()) for n in names]
    charges = [get(c,"E",0.0) for c in compositions]
    atom_names = sort!(unique!(String[k for c in compositions for k in keys(c) if String(k)!="E"]))
    for j in 1:nr
        lhs=arrays["reactant_stoich_coeffs"][:,j]; rhs=arrays["product_stoich_coeffs"][:,j]
        union!(used,[names[i] for i in eachindex(names) if lhs[i]!=0 || rhs[i]!=0])
        for elem in atom_names
            bal=sum((rhs[i]-lhs[i])*get(compositions[i],elem,0.0) for i in eachindex(names))
            abs(bal)>1e-10 && push!(issues,(level=:error,check=:element_balance,reaction=j,detail="$elem imbalance $bal"))
        end
        qbal=sum((rhs[i]-lhs[i])*charges[i] for i in eachindex(names))
        abs(qbal)>1e-10 && push!(issues,(level=:error,check=:charge_balance,reaction=j,detail="electron-count imbalance $qbal"))
    end
    for n in setdiff(Set(names),used)
        push!(issues,(level=:warning,check=:unused_species,species=n,detail="not referenced by selected reactions"))
    end
    # Equation fingerprints are based on stoichiometry and direction, independent of spacing.
    groups=Dict{String,Vector{Int}}()
    for j in 1:nr
        l=arrays["reactant_stoich_coeffs"][:,j]; r=arrays["product_stoich_coeffs"][:,j]
        rev=arrays["is_reversible"][j]
        key1=join(("$(names[i]):$(l[i]):$(r[i])" for i in eachindex(names) if l[i]!=0 || r[i]!=0),"|")*"/"*string(rev)
        key2=rev ? join(("$(names[i]):$(r[i]):$(l[i])" for i in eachindex(names) if l[i]!=0 || r[i]!=0),"|")*"/true" : key1
        push!(get!(groups,min(key1,key2),Int[]),j)
    end
    reactions=_gas_select_reactions(root,ph,String(path),data_paths)
    for ids in values(groups)
        length(ids)<2 && continue
        flagged=[i for i in ids if get(reactions[i],"duplicate",false)===true]
        if isempty(flagged)
            push!(issues,(level=:warning,check=:duplicate_unflagged_candidates,reactions=ids,detail="same stoichiometry candidates lack duplicate: true; Cantera type/third-body checks are not fully reproduced"))
        elseif length(flagged)<2
            push!(issues,(level=:warning,check=:duplicate_orphan_flag_candidates,reactions=ids,detail="duplicate flag has no same-stoichiometry flagged partner in this screen"))
        end
    end
    for (j,r) in enumerate(reactions)
        rc=get(r,"rate-constant",nothing)
        rc isa AbstractDict && get(rc,"A",0) isa Real && rc["A"]<0 && get(r,"negative-A",false)!==true &&
            push!(issues,(level=:error,check=:negative_A,reaction=j,detail="negative A lacks negative-A: true"))
    end
    thermo_issues = _audit_nasa_continuity(root,ph,names,defs)
    append!(issues,thermo_issues)
    missing_transport=String[]
    for n in names
        spec=defs[n]; tr=get(spec,"transport",nothing)
        tr===nothing && push!(missing_transport,n)
        tr isa AbstractDict || continue
        geom=String(get(tr,"geometry",""))
        geom in ("atom","linear","nonlinear") || push!(issues,(level=:warning,check=:transport_sanity,species=n,detail="unrecognized geometry $geom"))
        for (key,lo,hi) in (("well-depth",10.,1000.),("diameter",1.,10.),("polarizability",0.,100.))
            if haskey(tr,key)
                val=Float64(tr[key]); lo<=val<=hi || push!(issues,(level=:info,check=:transport_sanity,species=n,parameter=key,value=val,screen=(lo,hi),detail="outside broad heuristic screening range"))
            end
        end
    end
    isempty(missing_transport) || push!(issues,(level=:info,check=:missing_transport,species=missing_transport,detail="transport data absent from YAML; needed only by transport calculations"))
    # Binary hard-sphere collision screen; only ordinary elementary bimolecular rates are comparable.
    diam=Dict(n=>Float64(get(get(defs[n],"transport",Dict()),"diameter",NaN))*1e-10 for n in names)
    for j in 1:nr
        r=reactions[j]; get(r,"type","elementary")=="elementary" || continue
        haskey(r,"rate-constant") || continue
        sum(arrays["reactant_orders"][:,j])==2 || continue
        stoich=arrays["reactant_stoich_coeffs"][:,j]
        sum(stoich)==2 && all(x->x>=0 && isinteger(x),stoich) || continue
        participants=findall(x->!iszero(x),stoich)
        length(participants)<=2 || continue
        react=Int[]
        for i in participants; append!(react,fill(i,Int(stoich[i]))); end
        length(react)==2 || continue
        all(isfinite(diam[names[i]]) && diam[names[i]]>0 for i in react) || continue
        A,b,Ea=arrays["Arrhenius_coeffs"][j,:]
        mu=(arrays["molecular_weights"][react[1]]*arrays["molecular_weights"][react[2]] /
            (arrays["molecular_weights"][react[1]]+arrays["molecular_weights"][react[2]]))/1000/6.02214076e26
        d=(diam[names[react[1]]]+diam[names[react[2]]])/2
        worst=(ratio=0.,temperature=0.)
        for T in 300.:300.:3000.
            k=A*T^b*exp(-Ea*4184/R/T)
            kcol=pi*d^2*sqrt(8*1.380649e-23*T/(pi*mu))*6.02214076e26*1000
            ratio=k/kcol
            ratio>worst.ratio && (worst=(ratio=ratio,temperature=T))
        end
        worst.ratio>100 && push!(issues,(level=:warning,check=:collision_limit_screen,reaction=j,ratio=worst.ratio,temperature=worst.temperature,detail="rate exceeds hard-sphere collision estimate by >100x; screening threshold is heuristic"))
    end
    (path=String(path),phase=get(ph,"name",""),species=length(names),reactions=nr,issues=issues,
     checks=(element_balance=true,charge_balance=true,duplicate_reactions=:stoichiometric_candidates,unused_species=true,
             nasa_continuity=true,negative_A=true,missing_transport=true,
             collision_limit=:heuristic_hard_sphere,transport_sanity=:heuristic_ranges))
end

function _audit_nasa_continuity(root,ph,names,defs)
    issues=NamedTuple[]
    for name in names
        th=get(defs[name],"thermo",nothing); th isa AbstractDict || continue
        model=lowercase(String(get(th,"model",""))); model in ("nasa7","nasa9") || continue
        ranges=get(th,"temperature-ranges",Float64[]); data=get(th,"data",Any[])
        length(data)<2 && continue
        mid=Float64(ranges[2]); lo=Float64.(data[1]); hi=Float64.(data[2])
        cp(a,T)=model=="nasa7" ? a[1]+a[2]*T+a[3]*T^2+a[4]*T^3+a[5]*T^4 : a[1]/T^2+a[2]/T+a[3]+a[4]*T+a[5]*T^2+a[6]*T^3+a[7]*T^4
        h(a,T)=model=="nasa7" ? a[1]+a[2]*T/2+a[3]*T^2/3+a[4]*T^3/4+a[5]*T^4/5+a[6]/T : -a[1]/T^2+a[2]*log(T)/T+a[3]+a[4]*T/2+a[5]*T^2/3+a[6]*T^3/4+a[7]*T^4/5+a[8]/T
        s(a,T)=model=="nasa7" ? a[1]*log(T)+a[2]*T+a[3]*T^2/2+a[4]*T^3/3+a[5]*T^4/4+a[7] : -a[1]/(2T^2)-a[2]/T+a[3]*log(T)+a[4]*T+a[5]*T^2/2+a[6]*T^3/3+a[7]*T^4/4+a[9]
        vals=(cp(lo,mid)-cp(hi,mid),h(lo,mid)-h(hi,mid),s(lo,mid)-s(hi,mid))
        dcp(a,T)=model=="nasa7" ? a[2]+2a[3]*T+3a[4]*T^2+4a[5]*T^3 : -2a[1]/T^3-a[2]/T^2+a[4]+2a[5]*T+3a[6]*T^2+4a[7]*T^3
        valslope=dcp(lo,mid)-dcp(hi,mid)
        scale=(max(abs(cp(lo,mid)),abs(cp(hi,mid)),1.0),max(abs(h(lo,mid)),abs(h(hi,mid)),1.0),max(abs(s(lo,mid)),abs(s(hi,mid)),1.0))
        maxrel=maximum(abs(vals[i])/scale[i] for i in 1:3)
        (maxrel>1e-6 || abs(valslope)*mid/max(scale[1],1.0)>1e-6) && push!(issues,(level=:warning,check=:nasa_continuity,species=name,Tmid=mid,delta_cp_R=vals[1],delta_dcp_R_dT=valslope,delta_h_RT=vals[2],delta_s_R=vals[3],max_relative=maxrel))
    end
    issues
end

"""Refit the high-temperature NASA7 cp segment to be C0/C1 continuous.
Returns `(report, yaml_path)`; original coefficients are left unchanged.
The high cp polynomial is the least-squares fit to its original cp values on
an evenly spaced temperature grid subject to matching cp and dcp/dT at Tmid;
h and s integration constants are then anchored continuously at Tmid.
"""
function refit_nasa7_continuity(path, output; phase=nothing, samples=80)
    samples>=8 || throw(ArgumentError("samples must be at least 8"))
    root=YAML.load_file(path); phidx=phase===nothing ? 1 : findfirst(p->get(p,"name","")==String(phase),root["phases"])
    phidx===nothing && throw(ArgumentError("phase $phase not found")); ph=root["phases"][phidx]
    root=_solution_species_yaml(root,String(path),String[])
    names=String.(root["phases"][1]["species"])
    defs=Dict(String(s["name"])=>s for s in root["species"]); changes=NamedTuple[]
    for name in names
        th=get(defs[name],"thermo",nothing); th isa AbstractDict || continue
        lowercase(String(get(th,"model","")))=="nasa7" || continue
        data=get(th,"data",Any[]); ranges=get(th,"temperature-ranges",Float64[])
        length(data)==2 && length(ranges)==3 || continue
        low=Float64.(data[1]); old=Float64.(data[2]); mid=Float64(ranges[2]); tmax=Float64(ranges[3])
        x=collect(range(mid,tmax,length=samples)); z=x./mid
        basis=hcat(ones(samples),z,z.^2,z.^3,z.^4)
        target=[old[1]+old[2]*T+old[3]*T^2+old[4]*T^3+old[5]*T^4 for T in x]
        G=[1. 1. 1. 1. 1.; 0. 1. 2. 3. 4.]
        cp_lo=low[1]+low[2]*mid+low[3]*mid^2+low[4]*mid^3+low[5]*mid^4
        dcp_lo=mid*(low[2]+2low[3]*mid+3low[4]*mid^2+4low[5]*mid^3)
        b=[cp_lo,dcp_lo]; x0=G'*((G*G')\b); Z=nullspace(G)
        c=x0+Z*((basis*Z)\(target-basis*x0))
        high=vcat(c[1],c[2]/mid,c[3]/mid^2,c[4]/mid^3,c[5]/mid^4,old[6],old[7])
        h_lo=low[1]+low[2]*mid/2+low[3]*mid^2/3+low[4]*mid^3/4+low[5]*mid^4/5+low[6]/mid
        s_lo=low[1]*log(mid)+low[2]*mid+low[3]*mid^2/2+low[4]*mid^3/3+low[5]*mid^4/4+low[7]
        high[6]=(h_lo-(high[1]+high[2]*mid/2+high[3]*mid^2/3+high[4]*mid^3/4+high[5]*mid^4/5))*mid
        high[7]=s_lo-(high[1]*log(mid)+high[2]*mid+high[3]*mid^2/2+high[4]*mid^3/3+high[5]*mid^4/4)
        oldcp(T)=old[1]+old[2]*T+old[3]*T^2+old[4]*T^3+old[5]*T^4
        newcp(T)=high[1]+high[2]*T+high[3]*T^2+high[4]*T^3+high[5]*T^4
        oldh(T)=old[1]+old[2]*T/2+old[3]*T^2/3+old[4]*T^3/4+old[5]*T^4/5+old[6]/T
        newh(T)=high[1]+high[2]*T/2+high[3]*T^2/3+high[4]*T^3/4+high[5]*T^4/5+high[6]/T
        olds(T)=old[1]*log(T)+old[2]*T+old[3]*T^2/2+old[4]*T^3/3+old[5]*T^4/4+old[7]
        news(T)=high[1]*log(T)+high[2]*T+high[3]*T^2/2+high[4]*T^3/3+high[5]*T^4/4+high[7]
        Ts=range(mid,tmax,length=4*samples)
        dcp=maximum(abs(newcp(T)-oldcp(T)) for T in Ts)
        dh=maximum(abs(newh(T)-oldh(T)) for T in Ts)
        ds=maximum(abs(news(T)-olds(T)) for T in Ts)
        dg=maximum(abs((newh(T)-news(T))-(oldh(T)-olds(T))) for T in Ts)
        data[2]=high; defs[name]["thermo"]["data"]=data
        push!(changes,(species=name,Tmid=mid,max_abs_delta_cp_R=dcp,max_abs_delta_h_RT=dh,max_abs_delta_s_R=ds,max_abs_delta_g_RT=dg))
    end
    # YAML.jl emits malformed indentation for Cantera's multiline `note` fields;
    # these are comments/annotations, so remove them in the refitted copy.
    function strip_notes!(x)
        if x isa AbstractDict
            pop!(x,"note",nothing)
            for v in values(x); strip_notes!(v); end
        elseif x isa AbstractVector
            foreach(strip_notes!,x)
        end
        x
    end
    strip_notes!(root)
    YAML.write_file(output,root)
    (species_changed=length(changes),changes=changes,output=String(output)),String(output)
end
