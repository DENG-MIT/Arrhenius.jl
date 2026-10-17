# Native gas-phase YAML to the arrays consumed by `Reaction`.
# The representation deliberately mirrors the legacy sidecar schema so the
# kinetics evaluator remains independent of the input format.

# Values from the pinned Cantera Elements.cpp atomicWeightTable. Elements for
# which Cantera uses -1.0 are intentionally excluded: they cannot define a
# positive molecular weight without an isotope-specific composition.
const _GAS_ATOMIC_WEIGHTS = Dict(
    "H"=>1.008, "He"=>4.002602, "Li"=>6.94, "Be"=>9.0121831, "B"=>10.81,
    "C"=>12.011, "N"=>14.007, "O"=>15.999, "F"=>18.998403163, "Ne"=>20.1797,
    "Na"=>22.98976928, "Mg"=>24.305, "Al"=>26.9815384, "Si"=>28.085,
    "P"=>30.973761998, "S"=>32.06, "Cl"=>35.45, "Ar"=>39.95, "K"=>39.0983,
    "Ca"=>40.078, "Sc"=>44.955908, "Ti"=>47.867, "V"=>50.9415, "Cr"=>51.9961,
    "Mn"=>54.938043, "Fe"=>55.845, "Co"=>58.933194, "Ni"=>58.6934,
    "Cu"=>63.546, "Zn"=>65.38, "Ga"=>69.723, "Ge"=>72.630, "As"=>74.921595,
    "Se"=>78.971, "Br"=>79.904, "Kr"=>83.798, "Rb"=>85.4678, "Sr"=>87.62,
    "Y"=>88.90584, "Zr"=>91.224, "Nb"=>92.90637, "Mo"=>95.95, "Ru"=>101.07,
    "Rh"=>102.90549, "Pd"=>106.42, "Ag"=>107.8682, "Cd"=>112.414,
    "In"=>114.818, "Sn"=>118.710, "Sb"=>121.760, "Te"=>127.60,
    "I"=>126.90447, "Xe"=>131.293, "Cs"=>132.90545196, "Ba"=>137.327,
    "La"=>138.90547, "Ce"=>140.116, "Pr"=>140.90766, "Nd"=>144.242,
    "Sm"=>150.36, "Eu"=>151.964, "Gd"=>157.25, "Tb"=>158.925354,
    "Dy"=>162.500, "Ho"=>164.930328, "Er"=>167.259, "Tm"=>168.934218,
    "Yb"=>173.045, "Lu"=>174.9668, "Hf"=>178.49, "Ta"=>180.94788,
    "W"=>183.84, "Re"=>186.207, "Os"=>190.23, "Ir"=>192.217,
    "Pt"=>195.084, "Au"=>196.966570, "Hg"=>200.592, "Tl"=>204.38,
    "Pb"=>207.2, "Bi"=>208.98040, "Th"=>232.0377, "Pa"=>231.03588,
    "U"=>238.02891,
    # Cantera stores the electron separately in its isotope table.
    "E"=>0.0005485799088728283)

_gasfloat(x, what) = x isa Real && !(x isa Bool) && isfinite(Float64(x)) ?
    Float64(x) : throw(ArgumentError("$what must be a finite number"))

# Cantera UnitSystem uses std::pow for dimensional conversion factors.
@inline _gas_pow(x,y)=ccall(:pow,Cdouble,(Cdouble,Cdouble),Float64(x),Float64(y))

function _gas_units(root, local_units=Dict{Any,Any}())
    u = merge(Dict{Any,Any}(),get(root, "units", Dict{Any,Any}()),local_units)
    u isa AbstractDict || throw(ArgumentError("units must be a mapping"))
    length = lowercase(String(get(u,"length","m")))
    quantity = lowercase(String(get(u,"quantity","kmol")))
    time = lowercase(String(get(u,"time","s")))
    lf = get(Dict("m"=>1.0,"cm"=>1e-2,"mm"=>1e-3,"km"=>1e3),length,nothing)
    qf = get(Dict("kmol"=>1.0,"mol"=>1e-3,"molec"=>1/6.02214076e26,
        "molecule"=>1/6.02214076e26,"molecules"=>1/6.02214076e26),quantity,nothing)
    tf = get(Dict("s"=>1.0,"ms"=>1e-3,"us"=>1e-6,"ns"=>1e-9,"min"=>60.0),time,nothing)
    (lf !== nothing && qf !== nothing && tf !== nothing) ||
        throw(ArgumentError("unsupported rate units: length=$length quantity=$quantity time=$time"))
    ae = lowercase(replace(String(get(u,"activation-energy","J/kmol"))," "=>""))
    aef = get(Dict("j/kmol"=>1.0,"kj/kmol"=>1e3,"j/mol"=>1e3,"kj/mol"=>1e6,
        "cal/mol"=>4184.0,"kcal/mol"=>4.184e6,"k"=>R,"kelvin"=>R),ae,nothing)
    aef === nothing && throw(ArgumentError("unsupported activation-energy unit $ae"))
    lf,qf,tf,aef
end

function _gas_quantity(x, default_factor, what)
    x isa Real && return Float64(x)*default_factor
    x isa AbstractString || throw(ArgumentError("$what must be a number or quantity string"))
    m=match(r"^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(.*?)\s*$",x)
    m === nothing && throw(ArgumentError("invalid $what: $(repr(x))"))
    value=parse(Float64,m.captures[1]); unit=lowercase(replace(strip(m.captures[2])," "=>""))
    isempty(unit) && return value*default_factor
    factors=Dict("m"=>1.0,"cm"=>1e-2,"mm"=>1e-3,"km"=>1e3,
        "s"=>1.0,"ms"=>1e-3,"mol"=>1e-3,"kmol"=>1.0,
        "j/kmol"=>1.0,"kj/kmol"=>1e3,"j/mol"=>1e3,"kj/mol"=>1e6,
        "cal/mol"=>4184.0,"kcal/mol"=>4.184e6,"k"=>R,"kelvin"=>R,
        "pa"=>1.0,"kpa"=>1e3,"mpa"=>1e6,"bar"=>1e5,"atm"=>101325.0,"torr"=>101325.0/760)
    # Ea is the only quantity-string used in supported Arrhenius entries here.
    haskey(factors,unit) && return value*factors[unit]
    throw(ArgumentError("unsupported $what unit $unit"))
end

function _gas_arrhenius(data, order, units; name="rate-constant", negative_A=false, validate=true)
    data isa AbstractDict || throw(ArgumentError("$name must be a mapping"))
    Araw=get(data,"A",nothing); Araw === nothing && throw(ArgumentError("$name is missing A"))
    A=Araw isa AbstractString ? 0.0 : _gasfloat(Araw,"$name A")
    # Cantera's native rate basis is m, kmol, s. A quantity string carries its
    # own units; bare values inherit the file's defaults.
    if Araw isa AbstractString
        # Parse common Cantera forms such as cm^3/mol/s, cm6/mol2/s and m^3/kmol/s.
        m=match(r"^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(.*?)\s*$",Araw)
        m === nothing && throw(ArgumentError("invalid $name A quantity"))
        A=parse(Float64,m.captures[1]); unit=lowercase(replace(strip(m.captures[2])," "=>""))
        unit=replace(unit,"**"=>"^")
        # Only rate units consistent with the reaction order are accepted.
        target=3*(order-1)
        len=match(r"^(m|cm|mm)(?:\^(\d+))?/(kmol|mol)(?:\^(\d+))?/(s|ms|us)$",unit)
        if target==0 && unit in ("1/s","1/ms","1/us")
            tf=unit=="1/s" ? 1.0 : unit=="1/ms" ? 1e-3 : 1e-6
            A/=tf
        elseif len !== nothing
            lu,lp,qu,qp,tu=len.captures
            lp=lp===nothing ? 1 : parse(Int,lp); qp=qp===nothing ? 1 : parse(Int,qp)
            lp==target && qp==order-1 || throw(ArgumentError("$name A units do not match reaction order $order"))
            lf=lu=="m" ? 1.0 : lu=="cm" ? 1e-2 : 1e-3
            qf=qu=="kmol" ? 1.0 : 1e-3; tf=tu=="s" ? 1.0 : tu=="ms" ? 1e-3 : 1e-6
            source_factor=1.0
            source_factor *= _gas_pow(lf,lp)
            source_factor *= _gas_pow(qf,-qp)
            source_factor *= _gas_pow(tf,-1.0)
            A = A * source_factor / 1.0
        else
            throw(ArgumentError("unsupported $name A units $(repr(unit))"))
        end
    else
        lf,qf,tf,_=units
        # Match Cantera UnitSystem::convertTo order for the dimensional rate
        # unit: output factor, length, time, then quantity.
        A = A / 1.0
        A *= _gas_pow(lf,3*(order-1))
        A *= _gas_pow(tf,-1.0)
        A *= _gas_pow(qf,1-order)
    end
    _gas_validation(A >= 0 || negative_A, validate, "$name has negative A without negative-A: true")
    b=_gasfloat(get(data,"b",0.0),"$name b")
    # Cantera stores activation energy as Ea/R and reconstructs J/kmol in
    # input_data before the sidecar converts it back to cal/mol.
    Ea_j_kmol=_gas_quantity(get(data,"Ea",0.0),units[4],"$name Ea")
    # Match AnyValue::applyUnits operation order for sidecar parity.
    Ea=(Ea_j_kmol/R)/(1.0/R)/4184.0
    (A,b,Ea)
end

function _gas_pressure(x, root, local_units=Dict{Any,Any}())
    units=merge(Dict{Any,Any}(),get(root,"units",Dict()),local_units)
    default=lowercase(String(get(units,"pressure","Pa")))
    factors=Dict("pa"=>1.0,"kpa"=>1e3,"mpa"=>1e6,"bar"=>1e5,"atm"=>101325.0,"torr"=>101325.0/760)
    if x isa Real
        haskey(factors,default) || throw(ArgumentError("unsupported pressure unit $default"))
        return Float64(x)/1.0*_gas_pow(factors[default],1.0)
    end
    x isa AbstractString || throw(ArgumentError("pressure must be a number or quantity string"))
    m=match(r"^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(.*?)\s*$",x)
    m===nothing && throw(ArgumentError("invalid pressure $(repr(x))"))
    unit=lowercase(strip(m.captures[2])); unit=isempty(unit) ? default : unit
    haskey(factors,unit) || throw(ArgumentError("unsupported pressure unit $unit"))
    parse(Float64,m.captures[1])*factors[unit]
end

@inline function _gas_validation(valid, enabled, message)
    valid && return true
    if enabled
        throw(ArgumentError(message))
    else
        @warn message
        return false
    end
end

@inline function _gas_arrhenius_value(coefficients, T)
    A,b,Ea=coefficients
    A*exp(b*log(T)-Ea*(4184.0/R)/T)
end

function _gas_side(text, species_index, context)
    side=strip(text); isempty(side) && throw(ArgumentError("empty side in $context"))
    side=replace(side,r"\(\s*\+\s*([^)]+)\s*\)"=>"")
    out=zeros(Float64,length(species_index)); collider=nothing
    for raw in split(side,r"\s+\+\s+";keepempty=false)
        term=strip(raw); isempty(term) && continue
        term in ("M","m") && (collider="M"; continue)
        m=match(r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+(.+)$",term)
        coeff=m===nothing ? 1.0 : parse(Float64,m.captures[1]); name=m===nothing ? term : strip(m.captures[2])
        haskey(species_index,name) || throw(ArgumentError("undeclared species $name in $context"))
        out[species_index[name]] += coeff
    end
    out,collider
end

function _gas_reaction_equation(eq, index, context)
    eq isa AbstractString || throw(ArgumentError("$context has no equation"))
    delimiter=occursin("<=>",eq) ? "<=>" : occursin("=>",eq) ? "=>" : occursin("=",eq) ? "=" : nothing
    delimiter===nothing && throw(ArgumentError("invalid reaction equation in $context"))
    parts=split(eq,delimiter;limit=2); length(parts)==2 || throw(ArgumentError("invalid equation in $context"))
    function unwrap(s)
        m=match(r"\(\s*\+\s*([^)]+)\s*\)",s)
        m===nothing ? (s,nothing) : (replace(s,r"\(\s*\+\s*([^)]+)\s*\)"=>""),strip(m.captures[1]))
    end
    left,lc=unwrap(parts[1]); right,rc=unwrap(parts[2])
    lc==rc || (lc===nothing && rc===nothing) || throw(ArgumentError("third-body marker mismatch in $context"))
    react,m1=_gas_side(left,index,context); prod,m2=_gas_side(right,index,context)
    m1==m2 || (m1===nothing && m2===nothing) || throw(ArgumentError("third-body marker mismatch in $context"))
    collider=lc===nothing ? m1 : lc
    react,prod,delimiter!="=>",collider
end

function _gas_select_reactions(root, phase, mechanism, data_paths)
    selectors=get(phase,"reactions","all")
    selected=Any[]
    function append_entries!(entries,wanted,section,import_units=nothing)
        entries isa AbstractVector || throw(ArgumentError("reaction section $section is missing"))
        section_units=length(entries)>0 && entries[1] isa AbstractDict && length(entries[1])==1 && haskey(entries[1],"units")
        defaults=section_units ? entries[1]["units"] : Dict{Any,Any}()
        import_units===nothing && (import_units=Dict{Any,Any}())
        for (i,r) in enumerate(entries)
            section_units && i==1 && continue
            include = wanted=="all" || (wanted isa AbstractVector && get(r,"id",nothing) in wanted)
            include || continue
            merged_units=merge(Dict{Any,Any}(),import_units,defaults,get(r,"units",Dict{Any,Any}()))
            isempty(merged_units) ? push!(selected,r) : push!(selected,merge(r,Dict("units"=>merged_units)))
        end
        if wanted isa AbstractVector
            found=Set(get(r,"id",nothing) for r in entries if r isa AbstractDict)
            all(id->id in found,wanted) || throw(ArgumentError("reaction id not found in $section"))
        end
    end
    if selectors=="all"
        append_entries!(get(root,"reactions",Any[]),"all","reactions")
        return selected
    end
    selectors isa AbstractVector || throw(ArgumentError("phase reactions must be 'all' or a list"))
    for item in selectors
        if item isa AbstractString
            section=String(item); source=root
        else
            item isa AbstractDict && length(item)==1 || throw(ArgumentError("reaction selector must name one section"))
            section,wanted=first(item); section=String(section); source=root
        end
        wanted = item isa AbstractString ? "all" : wanted
        if occursin('/',section)
            file,section=split(section,'/';limit=2); source=YAML.load_file(_plasma_import(file,mechanism,data_paths))
        end
        append_entries!(get(source,section,nothing),wanted,section,source===root ? nothing : get(source,"units",Dict{Any,Any}()))
    end
    selected
end

function _native_gas_arrays(root, phase, names; mechanism="", data_paths=String[], validate=true)
    ns=length(names); index=Dict(String(n)=>i for (i,n) in enumerate(names))
    reactions=_gas_select_reactions(root,phase,mechanism,data_paths); nr=length(reactions)
    react=zeros(ns,nr); prod=zeros(ns,nr); orders=zeros(ns,nr); reversible=falses(nr)
    arr=zeros(nr,3); eff=zeros(ns,nr); low=NTuple{3,Float64}[]; troe=NTuple{4,Float64}[]
    three=Int[]; fall=Int[]; falltroe=Int[]; bmidx=Int[]; bm=zeros(0,4)
    pri=Int[]; pci=Int[]; pgo=[1]; pressures=Float64[]; pro=[1]; parr=NTuple{3,Float64}[]
    lmr_offsets=[1]; lmr_species=Int[]; lmr_eps=NTuple{3,Float64}[]
    for (j,r) in enumerate(reactions)
        ctx="reaction $j ($(get(r,"equation","?")))"
        lhs,rhs,rev,collider=_gas_reaction_equation(get(r,"equation",nothing),index,ctx)
        typ=String(get(r,"type","elementary")); units=_gas_units(root,get(r,"units",Dict{Any,Any}()))
        supported_types=("elementary","three-body","falloff","chemically-activated",
            "linear-Burke","pressure-dependent-Arrhenius","Blowers-Masel","Blowers-Masel-rate")
        typ in supported_types || throw(ArgumentError("unsupported reaction rate type $typ in $ctx"))
        explicit=collider!==nothing
        typ=="chemically-activated" && throw(ArgumentError("chemically-activated rates are not supported in $ctx"))
        any(haskey(r,k) for k in ("SRI","Tsang","Chebyshev")) &&
            throw(ArgumentError("SRI, Tsang, and Chebyshev rates are not supported in $ctx"))
        falloff=typ=="falloff"
        shared=findall((lhs .> 0) .& (rhs .> 0))
        if !explicit && !falloff && length(shared)==1 && lhs[only(shared)]>=1 && rhs[only(shared)]>=1 && sum(lhs)!=sum(rhs)
            collider=names[only(shared)]
            lhs[only(shared)]-=1; rhs[only(shared)]-=1
        elseif !explicit && typ=="three-body" && isempty(shared)
            throw(ArgumentError("cannot identify explicit third body in $ctx"))
        end
        plog_pseudo_collider=typ=="pressure-dependent-Arrhenius" && explicit && collider!==nothing && collider!="M"
        invalid_plog_m = typ=="pressure-dependent-Arrhenius" && collider=="M"
        _gas_validation(!invalid_plog_m, validate,
            "pressure-dependent-Arrhenius reaction has a superfluous M collider in $ctx")
        plog_pseudo_collider |= invalid_plog_m && !validate
        thirdbody=typ=="three-body" || (explicit && !plog_pseudo_collider) ||
            (!falloff && collider!==nothing && !plog_pseudo_collider)
        react[:,j]=lhs; prod[:,j]=rhs; reversible[j]=rev
        orders[:,j]=lhs
        raworders=get(r,"orders",Dict())
        raworders isa AbstractDict || throw(ArgumentError("$ctx orders must be a mapping"))
        rev && !isempty(raworders) && throw(ArgumentError("Reaction orders may only be given for irreversible reactions in $ctx"))
        for (n,o) in raworders
            haskey(index,String(n)) || throw(ArgumentError("unknown order species $n in $ctx"))
            lhs[index[String(n)]]>0 || get(r,"nonreactant-orders",false)===true ||
                throw(ArgumentError("nonreactant order for $n in $ctx requires nonreactant-orders: true"))
            value=_gasfloat(o,"$ctx order")
            value>=0 || get(r,"negative-orders",false)===true ||
                throw(ArgumentError("negative order for $n in $ctx requires negative-orders: true"))
            orders[index[String(n)],j]=value
        end
        if thirdbody || falloff
            default=collider!==nothing && collider!="M" ? 0.0 : _gasfloat(get(r,"default-efficiency",1.0),"$ctx default-efficiency")
            eff[:,j].=default
            collider!==nothing && collider!="M" && (haskey(index,collider) || throw(ArgumentError("unknown collider $collider in $ctx")); eff[index[collider],j]=1.0)
            for (n,x) in get(r,"efficiencies",Dict())
                haskey(index,String(n)) || throw(ArgumentError("unknown efficiency species $n in $ctx"))
                eff[index[String(n)],j]=_gasfloat(x,"$ctx efficiency")
            end
        end
        order=sum(orders[:,j])+(thirdbody && !falloff && typ!="linear-Burke" ? 1 : 0)
        if falloff
            push!(fall,j)
            negative_A=get(r,"negative-A",false)===true
            arr[j,:].=_gas_arrhenius(get(r,"high-P-rate-constant",nothing),order,units;name="high-P-rate-constant in $ctx",negative_A,validate)
            lowrate=_gas_arrhenius(get(r,"low-P-rate-constant",nothing),order+1,units;name="low-P-rate-constant in $ctx",negative_A,validate)
            _gas_validation(arr[j,1]*lowrate[1]>=0, validate,
                "falloff high/low pre-exponential factors must have the same sign in $ctx")
            push!(low,lowrate)
            if haskey(r,"Troe")
                t=r["Troe"]
                t1=_gasfloat(t["T1"],"Troe T1"); t3=_gasfloat(t["T3"],"Troe T3")
                # Cantera stores reciprocal T1/T3 internally and exports their reciprocals.
                push!(troe,(_gasfloat(t["A"],"Troe A"),1.0/(1.0/t1),_gasfloat(get(t,"T2",1e30),"Troe T2"),1.0/(1.0/t3)))
                push!(falltroe,length(troe))
            else push!(falltroe,-1) end
        elseif typ=="linear-Burke" || typ=="pressure-dependent-Arrhenius"
            rates=Any[]
            if typ=="linear-Burke"
                colliders=get(r,"colliders",Any[])
                !isempty(colliders) || throw(ArgumentError("$ctx linear-Burke reaction needs colliders"))
                base_collider=colliders[1]
                get(base_collider,"name",nothing)=="M" && get(base_collider,"type",nothing)=="pressure-dependent-Arrhenius" && !haskey(base_collider,"eig0") || throw(ArgumentError("unsupported linear-Burke base collider in $ctx"))
                append!(rates,get(base_collider,"rate-constants",Any[]))
                base=get(base_collider,"efficiency",Dict("A"=>1.0,"b"=>0.0,"Ea"=>0.0))
                baseA=_gasfloat(get(base,"A",1.0),"$ctx base efficiency A")
                baseb=_gasfloat(get(base,"b",0.0),"$ctx base efficiency b")
                baseEa=_gas_quantity(get(base,"Ea",0.0),units[4],"$ctx base efficiency Ea")/4184.0
                for collider_item in colliders[2:end]
                    get(collider_item,"type",nothing)===nothing && haskey(collider_item,"efficiency") || throw(ArgumentError("unsupported linear-Burke collider in $ctx"))
                    haskey(index,String(collider_item["name"])) || throw(ArgumentError("unknown linear-Burke collider in $ctx"))
                    fit=collider_item["efficiency"]
                    push!(lmr_species,index[String(collider_item["name"])])
                    push!(lmr_eps,(_gasfloat(get(fit,"A",1.0),"$ctx efficiency A")/baseA,
                        _gasfloat(get(fit,"b",0.0),"$ctx efficiency b")-baseb,
                        _gas_quantity(get(fit,"Ea",0.0),units[4],"$ctx efficiency Ea")/4184.0-baseEa))
                end
            else
                append!(rates,get(r,"rate-constants",Any[]))
            end
            grouped=Dict{Float64,Vector{NTuple{3,Float64}}}()
            for item in rates
                p=_gas_pressure(item["P"],root,get(r,"units",Dict{Any,Any}()))
                isfinite(p) && p>0 || throw(ArgumentError("$ctx PLOG pressure must be finite and positive"))
                # Cantera groups PLOG rates by log(P), then exports exp(log(P)).
                p=exp(log(p))
                push!(get!(grouped,p,NTuple{3,Float64}[]),
                    _gas_arrhenius(item,order,units;name="PLOG rate in $ctx",negative_A=true,validate))
            end
            isempty(grouped) && throw(ArgumentError("$ctx PLOG has no rates"))
            for (pressure,group) in grouped
                bad_temperatures=filter(T->!(sum(_gas_arrhenius_value(rate,T) for rate in group)>0),
                    (300.0,500.0,1000.0,2000.0,5000.0,10000.0))
                isempty(bad_temperatures) || _gas_validation(false,validate,
                    "PLOG rate sum must be positive at P=$pressure Pa and T=$(join(bad_temperatures, ", ")) K in $ctx")
            end
            ps=sort!(collect(keys(grouped))); arr[j,:].=first(grouped[first(ps)])
            push!(pri,j)
            if typ=="linear-Burke" || plog_pseudo_collider
                push!(pci,0)
            elseif collider!==nothing && collider!="M"
                push!(pci,index[collider])
            else
                efficiencies=get(r,"efficiencies",Dict()); defaults=get(r,"default-efficiency",1.0)
                named=findfirst(pair->last(pair)==1.0,collect(pairs(efficiencies)))
                push!(pci,defaults==0.0 && length(efficiencies)==1 && named!==nothing ? index[String(first(collect(pairs(efficiencies))[named]))] : 0)
            end
            for p in ps
                push!(pressures,p); append!(parr,grouped[p]); push!(pro,length(parr)+1)
            end
            push!(pgo,length(pressures)+1)
            push!(lmr_offsets,length(lmr_species)+1)
        elseif haskey(r,"rate-constant") && typ in ("elementary","three-body","Blowers-Masel","Blowers-Masel-rate")
            rc=r["rate-constant"]
            arr[j,:].=_gas_arrhenius(rc,order,units;name="rate-constant in $ctx",negative_A=get(r,"negative-A",false)===true,validate)
            if typ=="Blowers-Masel" || typ=="Blowers-Masel-rate"
                push!(bmidx,j); bm=vcat(bm,reshape([arr[j,1],arr[j,2],_gas_quantity(get(rc,"Ea0",0.0),units[4],"$ctx Ea0"),_gas_quantity(get(rc,"w",0.0),units[4],"$ctx w")],1,4))
            end
            thirdbody && push!(three,j)
        else
            throw(ArgumentError("unsupported reaction rate type $typ in $ctx"))
        end
    end
    Dict{String,Any}(
        "molecular_weights"=>[_gas_molecular_weight(root,phase,String(n)) for n in names],
        "reactant_stoich_coeffs"=>react,"product_stoich_coeffs"=>prod,"reactant_orders"=>orders,
        "is_reversible"=>reversible,"Arrhenius_coeffs"=>arr,"efficiencies_coeffs"=>eff,
        "Arrhenius_A0"=>[x[1] for x in low],"Arrhenius_b0"=>[x[2] for x in low],"Arrhenius_Ea0"=>[x[3] for x in low],
        "Troe_A"=>[x[1] for x in troe],"Troe_T1"=>[x[2] for x in troe],"Troe_T2"=>[x[3] for x in troe],"Troe_T3"=>[x[4] for x in troe],
        "index_three_body"=>three,"index_falloff"=>fall,"index_falloff_Troe"=>falltroe,
        "BlowersMasel_reaction_indices"=>bmidx,"BlowersMasel_coefficients"=>bm,
        "Plog_reaction_indices"=>pri,"Plog_collider_indices"=>pci,"Plog_group_offsets"=>pgo,
        "Plog_pressures"=>pressures,"Plog_rate_offsets"=>pro,"Plog_Arrhenius"=>isempty(parr) ? zeros(0,3) : reduce(vcat,(permutedims(collect(x)) for x in parr)),
        "Plog_LMR_offsets"=>lmr_offsets,"Plog_LMR_species"=>lmr_species,
        "Plog_LMR_eps"=>isempty(lmr_eps) ? zeros(0,3) : reduce(vcat,(permutedims(collect(x)) for x in lmr_eps)))
end

function _gas_molecular_weight(root,phase,name)
    defs=get(root,"species",Any[]); i=findfirst(s->get(s,"name","")==name,defs)
    i===nothing && throw(ArgumentError("species definition not found for $name"))
    comp=defs[i]["composition"]
    total=0.0
    # Preserve the reference summation sequence by using the selected phase's element order.
    element_order=get(phase,"elements",get(root,"elements",sort!(collect(keys(comp)))))
    element_order isa AbstractString && lowercase(element_order)=="all" &&
        (element_order=get(root,"elements",sort!(collect(keys(comp)))))
    for e in element_order
        element=String(e); haskey(comp,element) || continue
        haskey(_GAS_ATOMIC_WEIGHTS,element) || throw(ArgumentError("unknown element $element; supply an atomic weight table"))
        total += _gasfloat(comp[e],"$name composition")*_GAS_ATOMIC_WEIGHTS[element]
    end
    total
end

_read_gas_kinetics(root,phase,names;kwargs...)=_native_gas_arrays(root,phase,names;kwargs...)
