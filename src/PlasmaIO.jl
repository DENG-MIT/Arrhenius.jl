const _PLASMA_NA = 6.02214076e26
const _PLASMA_ME = 9.1093837015e-31
const _PLASMA_QE = 1.602176634e-19
const _PLASMA_AW = Dict{String,Float64}(
    "H"=>1.008, "He"=>4.002602, "C"=>12.011, "N"=>14.007,
    "O"=>15.999, "Ar"=>39.95, "E"=>_PLASMA_ME*_PLASMA_NA)

_plasma_error(s) = throw(ArgumentError(s))
function _plasma_float(x, what)
    x isa Real || _plasma_error("$what must be a real number, got $(repr(x))")
    y = Float64(x)
    isfinite(y) || _plasma_error("$what must be finite")
    y
end
function _plasma_quantity(x, default_unit, what)
    if x isa Real
        return _plasma_float(x, what), String(default_unit)
    end
    x isa AbstractString || _plasma_error("$what must be a scalar quantity")
    m = match(r"^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(.*?)\s*$", x)
    m === nothing && _plasma_error("invalid $what: $(repr(x))")
    value = try parse(Float64, m.captures[1]) catch; _plasma_error("invalid $what: $(repr(x))") end
    isfinite(value) || _plasma_error("$what must be finite")
    value, isempty(m.captures[2]) ? String(default_unit) : m.captures[2]
end
_plasma_unit(x) = lowercase(replace(strip(String(x)), " "=>""))
function _plasma_temperature(x, what)
    value, unit = _plasma_quantity(x, "K", what)
    _plasma_unit(unit) in ("k", "kelvin") || _plasma_error("$what must use K")
    value > 0 || _plasma_error("$what must be positive")
    value
end
function _plasma_pressure(x, what)
    value, unit = _plasma_quantity(x, "Pa", what)
    factors = Dict("pa"=>1.0, "kpa"=>1e3, "mpa"=>1e6, "bar"=>1e5,
                   "atm"=>101325.0, "torr"=>101325.0/760)
    key = _plasma_unit(unit)
    haskey(factors,key) || _plasma_error("unsupported $what unit $(repr(unit))")
    value *= factors[key]
    isfinite(value) && value > 0 || _plasma_error("$what must be positive and finite")
    value
end
function _plasma_electron_energy(x, what)
    value, unit = _plasma_quantity(x, "eV", what)
    _plasma_unit(unit) == "ev" || _plasma_error("$what must use eV")
    value > 0 || _plasma_error("$what must be positive")
    value
end
function _plasma_activation(x, default_unit, what)
    value, unit = _plasma_quantity(x, default_unit, what)
    factors = Dict("k"=>R, "kelvin"=>R, "j/kmol"=>1.0, "j/mol"=>1e3,
                   "kj/mol"=>1e6, "cal/mol"=>4184.0, "kcal/mol"=>4.184e6,
                   "ev"=>_PLASMA_QE*_PLASMA_NA)
    key = _plasma_unit(unit)
    haskey(factors,key) || _plasma_error("unsupported $what unit $(repr(unit))")
    value *= factors[key]
    isfinite(value) || _plasma_error("$what converts to a nonfinite value")
    value
end
function _plasma_units(root)
    u = get(root,"units",Dict{Any,Any}()); u isa AbstractDict || _plasma_error("units must be a mapping")
    lu, qu, tu = _plasma_unit(get(u,"length","m")), _plasma_unit(get(u,"quantity","kmol")), _plasma_unit(get(u,"time","s"))
    ae = String(get(u,"activation-energy","J/kmol"))
    lengths=Dict("m"=>1.0,"cm"=>1e-2,"mm"=>1e-3,"km"=>1e3)
    quantities=Dict("kmol"=>1.0,"mol"=>1e-3,"molec"=>1/_PLASMA_NA,"molecule"=>1/_PLASMA_NA,"molecules"=>1/_PLASMA_NA)
    times=Dict("s"=>1.0,"ms"=>1e-3,"us"=>1e-6,"μs"=>1e-6,"ns"=>1e-9,"min"=>60.0)
    haskey(lengths,lu) || _plasma_error("unsupported length unit $(repr(lu))")
    haskey(quantities,qu) || _plasma_error("unsupported quantity unit $(repr(qu))")
    haskey(times,tu) || _plasma_error("unsupported time unit $(repr(tu))")
    _plasma_activation(0.0,ae,"activation energy")
    lengths[lu], quantities[qu], times[tu], ae
end

function _plasma_table(root, context)
    entries=get(root,"species",nothing); entries isa AbstractVector || _plasma_error("$context has no species list")
    table=Dict{String,Any}(); order=String[]
    for (i,e) in enumerate(entries)
        e isa AbstractDict || _plasma_error("$context species $i must be a mapping")
        n=get(e,"name",nothing); n isa AbstractString || _plasma_error("$context species $i has no name")
        n=String(n); haskey(table,n) && _plasma_error("duplicate species $n in $context")
        table[n]=e; push!(order,n)
    end
    table,order
end
function _plasma_import(file,path,data_paths)
    all(d->d isa AbstractString,data_paths) || _plasma_error("data_paths entries must be strings")
    candidates=isabspath(file) ? [normpath(file)] : [normpath(joinpath(dirname(abspath(path)),file)); [normpath(joinpath(String(d),file)) for d in data_paths]...]
    i=findfirst(isfile,candidates)
    i===nothing && _plasma_error("imported species file $(repr(file)) not found; searched "*join(candidates,", "))
    candidates[i]
end
function _plasma_species(root,phase,path,data_paths; source_paths=nothing)
    data_paths isa AbstractVector || _plasma_error("data_paths must be a list")
    selection=get(phase,"species",nothing)
    selection == "all" && (selection = [Dict("species"=>"all")])
    selection isa AbstractVector && !isempty(selection) || _plasma_error("phase needs a species list")
    local_table,local_order=_plasma_table(get(root,"species",nothing) === nothing ? Dict("species"=>Any[]) : root,"mechanism"); names=String[]; defs=Any[]
    function add!(table,order,wanted,context,source=path)
        chosen=wanted=="all" ? order : wanted
        chosen isa AbstractVector || _plasma_error("$context selection must be 'all' or a name list")
        for raw in chosen
            raw isa AbstractString || _plasma_error("species names must be strings")
            name=String(raw); haskey(table,name) || _plasma_error("species $name not found in $context")
            name in names && _plasma_error("species $name selected more than once")
            push!(names,name); push!(defs,table[name])
            source_paths === nothing || push!(source_paths,String(source))
        end
    end
    for (i,item) in enumerate(selection)
        if item isa AbstractString
            add!(local_table,local_order,[item],"mechanism")
            continue
        end
        item isa AbstractDict && length(item)==1 || _plasma_error("phase species item $i must select one section")
        section,wanted=first(item); section=String(section)
        if section=="species"
            add!(local_table,local_order,wanted,"mechanism"); continue
        end
        endswith(section,"/species") || _plasma_error("unsupported species section $section; expected <file>/species")
        file=section[1:end-length("/species")]; imported=_plasma_import(file,path,data_paths)
        imported_root=YAML.load_file(imported); imported_root isa AbstractDict || _plasma_error("$imported is not a YAML mapping")
        table,order=_plasma_table(imported_root,imported); add!(table,order,wanted,imported,imported)
    end
    names,defs
end
function _plasma_species_data(names,defs,elements,atomic_weights)
    atomic_weights isa AbstractDict || _plasma_error("atomic_weights must be a mapping")
    aw=copy(_PLASMA_AW)
    for (k,v) in atomic_weights
        w=_plasma_float(v,"atomic weight $k"); w>0 || _plasma_error("atomic weight $k must be positive"); aw[String(k)]=w
    end
    aw["E"]=_PLASMA_ME*_PLASMA_NA
    eindex=Dict(e=>i for (i,e) in enumerate(elements)); length(eindex)==length(elements) || _plasma_error("duplicate phase elements")
    matrix=zeros(Float64,length(elements),length(names)); mw=zeros(length(names)); electrons=Int[]
    for (j,d) in enumerate(defs)
        comp=get(d,"composition",nothing); comp isa AbstractDict && !isempty(comp) || _plasma_error("species $(names[j]) has no composition")
        nonzero=0
        for (raw,c0) in comp
            e=String(raw); haskey(eindex,e) || _plasma_error("species $(names[j]) uses phase-absent element $e")
            haskey(aw,e) || _plasma_error("missing atomic weight for $e")
            c=_plasma_float(c0,"$(names[j])/$e composition"); e=="E" || c>=0 || _plasma_error("only E composition may be negative")
            matrix[eindex[e],j]=c; mw[j]+=c*aw[e]; c!=0 && (nonzero+=1)
        end
        nonzero==1 && get(comp,"E",0)==1 && push!(electrons,j)
        isfinite(mw[j]) && mw[j]>0 || _plasma_error("species $(names[j]) has nonpositive molecular weight")
    end
    length(electrons)==1 || _plasma_error("exactly one pure {E: 1} electron species is required")
    matrix,mw,only(electrons)
end

function _plasma_name(name,index,electron,context)
    name=="M" && return "M"
    name=="E" && !haskey(index,"E") && return electron
    haskey(index,name) || _plasma_error("$context has unknown species $(repr(name))")
    name
end
function _plasma_side(text,index,electron,context)
    occursin("(+",text) && _plasma_error("$context uses unsupported falloff syntax")
    counts=Dict{String,Int}()
    for raw in split(strip(String(text)),r"\s+\+\s+";keepempty=false)
        term=strip(raw); isempty(term) && _plasma_error("empty term in $context")
        m=match(r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+))\s+(.+)$",term)
        if m===nothing; n=1; rawname=term
        else
            v=parse(Float64,m.captures[1]); isinteger(v) || _plasma_error("fractional stoichiometry in $context")
            v>0 && v<=typemax(Int) || _plasma_error("invalid stoichiometry in $context")
            n=Int(v); rawname=strip(m.captures[2])
        end
        name=_plasma_name(rawname,index,electron,context); counts[name]=get(counts,name,0)+n
    end
    isempty(counts) && _plasma_error("empty reaction side in $context"); counts
end
function _plasma_effmap(r,index,electron,context)
    raw=get(r,"efficiencies",nothing); raw===nothing && return nothing
    raw isa AbstractDict || _plasma_error("$context efficiencies must be a mapping")
    out=Dict{String,Float64}()
    for (k,v) in raw
        name=_plasma_name(String(k),index,electron,context); name=="M" && _plasma_error("M cannot have an efficiency")
        value=_plasma_float(v,"$context efficiency $name"); value>=0 || _plasma_error("negative efficiency in $context"); out[name]=value
    end
    out
end
function _plasma_remove!(side,name,context)
    n=get(side,name,0); n>0 || _plasma_error("missing collider $name in $context")
    n==1 ? delete!(side,name) : (side[name]=n-1)
end
function _plasma_thirdbody!(lhs,rhs,r,kind,index,electron,context)
    emap=_plasma_effmap(r,index,electron,context); hasdefault=haskey(r,"default-efficiency")
    lm,rm=get(lhs,"M",0),get(rhs,"M",0); collider=nothing; generic=false
    if lm!=0 || rm!=0
        kind in ("elementary","electron-collision-plasma") && _plasma_error("$context cannot use M")
        lm==1 && rm==1 || _plasma_error("$context needs exactly one M on both sides")
        collider="M"; generic=true
    else
        shared=[n for n in keys(lhs) if haskey(rhs,n)]
        shared_count=sum((1 + (lhs[n]>1 && rhs[n]>1 ? 1 : 0) for n in shared); init=0)
        if kind=="three-body"
            isempty(shared) && _plasma_error("$context has no shared third-body collider")
            if shared_count==1; collider=only(shared)
            elseif emap!==nothing && length(emap)==1 && first(keys(emap)) in shared; collider=first(keys(emap))
            else _plasma_error("ambiguous third body in $context; specify one efficiency or M") end
        elseif !(kind in ("elementary","electron-collision-plasma")) && shared_count==1 && (sum(values(lhs))==3 || sum(values(rhs))==3)
            collider=only(shared)
        elseif !(kind in ("elementary","electron-collision-plasma")) && shared_count>1 && emap!==nothing
            length(emap)==1 || _plasma_error("ambiguous explicit third body in $context")
            candidate=first(keys(emap))
            candidate in shared || _plasma_error("efficiency does not identify a shared collider in $context")
            collider=candidate
        elseif !(kind in ("elementary","electron-collision-plasma")) && shared_count>1 && hasdefault
            _plasma_error("ambiguous explicit third body in $context")
        end
    end
    if collider===nothing
        (emap!==nothing || hasdefault) && _plasma_error("efficiencies without a resolvable third body in $context")
        return false,ones(Float64,length(index))
    end
    _plasma_remove!(lhs,collider,context); _plasma_remove!(rhs,collider,context)
    if generic
        default=hasdefault ? _plasma_float(r["default-efficiency"],"$context default efficiency") : 1.0
        default>=0 || _plasma_error("negative default efficiency in $context")
        eff=fill(default,length(index)); emap!==nothing && foreach(p->(eff[index[p.first]]=p.second),emap)
    else
        hasdefault && _plasma_float(r["default-efficiency"],"$context default efficiency")!=0 && _plasma_error("specific collider default efficiency must be zero")
        emap!==nothing && any(n!=collider for n in keys(emap)) && _plasma_error("incompatible efficiencies for specific collider $collider")
        eff=zeros(Float64,length(index)); eff[index[collider]]=emap===nothing ? 1.0 : get(emap,collider,1.0)
    end
    true,eff
end
function _plasma_orders(lhs,r,index,electron,context)
    out=zeros(Int,length(index)); foreach(p->(out[index[p.first]]=p.second),lhs)
    raw=get(r,"orders",nothing); raw===nothing && return out
    raw isa AbstractDict || _plasma_error("$context orders must be a mapping")
    allow=get(r,"nonreactant-orders",false); allow isa Bool || _plasma_error("nonreactant-orders must be boolean")
    for (k,v0) in raw
        name=_plasma_name(String(k),index,electron,context); name=="M" && _plasma_error("cannot override M order")
        v=_plasma_float(v0,"$context order $name"); v>=0 || _plasma_error("negative reaction order in $context")
        isinteger(v) || _plasma_error("fractional reaction order in $context")
        haskey(lhs,name) || allow || _plasma_error("nonreactant order $name requires nonreactant-orders: true")
        out[index[name]]=Int(v)
    end
    out
end
function _plasma_parameters(r,kind,order,units,context)
    p=zeros(6); kind=="electron-collision-plasma" && return p
    rc=get(r,"rate-constant",nothing); rc isa AbstractDict || _plasma_error("$context needs rate-constant")
    haskey(rc,"A") || _plasma_error("$context rate is missing A")
    lf,qf,tf,ae=units; A=_plasma_float(rc["A"],"$context A"); A>=0 || _plasma_error("negative A in $context")
    p[1]=A*(lf^3/qf)^(order-1)/tf; isfinite(p[1]) || _plasma_error("nonfinite converted A in $context")
    p[2]=_plasma_float(get(rc,"b",0),"$context b")
    if kind=="two-temperature-plasma"
        haskey(rc,"Ea") && _plasma_error("two-temperature $context must use Ea-gas/Ea-electron")
        p[3]=_plasma_activation(get(rc,"Ea-gas",0),ae,"$context Ea-gas")
        p[4]=_plasma_activation(get(rc,"Ea-electron",0),ae,"$context Ea-electron")
        p[5]=_plasma_float(get(rc,"b-gas",0),"$context b-gas")
        if haskey(rc,"T-inv")
            p[6]=inv(_plasma_temperature(rc["T-inv"],"$context T-inv"))
        end
    else
        any(haskey(rc,k) for k in ("Ea-gas","Ea-electron","b-gas","T-inv")) && _plasma_error("two-temperature fields on ordinary $context")
        p[3]=_plasma_activation(get(rc,"Ea",0),ae,"$context Ea")
    end
    all(isfinite,p) || _plasma_error("nonfinite rate parameters in $context"); p
end

function _plasma_reactions(root,phase)
    selected=get(phase,"reactions",nothing); selected isa AbstractVector || _plasma_error("phase reactions must select named sections")
    out=Any[]
    for (i,item) in enumerate(selected)
        item isa AbstractDict && length(item)==1 || _plasma_error("phase reactions item $i must select one section")
        section,wanted=first(item); entries=get(root,String(section),nothing); entries isa AbstractVector || _plasma_error("missing reaction section $section")
        if wanted=="all"; append!(out,entries); continue end
        wanted isa AbstractVector || _plasma_error("reaction selection must be 'all' or IDs")
        for id in wanted
            hit=findfirst(e->e isa AbstractDict && get(e,"id",nothing)==id,entries)
            hit===nothing && _plasma_error("reaction ID $id not found in $section"); push!(out,entries[hit])
        end
    end
    all(e->e isa AbstractDict,out) || _plasma_error("reaction entries must be mappings"); out
end
function _plasma_initial_x(state,names,index,electron,mw)
    xkeys=[k for k in ("X","mole-fractions") if haskey(state,k)]; ykeys=[k for k in ("Y","mass-fractions") if haskey(state,k)]
    length(xkeys)<=1 && length(ykeys)<=1 && (isempty(xkeys)||isempty(ykeys)) || _plasma_error("state has multiple compositions")
    isempty(xkeys) && isempty(ykeys) && return [i==1 ? 1.0 : 0.0 for i in eachindex(names)]
    mass=!isempty(ykeys); raw=state[only(mass ? ykeys : xkeys)]
    if raw isa AbstractDict
        entries=collect(pairs(raw))
    elseif raw isa AbstractString
        entries=Pair{String,String}[]
        for item in split(raw,',')
            fields=split(item,':';limit=2)
            length(fields)==2 || _plasma_error("state composition item $(repr(item)) must be name:value")
            push!(entries,strip(fields[1])=>strip(fields[2]))
        end
    else
        _plasma_error("state composition must be a mapping or string")
    end
    out=zeros(length(names))
    for (k,v0) in entries
        name=_plasma_name(String(k),index,electron,"state composition")
        if v0 isa AbstractString
            v=try parse(Float64,strip(v0)) catch; _plasma_error("nonnumeric state composition for $name") end
        else
            v=_plasma_float(v0,"state composition for $name")
        end
        isfinite(v) && v>=0 || _plasma_error("invalid state fraction"); out[index[name]]+=v
    end
    sum(out)>0 || _plasma_error("state composition sum must be positive"); mass && (out./=mw); out./=sum(out); out
end
function _plasma_validate(em,mw,S,rp,grid,shape,x)
    all(isfinite,rp) || _plasma_error("nonfinite rate parameters")
    length(grid)>=2 && all(isfinite,grid) && all(>=(0),grid) && all(>(0),diff(grid)) || _plasma_error("EEDF grid must be finite, nonnegative, and strictly increasing")
    isfinite(shape) && shape>0 || _plasma_error("shape-factor must be positive and finite")
    all(isfinite,x) && all(>=(0),x) && isapprox(sum(x),1;atol=1e-14,rtol=0) || _plasma_error("invalid initial mole fractions")
    residual=em*S
    for j in axes(S,2)
        maximum(abs,view(residual,:,j);init=0.0)<=1e-12 || _plasma_error("reaction $j does not conserve elements")
        sj=view(S,:,j); abs(dot(mw,sj))<=1e-12*max(1,sum(abs.(mw.*sj))) || _plasma_error("reaction $j does not conserve mass")
    end
end

function PlasmaMechanism(path::AbstractString; phase=nothing, data_paths=String[], atomic_weights=Dict())
    isfile(path) || _plasma_error("plasma mechanism not found: $path")
    root=YAML.load_file(path); root isa AbstractDict || _plasma_error("mechanism root must be a mapping")
    ph=_eedf_phase(root,phase); get(ph,"thermo",nothing)=="plasma" || _plasma_error("selected phase is not plasma")
    eedf=get(ph,"electron-energy-distribution",nothing); eedf isa AbstractDict || _plasma_error("phase has no EEDF")
    get(eedf,"type",nothing)=="Boltzmann-two-term" && return _plasma_boltzmann_mechanism(root,ph,path,data_paths,atomic_weights)
    get(eedf,"type",nothing)=="isotropic" || _plasma_error("supported EEDF types are isotropic and Boltzmann-two-term")
    raw_elements=get(ph,"elements",nothing); raw_elements isa AbstractVector && !isempty(raw_elements) || _plasma_error("phase needs elements")
    elements=String.(raw_elements); names,defs=_plasma_species(root,ph,path,data_paths)
    em,mw,electron_index=_plasma_species_data(names,defs,elements,atomic_weights)
    index=Dict(n=>i for (i,n) in enumerate(names)); electron=names[electron_index]
    entries=_plasma_reactions(root,ph); ns,nr=length(names),length(entries)
    react=zeros(Int,ns,nr); prod=zeros(Int,ns,nr); orders=zeros(Int,ns,nr); third=falses(nr); eff=ones(ns,nr)
    kinds=zeros(UInt8,nr); params=zeros(6,nr); energies=[Float64[] for _=1:nr]; cross=[Float64[] for _=1:nr]
    units=_plasma_units(root); supported=("arrhenius","elementary","three-body","two-temperature-plasma","electron-collision-plasma")
    for (j,r) in enumerate(entries)
        eq=get(r,"equation",nothing); eq isa AbstractString || _plasma_error("reaction $j has no equation"); ctx="reaction $j $(repr(eq))"
        reversible=get(r,"reversible",false); reversible isa Bool || _plasma_error("$ctx reversible must be boolean")
        (reversible || occursin("<=>",eq)) && _plasma_error("reversible $ctx is unsupported")
        sides=split(eq,"=>"); length(sides)==2 || _plasma_error("$ctx needs exactly one =>")
        rawkind=get(r,"type","arrhenius"); rawkind isa AbstractString || _plasma_error("$ctx type must be a string")
        kind=lowercase(String(rawkind)); kind in supported || _plasma_error("unsupported type $kind in $ctx")
        lhs=_plasma_side(sides[1],index,electron,ctx); rhs=_plasma_side(sides[2],index,electron,ctx)
        third[j],eff[:,j]=_plasma_thirdbody!(lhs,rhs,r,kind,index,electron,ctx)
        orders[:,j]=_plasma_orders(lhs,r,index,electron,ctx)
        foreach(p->(react[index[p.first],j]=p.second),lhs); foreach(p->(prod[index[p.first],j]=p.second),rhs)
        effective=kind in ("arrhenius","elementary","three-body") ? "arrhenius" : kind
        kinds[j]=effective=="arrhenius" ? 0x01 : effective=="two-temperature-plasma" ? 0x02 : 0x03
        params[:,j]=_plasma_parameters(r,effective,sum(orders[:,j])+(third[j] ? 1 : 0),units,ctx)
        if effective=="electron-collision-plasma"
            react[electron_index,j]>0 || _plasma_error("$ctx needs an incident electron")
            energies[j],cross[j]=_eedf_table(r,ctx)
        elseif haskey(r,"energy-levels") || haskey(r,"cross-sections")
            _plasma_error("thermal $ctx cannot contain collision tables")
        end
    end
    S=Float64.(prod-react); rawgrid=get(eedf,"energy-levels",nothing); rawgrid isa AbstractVector || _plasma_error("EEDF has no energy-levels")
    grid=[_plasma_float(v,"EEDF energy level") for v in rawgrid]; shape=_plasma_float(get(eedf,"shape-factor",2),"shape-factor")
    E0=_plasma_electron_energy(get(eedf,"mean-electron-energy",1),"mean-electron-energy")
    state=get(ph,"state",Dict{Any,Any}()); state isa AbstractDict || _plasma_error("phase state must be a mapping")
    T0=_plasma_temperature(get(state,"T",300),"initial temperature"); P0=_plasma_pressure(get(state,"P",101325),"initial pressure")
    x0=_plasma_initial_x(state,names,index,electron,mw); _plasma_validate(em,mw,S,params,grid,shape,x0)
    PlasmaMechanism(String(get(ph,"name","plasma")),ns,nr,names,elements,mw,em,electron_index,
        react,prod,orders,S,third,eff,kinds,params,energies,cross,grid,shape,T0,P0,E0,x0)
end
