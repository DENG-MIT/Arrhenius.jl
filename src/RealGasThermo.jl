# Redlich–Kwong expressions adapted from Cantera RedlichKwongMFTP.cpp,
# revision 726522be4e2a13454d8415b7ef799d621f665cf3.
# https://github.com/Cantera/cantera/blob/726522be4e2a13454d8415b7ef799d621f665cf3/src/thermo/RedlichKwongMFTP.cpp
#
# Copyright (c) 2001-2009, California Institute of Technology
# All rights reserved.
# Copyright (c) 2009 Sandia Corporation. Under the terms of Contract
# AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain
# rights in this software.
# Copyright (c) 2011-2026, Cantera Developers. All rights reserved.
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# - Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# - Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# - Neither the name of the California Institute of Technology, Sandia
#   Corporation nor the names of other contributors may be used to endorse or
#   promote products derived from this software without specific prior written
#   permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""
    RedlichKwongThermo(yaml_or_path; phase=nothing, molecular_weights=nothing,
                      critical_properties=nothing)

Native single-phase Redlich–Kwong mixture model. Species use the existing ideal
reference thermodynamics and `a(T)=a0+a1*T`, quadratic attraction mixing and
linear covolume mixing. Explicit `binary-a` coefficients override geometric
means. YAML units and species critical parameters are accepted. An optional
critical-properties YAML file supplies missing critical data; none are guessed.
Pure-species a0 and b must be nonnegative; a1 and binary coefficients may be
signed. Unlike-species coefficients use the nonnegative geometric mean while
pure-species a1 retains its sign. Opposite-sign a1 pairs require explicit
`binary-a` coefficients because their geometric mean is not real.

The standalone constructor needs no kinetic sidecar or Cantera installation.
Default atomic weights cover H, He, C, N, O, S, Cl and Ar. For other elements or
isotopes, supply molecular weights in kg/kmol. No phase-equilibrium or
Peng–Robinson calculation is implied by this model.
"""
struct RedlichKwongThermo{T<:AbstractFloat} <: Thermo
    reference::IdealGasThermo{T}
    species_names::Vector{String}
    MW::Vector{T}
    a0::Matrix{T}
    a1::Matrix{T}
    b::Vector{T}
    sqrt_a0::Vector{T}
    sqrt_a1::Vector{T}
    geometric_mixing::Bool
end

function _rk_quantity(value, dimension, units)
    length_units = Dict("m"=>1.0,"cm"=>1e-2,"mm"=>1e-3)
    amount_units = Dict("kmol"=>1.0,"mol"=>1e-3,"mole"=>1e-3)
    if value isa AbstractString
        matched = match(r"^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(.*?)\s*$",value)
        isnothing(matched) && throw(ArgumentError("invalid RK quantity: $value"))
        number,unit = parse(Float64,matched[1]),replace(String(matched[2])," "=>"")
        if isempty(unit)
            return _rk_quantity(number,dimension,units)
        end
        allowed = dimension == :b ? ("m^3/kmol","cm^3/mol") : dimension == :a0 ?
            ("Pa*m^6/kmol^2*K^0.5","Pa*m^6*K^0.5/kmol^2","Pa*cm^6/mol^2*K^0.5","bar*cm^6/mol^2*K^0.5") :
            ("Pa*m^6/kmol^2/K^0.5","Pa*cm^6/mol^2/K^0.5","bar*cm^6/mol^2/K^0.5")
        unit in allowed || throw(ArgumentError("unsupported explicit RK unit: $unit"))
        return number * (occursin("cm",unit) ? (dimension == :b ? 1e-3 : 1e-6) : 1.0) *
            (startswith(unit,"bar") ? 1e5 : 1.0)
    end
    lunit,q = String(get(units,"length","m")),String(get(units,"quantity","kmol"))
    haskey(length_units,lunit) && haskey(amount_units,q) ||
        throw(ArgumentError("unsupported RK length or quantity unit"))
    get(units,"temperature","K") == "K" || throw(ArgumentError("RK temperature must use kelvin"))
    volume = length_units[lunit]^3/amount_units[q]
    return Float64(value)*(dimension == :b ? volume :
        _thermo_quantity(1.0,:pressure,units)*volume^2)
end

function _rk_pair(value,units)
    if value isa AbstractVector
        length(value) == 2 || throw(ArgumentError("RK a must be scalar or [a0,a1]"))
        return _rk_quantity(value[1],:a0,units),_rk_quantity(value[2],:a1,units)
    end
    return _rk_quantity(value,:a0,units),0.0
end

function RedlichKwongThermo(source;phase=nothing,molecular_weights=nothing,critical_properties=nothing)
    yaml = source isa AbstractString ? YAML.load_file(source) : source
    phases = yaml["phases"]
    index = isnothing(phase) ? findfirst(p -> p["thermo"] == "Redlich-Kwong",phases) :
        findfirst(p -> p["name"] == phase,phases)
    isnothing(index) && throw(ArgumentError("Redlich-Kwong phase not found"))
    selected = phases[index]
    selected["thermo"] == "Redlich-Kwong" || throw(ArgumentError("selected phase must use Redlich-Kwong"))
    names = String.(selected["species"])
    length(unique(names)) == length(names) || throw(ArgumentError("duplicate species names"))
    specs = Dict(s["name"]=>s for s in yaml["species"])
    chosen = [specs[name] for name in names]
    units = get(yaml,"units",Dict())
    reference_doc = copy(yaml)
    reference_doc["phases"] = [selected]
    reference = IdealGasThermo(reference_doc)
    n = length(names)
    weights = Dict("H"=>1.008,"He"=>4.002602,"C"=>12.011,"N"=>14.007,
        "O"=>15.999,"S"=>32.06,"Cl"=>35.45,"Ar"=>39.95)
    if isnothing(molecular_weights)
        all(all(haskey(weights,e) for e in keys(s["composition"])) for s in chosen) ||
            throw(ArgumentError("supply molecular_weights for elements outside H/He/C/N/O/S/Cl/Ar"))
        mw = [sum(weights[e]*count for (e,count) in s["composition"]) for s in chosen]
    else
        mw = Float64.(molecular_weights)
    end
    length(mw) == n && all(isfinite,mw) && all(>(0),mw) || throw(ArgumentError("invalid molecular weights"))
    critical = isnothing(critical_properties) ? Dict() : critical_properties isa AbstractString ?
        YAML.load_file(critical_properties) : critical_properties
    critspec = Dict(uppercase(s["name"])=>s for s in get(critical,"species",[]))
    diag0,diag1,b = zeros(n),zeros(n),zeros(n)
    eos_data = Vector{Any}(undef,n)
    for (i,spec) in enumerate(chosen)
        eos = get(spec,"equation-of-state",Dict())
        eoslist = eos isa AbstractVector ? eos : [eos]
        ei = findfirst(e -> get(e,"model","") == "Redlich-Kwong",eoslist)
        eos = isnothing(ei) ? Dict() : eoslist[ei]
        eos_data[i] = eos
        local_units = merge(units,get(spec,"units",Dict()),get(eos,"units",Dict()))
        if haskey(eos,"a") && haskey(eos,"b")
            diag0[i],diag1[i] = _rk_pair(eos["a"],local_units)
            b[i] = _rk_quantity(eos["b"],:b,local_units)
        else
            from_db = !haskey(spec,"critical-parameters")
            node = from_db ? get(critspec,uppercase(names[i]),Dict()) : spec
            haskey(node,"critical-parameters") || throw(ArgumentError("missing RK or critical parameters for $(names[i])"))
            crit = node["critical-parameters"]
            cu = merge(from_db ? get(critical,"units",Dict()) : units,
                get(node,"units",Dict()),get(crit,"units",Dict()))
            Tc = _thermo_quantity(crit["critical-temperature"],:temperature,cu)
            Pc = _thermo_quantity(crit["critical-pressure"],:pressure,cu)
            Tc > 0 && Pc > 0 && isfinite(Tc) && isfinite(Pc) || throw(ArgumentError("invalid critical parameters"))
            diag0[i] = 0.427480233540*R^2*Tc^2.5/Pc
            b[i] = 0.0866403499650*R*Tc/Pc
        end
    end
    all(isfinite,vcat(diag0,diag1,b)) && all(>=(0),diag0) && all(>=(0),b) ||
        throw(ArgumentError("RK coefficients must be finite, with nonnegative pure-species a0 and b"))
    all(i -> b[i] > 0 || (diag0[i] == 0 && diag1[i] == 0),1:n) ||
        throw(ArgumentError("nonzero attraction requires positive covolume"))
    sqrt0,sqrt1 = sqrt.(diag0),sqrt.(abs.(diag1))
    a0,a1 = sqrt0*sqrt0',sqrt1*sqrt1'
    geometric = all(>=(0),diag1)
    if !geometric
        # Cantera preserves signed diagonal slopes, but uses sqrt(a1_i*a1_j)
        # for unlike species. Such a matrix no longer has the rank-one form
        # used by the nonnegative fast path. Defer non-real pairs until after
        # explicit binary coefficients have been applied.
        for j in 1:n, i in 1:n
            if i == j
                a1[i,j] = diag1[i]
            elseif !iszero(diag1[i]) && !iszero(diag1[j]) && signbit(diag1[i]) != signbit(diag1[j])
                a1[i,j] = NaN
            end
        end
    end
    seen = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}()
    for (i,eos) in enumerate(eos_data)
        bin = get(eos,"binary-a",Dict())
        bu = merge(units,get(chosen[i],"units",Dict()),get(eos,"units",Dict()),get(bin,"units",Dict()))
        for (name,value) in bin
            name == "units" && continue
            j = findfirst(==(name),names)
            isnothing(j) && throw(ArgumentError("unknown binary-a species: $name"))
            pair = _rk_pair(value,bu)
            all(isfinite,pair) || throw(ArgumentError("binary-a must be finite"))
            key = minmax(i,j)
            haskey(seen,key) && seen[key] != pair && throw(ArgumentError("conflicting binary-a coefficients"))
            seen[key] = pair
            a0[i,j] = a0[j,i] = pair[1]
            a1[i,j] = a1[j,i] = pair[2]
            geometric = false
        end
    end
    all(isfinite,a1) || throw(ArgumentError(
        "opposite-sign RK a1 pairs require explicit binary-a coefficients"))
    return RedlichKwongThermo(reference,names,mw,a0,a1,b,sqrt0,sqrt1,geometric)
end

"Reusable property arrays; use one workspace per concurrent caller."
struct RedlichKwongWorkspace{T}
    Ak::Vector{T}
    dAk::Vector{T}
    cp0::Vector{T}
    h0::Vector{T}
    s0::Vector{T}
    u_TV::Vector{T}
    lnphi::Vector{T}
    hbar::Vector{T}
    vbar::Vector{T}
end
RedlichKwongWorkspace(model::RedlichKwongThermo,::Type{T}=Float64) where {T} =
    RedlichKwongWorkspace((zeros(T,length(model.MW)) for _ in 1:9)...)

function _rk_mixing!(w,m,T,X)
    if m.geometric_mixing
        p0,p1 = dot(m.sqrt_a0,X),dot(m.sqrt_a1,X)
        @inbounds for i in eachindex(X)
            w.dAk[i] = m.sqrt_a1[i]*p1
            w.Ak[i] = m.sqrt_a0[i]*p0+T*w.dAk[i]
        end
    else
        mul!(w.dAk,m.a1,X)
        mul!(w.Ak,m.a0,X)
        @. w.Ak += T*w.dAk
    end
    return dot(X,w.Ak),dot(X,w.dAk),dot(X,m.b)
end

function _rk_composition(m,X,basis)
    basis in (:mole,:mass) || throw(ArgumentError("basis must be :mole or :mass"))
    if X isa AbstractString
        amounts = Dict{String,Float64}()
        for entry in split(X,',')
            parts = strip.(split(entry,':';limit=2))
            value = length(parts) == 1 ? 1.0 : parse(Float64,parts[2])
            isfinite(value) && value >= 0 || throw(ArgumentError("composition amounts must be finite and nonnegative"))
            amounts[parts[1]] = get(amounts,parts[1],0.0)+value
        end
        X = amounts
    end
    if X isa AbstractDict
        result = zeros(length(m.MW))
        for (name,value) in X
            i = findfirst(==(String(name)),m.species_names)
            isnothing(i) && throw(ArgumentError("unknown species: $name"))
            result[i] = value
        end
    else
        result = Float64.(X)
    end
    length(result) == length(m.MW) || throw(DimensionMismatch("one composition per species required"))
    all(isfinite,result) && all(>=(0),result) && maximum(result) > 0 || throw(ArgumentError("invalid composition"))
    result ./= maximum(result)
    basis == :mass && (result ./= m.MW)
    return result ./ sum(result)
end

@inline _rk_pressure(T,v,a,b) = R*T/(v-b)-a/(sqrt(T)*v*(v+b))
@inline _rk_dpdv(T,v,a,b) = -R*T/(v-b)^2+a*(2v+b)/(sqrt(T)*v^2*(v+b)^2)

function _rk_volume(T,P,a,b,root)
    root in (:gas,:liquid) || throw(ArgumentError("root must be :gas or :liquid"))
    A,B = a*P/(R^2*T^2.5),b*P/(R*T)
    c,d = A-B-B^2,-A*B
    p,q = c-1/3,-2/27+c/3+d
    delta = (q/2)^2+(p/3)^3
    if delta >= 0
        z = cbrt(-q/2+sqrt(delta))+cbrt(-q/2-sqrt(delta))+1/3
        candidates = [z]
    else
        radius = 2sqrt(-p/3)
        angle = acos(clamp(3q/(p*radius),-1,1))/3
        candidates = [radius*cos(angle-2pi*k/3)+1/3 for k in 0:2]
    end
    filter!(z -> z > B && isfinite(z),candidates)
    isempty(candidates) && throw(DomainError((T,P),"no physical RK volume root"))
    v = (root == :gas ? maximum(candidates) : minimum(candidates))*R*T/P
    # Polish the selected cubic root using the unexpanded EOS.
    for _ in 1:8
        residual = _rk_pressure(T,v,a,b)-P
        abs(residual) <= 2e-13*P && break
        step = residual/_rk_dpdv(T,v,a,b)
        vnext = v-step
        vnext > b && isfinite(vnext) || break
        v = vnext
    end
    abs(_rk_pressure(T,v,a,b)-P) <= 2e-9*P || throw(ErrorException("RK volume solve did not converge"))
    return v
end

"""
    redlich_kwong_properties!(work, model, T, rho, X)

Evaluate a homogeneous RK state from kelvin, kg/m³ and normalized mole
fractions. The returned h/u/s/cp/cv are J/kmol or J/kmol/K; h_mass etc. use kg.
`u_TV[k] = ∂U/∂n[k]|T,V,n[j≠k]` is the energy coefficient for reacting
constant-volume systems. It differs from ordinary partial molar internal energy.
Returned vector properties alias workspace storage until its next evaluation.
Only mechanically stable states with positive pressure and heat capacity are
accepted. Species polynomials are evaluated at the supplied temperature without
clamping, as in Cantera. Users must check their fit ranges and the accuracy of
the cubic EOS for their application.
"""
function redlich_kwong_properties!(w::RedlichKwongWorkspace,m::RedlichKwongThermo,T,rho,X)
    n = length(m.MW)
    length(X) == n && length(w.Ak) == n || throw(DimensionMismatch("RK workspace/composition size mismatch"))
    isfinite(T) && T > 0 && isfinite(rho) && rho > 0 || throw(ArgumentError("positive finite T and rho required"))
    all(isfinite,X) && all(>=(0),X) && abs(sum(X)-1) <= 1e-10 || throw(ArgumentError("normalized nonnegative mole fractions required"))
    a,aT,b = _rk_mixing!(w,m,T,X)
    MW = dot(X,m.MW)
    v = MW/rho
    v > b || throw(DomainError(v,"molar volume must exceed covolume"))
    P = _rk_pressure(T,v,a,b)
    dpdv = _rk_dpdv(T,v,a,b)
    P > 0 && dpdv < 0 || throw(DomainError((T,rho),"RK state must have positive pressure and mechanical stability"))
    sqt,RT = sqrt(T),R*T
    for i in 1:n
        if isnothing(m.reference.extra)
            coeff = _nasa7_coefficients(m.reference,i,T)
            w.cp0[i] = coeff[i,1]+T*(coeff[i,2]+T*(coeff[i,3]+T*(coeff[i,4]+T*coeff[i,5])))
            w.h0[i] = coeff[i,1]+T*(coeff[i,2]/2+T*(coeff[i,3]/3+T*(coeff[i,4]/4+T*coeff[i,5]/5)))+coeff[i,6]/T
            w.s0[i] = coeff[i,1]*log(T)+T*(coeff[i,2]+T*(coeff[i,3]/2+T*(coeff[i,4]/3+T*coeff[i,5]/4)))+coeff[i,7]
        else
            w.cp0[i],w.h0[i],w.s0[i] = _extended_thermo(m.reference,i,T)
        end
    end
    Z = P*v/RT
    L = log1p(b/v)
    Lb = b == 0 ? inv(v) : L/b
    F = T*aT-1.5a
    u_departure = F*Lb/sqt
    h_departure = u_departure+P*v-RT
    s_departure = R*log(Z*(1-b/v))+(aT-a/(2T))*Lb/sqt
    cv = R*(dot(X,w.cp0)-1)+Lb/sqt*(0.75a/T-aT)
    dpdT = R/(v-b)-(aT-a/(2T))/(sqt*v*(v+b))
    cp = cv-T*dpdT^2/dpdv
    cv > 0 && cp > 0 || throw(DomainError((T,rho),"RK heat capacities must be positive"))
    h = RT*dot(X,w.h0)+h_departure
    u = h-P*v
    s = R*(dot(X,w.s0)-sum(x == 0 ? zero(x) : x*log(x) for x in X)-log(P/one_atm))+s_departure
    @inbounds for i in 1:n
        bi = m.b[i]
        if b == 0
            w.u_TV[i] = RT*(w.h0[i]-1)
            w.lnphi[i] = 0
            w.hbar[i] = RT*w.h0[i]
            w.vbar[i] = v
            continue
        end
        Sk = 2T*w.dAk[i]-3w.Ak[i]
        C = -Lb+1/(v+b)
        w.u_TV[i] = RT*(w.h0[i]-1)+(L*Sk+bi*F*C)/(b*sqt)
        w.lnphi[i] = log(RT/((v-b)*P))+bi/(v-b)-2w.Ak[i]*Lb/(sqt*RT)+
            a*bi*Lb/(b*sqt*RT)-a*bi/(b*sqt*(v+b)*RT)
        dpdni = RT/(v-b)*(1+bi/(v-b))-2w.Ak[i]/(sqt*v*(v+b))+a*bi/(sqt*v*(v+b)^2)
        w.vbar[i] = -dpdni/dpdv
        w.hbar[i] = RT*w.h0[i]-RT-bi*Lb*F/(b*sqt)+Lb*Sk/sqt+
            bi*F/((v+b)*b*sqt)-T*dpdT/dpdv*dpdni
    end
    return (;T,P,rho,X,MW,v,Z,h,u,s,cp,cv,h_mass=h/MW,u_mass=u/MW,s_mass=s/MW,
        cp_mass=cp/MW,cv_mass=cv/MW,h_departure,u_departure,s_departure,dpdT,dpdv,
        compressibility=-1/(v*dpdv),expansion=-dpdT/(v*dpdv),
        sound_speed=v*sqrt(-cp/cv*dpdv/MW),u_TV=w.u_TV,lnphi=w.lnphi,
        partial_molar_enthalpies=w.hbar,partial_molar_volumes=w.vbar)
end

"""
    redlich_kwong_state(model; T, X, P=nothing, rho=nothing, basis=:mole, root=:gas)

Return native RK mixture and partial properties from T/P or T/density. The
largest physical cubic root is selected by default; `root=:liquid` explicitly
selects the smallest. Selection does not perform phase-equilibrium stability or
flash calculations. All arguments use SI units.
"""
function redlich_kwong_state(m::RedlichKwongThermo;T,X,P=nothing,rho=nothing,basis=:mole,root=:gas)
    xor(isnothing(P),isnothing(rho)) || throw(ArgumentError("specify exactly one of P or rho"))
    isfinite(T) && T > 0 || throw(ArgumentError("positive finite temperature required"))
    x = _rk_composition(m,X,Symbol(basis))
    work = RedlichKwongWorkspace(m)
    if !isnothing(P)
        isfinite(P) && P > 0 || throw(ArgumentError("positive finite pressure required"))
        a,_,b = _rk_mixing!(work,m,T,x)
        rho = dot(x,m.MW)/_rk_volume(T,P,a,b,Symbol(root))
    end
    return redlich_kwong_properties!(work,m,T,rho,x)
end

export RedlichKwongThermo,RedlichKwongWorkspace,redlich_kwong_state,redlich_kwong_properties!
