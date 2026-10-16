using Arrhenius, YAML, LinearAlgebra, Statistics, InteractiveUtils, SHA, NPZ
BLAS.set_num_threads(1)

function measure(label, f; reps=9, n=1000)
    f(); alloc = @allocated f()
    samples = Float64[]
    for _ in 1:reps
        t0 = time_ns()
        for _ in 1:n; f(); end
        push!(samples,(time_ns()-t0)/n)
    end
    println(label," median_ns=",round(Int,median(samples))," alloc_B=",alloc," samples_ns=",round.(Int,samples))
end

function main()
    spec = YAML.load_file(ARGS[1])
    mechanism = spec["mechanism"]
    gas = Arrhenius.CreateSolution(mechanism * ".yaml")
    Xdict = Dict(String(k) => parse(Float64, v) for (k,v) in split.(split(spec["X"], ","), ":"))
    X = Arrhenius.mole_fractions(gas, Xdict)
    T = 1500.0
    state_index = 0
    if length(ARGS) >= 2 && endswith(ARGS[2], ".npz")
        state = NPZ.npzread(ARGS[2])
        state_index = max(1, cld(length(state["T"]), 2))
        T = Float64(state["T"][state_index])
        Y = Float64.(state["Y"][:,state_index])
        X = Vector{Float64}(Y ./ gas.MW)
        X ./= sum(X)
    elseif length(ARGS) >= 2
        T = parse(Float64, ARGS[2])
    end
    P = Float64(spec["P_atm"] * Arrhenius.one_atm)
    C = Vector{Float64}(X .* (P / (Arrhenius.R*T)))
    n = gas.n_species
    h = zeros(n); cp = zeros(n); s = zeros(n)
    hm = zeros(n); sm = zeros(n)
    thermo!() = Arrhenius.cal_h_cp_s0!(h,cp,s,gas,T,P,X)
    thermo!()
    hm .= h .* (Arrhenius.R*T); sm .= s .* Arrhenius.R
    kw = Arrhenius.KineticsWorkspace(gas.reaction)
    out = zeros(n)
    cache = Arrhenius._KineticsTemperatureCache(gas.reaction)
    rates!() = Arrhenius.wdot!(out,gas.reaction,T,C,sm,hm,kw;temperature_cache=cache,pressure=P,get_rate_constants=true)
    qdot!() = Arrhenius.wdot!(out,gas.reaction,T,C,sm,hm,kw;temperature_cache=cache,pressure=P,get_qdot=true)
    net!() = Arrhenius.wdot!(out,gas.reaction,T,C,sm,hm,kw;temperature_cache=cache,pressure=P)
    rates!(); qdot!(); net!()
    c0 = C[1]
    flip = Ref(false)
    function composition_perturb!(f)
        flip[] = !flip[]
        C[1] = c0 * (flip[] ? 1.0 + 1e-10 : 1.0 - 1e-10)
        f()
    end
    report_measure(label,f) = measure(label, () -> composition_perturb!(f))
    uncached!() = Arrhenius.wdot!(out,gas.reaction,T,C,sm,hm,kw;pressure=P,get_rate_constants=true)

    println("L0 mechanism=",mechanism," nspecies=",n," nreactions=",gas.reaction.n_reactions," T=",T," P=",P," state_index=",state_index)
    measure("thermo_h_cp_s",thermo!)
    report_measure("rate_factors_cached",rates!)
    report_measure("mass_action_qdot_cached",qdot!)
    report_measure("net_production_cached",net!)
    measure("rate_factors_uncached",uncached!)
    println("checksums h=",sum(h)," cp=",sum(cp)," s=",sum(s)," kf=",sum(kw.kf)," kr=",sum(kw.kr)," wdot=",sum(out))
    println("output_sha256=",bytes2hex(sha256(reinterpret(UInt8,out))))
    @code_warntype Arrhenius._mass_action!(kw,gas.reaction,C)
end
main()
