using Arrhenius, YAML, LinearAlgebra, Statistics, InteractiveUtils, SHA, NPZ, SparseArrays
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

    reactant_species, reactant_orders = rowvals(gas.reaction.reactant_orders), nonzeros(gas.reaction.reactant_orders)
    product_species, product_orders = rowvals(gas.reaction.product_stoich_coeffs), nonzeros(gas.reaction.product_stoich_coeffs)
    forward_products!() = begin
        total = 0.0
        @inbounds for i in 1:gas.reaction.n_reactions
            total += Arrhenius._mass_action_product(C, reactant_species, reactant_orders,
                gas.reaction.reactant_stoich_coeffs, i, nzrange(gas.reaction.reactant_orders, i))
        end
        total
    end
    reverse_products!() = begin
        total = 0.0
        @inbounds for i in eachindex(gas.reaction.is_reversible)
            gas.reaction.is_reversible[i] || continue
            total += Arrhenius._mass_action_product(C, product_species, product_orders,
                gas.reaction.product_stoich_coeffs, i, nzrange(gas.reaction.product_stoich_coeffs, i))
        end
        total
    end
    scatter!() = begin
        fill!(out, 0.0)
        rows, coefficients = rowvals(gas.reaction.vk), nonzeros(gas.reaction.vk)
        @inbounds for i in eachindex(kw.rates_of_progress)
            q = kw.rates_of_progress[i]
            for slot in nzrange(gas.reaction.vk, i)
                out[rows[slot]] = muladd(coefficients[slot], q, out[rows[slot]])
            end
        end
        out
    end
    collider_efficiencies!() = begin
        total = 0.0
        @inbounds for i in gas.reaction.index_three_body
            total += dot(view(gas.reaction.efficiencies_coeffs, :, i), C)
        end
        @inbounds for i in gas.reaction.index_falloff
            total += dot(view(gas.reaction.efficiencies_coeffs, :, i), C)
        end
        total
    end

    println("L0 mechanism=",mechanism," nspecies=",n," nreactions=",gas.reaction.n_reactions," T=",T," P=",P," state_index=",state_index)
    measure("thermo_h_cp_s",thermo!)
    report_measure("rate_factors_cached",rates!)
    report_measure("mass_action_qdot_cached",qdot!)
    report_measure("net_production_cached",net!)
    measure("rate_factors_uncached",uncached!)
    qdot!()
    measure("forward_concentration_products",forward_products!)
    measure("reverse_concentration_products",reverse_products!)
    measure("qdot_to_wdot_scatter",scatter!)
    measure("third_body_and_falloff_efficiency_sums",collider_efficiencies!)
    println("checksums h=",sum(h)," cp=",sum(cp)," s=",sum(s)," kf=",sum(kw.kf)," kr=",sum(kw.kr)," wdot=",sum(out))
    println("output_sha256=",bytes2hex(sha256(reinterpret(UInt8,out))))
    @code_warntype Arrhenius._mass_action!(kw,gas.reaction,C)
end
main()
