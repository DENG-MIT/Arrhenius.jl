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
    value = median(samples)
    println(label," median_ns=",round(Int,value)," alloc_B=",alloc," samples_ns=",round.(Int,samples))
    return value
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
    if hasproperty(gas.reaction, :reactant_mass_action_terms)
        planned_forward_products!() = begin
            total = 0.0
            plans = gas.reaction.reactant_mass_action_terms
            @inbounds for i in 1:gas.reaction.n_reactions
                total += Arrhenius._mass_action_product(C, reactant_species, reactant_orders,
                    gas.reaction.reactant_stoich_coeffs, i,
                    nzrange(gas.reaction.reactant_orders, i), plans[i])
            end
            total
        end
        planned_reverse_products!() = begin
            total = 0.0
            plans = gas.reaction.product_mass_action_terms
            @inbounds for i in eachindex(gas.reaction.is_reversible)
                gas.reaction.is_reversible[i] || continue
                total += Arrhenius._mass_action_product(C, product_species, product_orders,
                    gas.reaction.product_stoich_coeffs, i,
                    nzrange(gas.reaction.product_stoich_coeffs, i), plans[i])
            end
            total
        end
        measure("forward_concentration_products_preclassified",planned_forward_products!)
        measure("reverse_concentration_products_preclassified",planned_reverse_products!)
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

# Profile the complete finite-difference path on a fixed 9-node, 125-species
# flame. This reuses the same residual and Jacobian entry points as solves while
# separating frozen-property work, residual assembly, and band insertion.
function profile_fd_jacobian()
    spec = YAML.load_file(ARGS[1])
    mechanism = spec["mechanism"]
    gas = Arrhenius.CreateSolution(mechanism * ".yaml")
    Xdict = Dict(String(k) => parse(Float64, v) for (k,v) in split.(split(spec["X"], ","), ":"))
    X = Arrhenius.mole_fractions(gas, Xdict)
    P = Float64(spec["P_atm"] * Arrhenius.one_atm)
    f = Arrhenius.FreeFlame(gas;T=300.0,P,X,width=.03,grid=collect(range(0,.03;length=9)))
    w = Arrhenius.FlameWorkspace(f)
    u = f.state
    base = similar(u)
    Arrhenius.flame_residual!(base,f,u,w;update_transport=true)
    N = size(u,2)
    nodes = 2:3:N
    rp = similar(base)
    perturbed = copy(u)
    properties_color!() = Arrhenius._flame_properties!(w,f,perturbed;
        update_transport=false,nodes)
    residual_color!() = Arrhenius.flame_residual!(rp,f,perturbed,w;
        update_transport=false,nodes)
    properties_ns = measure("fd_properties_and_frozen_flux_color",properties_color!;reps=7,n=20)
    residual_ns = measure("fd_residual_color",residual_color!;reps=7,n=20)
    B = size(u,1)
    band = zeros(6B-2,B*N)
    steps = ones(N)
    insertion!() = begin
        fill!(band,0.0)
        @inbounds for k in 1:B, color in 1:3
            for j in color:3:N
                column=(j-1)*B+k
                for jj in max(1,j-1):min(N,j+1), kk in 1:B
                    row=(jj-1)*B+kk
                    band[4B-1+row-column,column]=rp[kk,jj]-base[kk,jj]
                end
            end
        end
        band
    end
    insertion_ns = measure("fd_band_insertion_all_columns",insertion!;reps=7,n=1)
    full_jacobian!() = Arrhenius._flame_jacobian(f,u,w,base;analytic=false)
    jacobian_ns = measure("fd_jacobian_9node",full_jacobian!;reps=5,n=1)
    println("FD_PROFILE nspecies=",gas.n_species," nreactions=",gas.reaction.n_reactions,
        " nodes=",N," state_columns=",B*N," residual_calls=",3B)
    println("FD_PROFILE residual_assembly_est_ns_per_color=",round(Int,residual_ns-properties_ns),
        " residual_share_est_ns=",round(Int,3B*residual_ns),
        " band_insertion_ns=",round(Int,insertion_ns),
        " remaining_jacobian_ns=",round(Int,jacobian_ns-3B*residual_ns-insertion_ns))
end

if get(ENV,"KERNEL_L0_FDJAC","0") == "1"
    main()
    profile_fd_jacobian()
else
    main()
end
