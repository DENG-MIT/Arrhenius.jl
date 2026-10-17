using Arrhenius, LinearAlgebra
BLAS.set_num_threads(1)

if length(ARGS) >= 3
    ref = Arrhenius.npzread(ARGS[1])
    mechanism = ARGS[2]
    pressure = parse(Float64, ARGS[3]) * one_atm
    grid = vec(ref["grid"])
    gas = CreateSolution(mechanism)
    Tprofile = vec(ref["T"])
    Yprofile = ref["Y"]
    velocity = vec(ref["velocity"])
    X0 = Dict(name => Yprofile[k,1] for (k,name) in enumerate(gas.species_names)
              if Yprofile[k,1] > 0)
    f = FreeFlame(gas; T=Tprofile[1], P=pressure, X=X0, grid,
        width=grid[end]-grid[1], transport_model=:mixture_averaged,
        discretization=:finite_difference)
    n, N = gas.n_species, length(grid)
    f.state[1,:] .= Tprofile ./ 1000
    f.state[2:n+1,:] .= Yprofile
    for j in 1:N
        rho = Arrhenius._flame_density(f, 1000f.state[1,j], @view(f.state[2:n+1,j])).density
        f.state[end,j] = rho * velocity[j]
    end
else
    root = joinpath(@__DIR__, "..", "test", "data", "cantera")
    ref = Arrhenius.npzread(joinpath(root, "free_auto_stagni.npz"))
    mechanism = joinpath(@__DIR__, "..", "mechanism", "stagni_2023.yaml")
    grid = vec(ref["grid"])
    gas = CreateSolution(mechanism)
    f = FreeFlame(gas; T=300.0, P=one_atm,
        X="NH3:0.7,H2:0.3,O2:0.5625,N2:2.115", grid,
        width=grid[end]-grid[1], transport_model=:mixture_averaged,
        discretization=:finite_difference)
    n, N = gas.n_species, length(grid)
    f.state[1,:] .= vec(ref["values/T"]) ./ 1000
    for (k,name) in enumerate(gas.species_names)
        f.state[k+1,:] .= vec(ref["values/"*name])
    end
    for j in 1:N
        rho = Arrhenius._flame_density(f, 1000f.state[1,j], @view(f.state[2:n+1,j])).density
        f.state[end,j] = rho * ref["values/velocity"][j]
    end
end
w = Arrhenius.FlameWorkspace(f)
r = flame_residual!(similar(f.state), f, f.state, w)

for (mode, analytic) in ((:finite_difference, false), (:cantera_compatible, true)), local_property_updates in (false, true)
    stats = solver_stats(f)
    before = stats.property_nodes
    Arrhenius._flame_jacobian(f, f.state, w, r; analytic, local_property_updates) # compile/warm
    nodes_per_jacobian = stats.property_nodes - before
    samples = Float64[]
    for _ in 1:15
        push!(samples, @elapsed Arrhenius._flame_jacobian(f, f.state, w, r; analytic, local_property_updates))
    end
    println("mode=$mode local_property_updates=$local_property_updates points=$N property_nodes_per_jacobian=$nodes_per_jacobian median_s=$(sort(samples)[cld(length(samples),2)]) min_s=$(minimum(samples))")
end
