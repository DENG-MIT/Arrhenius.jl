"""
    SurfaceMechanism(path)

Load parameters prepared by `mechanism/export_surface.py`. Supports one ideal
surface, one ideal gas, and persistent fixed-stoichiometry solids; NASA7 and
constant-cp species; elementary and sticking Arrhenius reactions, explicit
orders, Motz–Wise corrections, and polynomial coverage activation energies.
Species are ordered surface, gas, then bulk. Rates use kmol, m, s, K, and J.
"""
struct SurfaceMechanism
    name::String
    species_names::Vector{String}
    element_names::Vector{String}
    n_surface::Int
    n_gas::Int
    site_density::Float64
    site_sizes::Vector{Float64}
    molecular_weights::Vector{Float64}
    elemental_matrix::Matrix{Float64}
    reactants::Matrix{Float64}
    products::Matrix{Float64}
    orders::Matrix{Float64}
    stoichiometry::Matrix{Float64}
    arrhenius::Matrix{Float64}
    reversible::Vector{Bool}
    coverage_a::Matrix{Float64}
    coverage_m::Matrix{Float64}
    coverage_energy::Array{Float64,3}
    sticking_species::Vector{Int}
    sticking_order::Vector{Float64}
    sticking_factor::Vector{Float64}
    motz_wise::Vector{Bool}
    thermo_type::Vector{Int}
    thermo_coefficients::Matrix{Float64}
    reference_pressure::Vector{Float64}
    bulk_molar_volumes::Vector{Float64}
    initial_coverages::Vector{Float64}
    initial_mole_fractions::Vector{Float64}
    initial_temperature::Float64
    initial_pressure::Float64
    gas_file::String
end

function SurfaceMechanism(path::AbstractString)
    data = npzread(path)
    str(key) = String(vec(UInt8.(data[key])))
    str("surface_format_utf8") == "arrhenius-surface-v1" ||
        throw(ArgumentError("unsupported surface parameter format"))
    vecf(key) = vec(Float64.(data[key]))
    mat(key) = Matrix{Float64}(data[key])
    ns, ng = Int(only(data["n_surface"])), Int(only(data["n_gas"]))
    names = String.(split(str("species_names_utf8"), '\n'))
    reactants, products, orders = mat("reactants"), mat("products"), mat("orders")
    nu = products - reactants
    nt, nr = size(nu)
    ns > 0 && ng > 0 && ns + ng <= nt && length(names) == nt ||
        throw(ArgumentError("invalid surface/gas species dimensions"))
    length(unique(names)) == nt || throw(ArgumentError("duplicate species names"))
    size(orders) == size(nu) && size(reactants) == size(nu) ||
        throw(ArgumentError("inconsistent reaction matrices"))
    all(isfinite, orders) && all(>=(0), orders) || throw(ArgumentError("invalid reaction orders"))
    site_density = only(vecf("site_density"))
    sizes = vecf("site_sizes")
    isfinite(site_density) && site_density > 0 && length(sizes) == ns && all(>(0), sizes) ||
        throw(ArgumentError("positive site density and species sizes required"))
    elements = mat("elemental_matrix")
    maximum(abs, transpose(sizes) * nu[1:ns, :]; init=0.0) < 1e-10 ||
        throw(ArgumentError("reaction site balance is not conserved"))
    maximum(abs, elements * nu; init=0.0) < 1e-10 ||
        throw(ArgumentError("reaction elemental balance is not conserved"))
    types = vec(Int.(data["thermo_type"]))
    all(t -> t in (1, 2), types) || throw(ArgumentError("unsupported surface species thermo"))
    gas_file = haskey(data, "gas_file_utf8") ? str("gas_file_utf8") : ""
    return SurfaceMechanism(str("phase_name_utf8"), names,
        String.(split(str("element_names_utf8"), '\n')), ns, ng, site_density, sizes,
        vecf("molecular_weights"), elements, reactants, products, orders, nu,
        mat("arrhenius"), vec(Bool.(data["reversible"])), mat("coverage_a"),
        mat("coverage_m"), Array{Float64,3}(data["coverage_energy"]),
        vec(Int.(data["sticking_species"])), vecf("sticking_order"),
        vecf("sticking_factor"), vec(Bool.(data["motz_wise"])), types,
        mat("thermo_coefficients"), vecf("reference_pressure"),
        haskey(data, "bulk_molar_volumes") ? vecf("bulk_molar_volumes") : Float64[],
        vecf("initial_coverages"), vecf("initial_mole_fractions"),
        only(vecf("initial_temperature")), only(vecf("initial_pressure")),
        isempty(gas_file) ? "" : joinpath(dirname(abspath(path)), gas_file))
end

function _surface_composition(names, value)
    x = zeros(length(names))
    if value isa AbstractDict || value isa NamedTuple
        for (name, amount) in pairs(value)
            k = findfirst(==(String(name)), names)
            isnothing(k) && throw(ArgumentError("unknown species $name"))
            x[k] = amount
        end
    else
        length(value) == length(names) || throw(DimensionMismatch("composition length"))
        x .= value
    end
    all(isfinite, x) && all(>=(0), x) && sum(x) > 0 ||
        throw(ArgumentError("composition must be nonnegative, finite, and nonzero"))
    return x ./ sum(x)
end

"""
    IdealSurface(mechanism; temperature, pressure, mole_fractions, coverages,
                 gas_temperature=temperature)

A surface with fixed temperature, gas pressure, and adjacent gas composition.
Coverages denote fractions of occupied sites and sum to one. An optional gas
temperature permits direct nonisothermal interface rate evaluation. For a
stationary surface equilibrated thermally with its gas, use the default equal
temperatures. Persistent solid phases have unit activity and positive inventory.
"""
struct IdealSurface
    mechanism::SurfaceMechanism
    temperature::Float64
    pressure::Float64
    gas_temperature::Float64
    mole_fractions::Vector{Float64}
    coverages::Vector{Float64}
end
function IdealSurface(mechanism::SurfaceMechanism;
                      temperature=mechanism.initial_temperature,
                      pressure=mechanism.initial_pressure, gas_temperature=temperature,
                      mole_fractions=mechanism.initial_mole_fractions,
                      coverages=mechanism.initial_coverages)
    all(x -> isfinite(x) && x > 0, (temperature, pressure, gas_temperature)) ||
        throw(ArgumentError("positive finite temperatures and pressure required"))
    ns, ng = mechanism.n_surface, mechanism.n_gas
    X = _surface_composition(mechanism.species_names[ns+1:ns+ng], mole_fractions)
    theta = _surface_composition(mechanism.species_names[1:ns], coverages)
    return IdealSurface(mechanism, temperature, pressure, gas_temperature, X, theta)
end

"A reusable surface rate workspace. Returned rate arrays alias this workspace."
struct SurfaceWorkspace
    concentrations::Vector{Float64}
    standard_mu_RT::Vector{Float64}
    enthalpy::Vector{Float64}
    forward_rate_constants::Vector{Float64}
    reverse_rate_constants::Vector{Float64}
    forward_rates::Vector{Float64}
    reverse_rates::Vector{Float64}
    net_rates::Vector{Float64}
    production_rates::Vector{Float64}
    coverage_rates::Vector{Float64}
end
function SurfaceWorkspace(m::SurfaceMechanism)
    nt, nr = size(m.stoichiometry)
    return SurfaceWorkspace(zeros(nt), zeros(nt), zeros(nt), zeros(nr), zeros(nr),
        zeros(nr), zeros(nr), zeros(nr), zeros(nt), zeros(m.n_surface))
end

# Standard enthalpy and entropy at the species reference pressure.
@inline function _surface_species_thermo(m, k, T)
    c = m.thermo_coefficients
    if m.thermo_type[k] == 1
        j = T <= c[k, 1] ? 9 : 2
        a1, a2, a3, a4, a5, a6, a7 = (c[k, j+i] for i in 0:6)
        h = R * (T*(a1 + T*(a2/2 + T*(a3/3 + T*(a4/4 + T*a5/5)))) + a6)
        s = R * (a1*log(T) + T*(a2 + T*(a3/2 + T*(a4/3 + T*a5/4))) + a7)
        return h, s
    end
    T0, h0, s0, cp0 = c[k, 1], c[k, 2], c[k, 3], c[k, 4]
    return h0 + cp0*(T-T0), s0 + cp0*log(T/T0)
end

"""
    surface_rates!(workspace, mechanism, T, P, X, coverages; gas_temperature=T)
    surface_rates!(workspace, surface, coverages=surface.coverages)

Calculate elementary reaction and species production rates per unit surface
area. `workspace.production_rates` is in kmol/m²/s in mechanism species order;
`coverage_rates` is in s⁻¹ and includes each species' site size. `X` and
`coverages` must already be normalized. Tiny negative trial coverages are
clipped only while evaluating rates; accepted ODE states are domain checked.
"""
function surface_rates!(w::SurfaceWorkspace, m::SurfaceMechanism, T, P, X, theta;
                        gas_temperature=T)
    ns, ng = m.n_surface, m.n_gas
    length(X) == ng && length(theta) == ns || throw(DimensionMismatch("surface composition"))
    isfinite(T) && isfinite(P) && isfinite(gas_temperature) && min(T, P, gas_temperature) > 0 ||
        throw(DomainError((T, P, gas_temperature), "positive finite surface state required"))
    nt, nr = size(m.stoichiometry)
    for k in 1:nt
        Tk = ns < k <= ns+ng ? gas_temperature : T
        h, s = _surface_species_thermo(m, k, Tk)
        mu = h - Tk*s
        c0 = 1.0
        if k <= ns
            c0 = m.site_density / m.site_sizes[k]
            w.concentrations[k] = c0 * max(theta[k], 0.0)
        elseif k <= ns+ng
            c0 = P / (R*Tk)
            mu += R*Tk*log(P/m.reference_pressure[k])
            w.concentrations[k] = c0 * max(X[k-ns], 0.0)
        else
            h += (P-m.reference_pressure[k]) * m.bulk_molar_volumes[k-ns-ng]
            mu = h - Tk*s
            w.concentrations[k] = 1.0
        end
        w.enthalpy[k] = h
        w.standard_mu_RT[k] = mu/(R*T) - log(c0)
    end
    logT = log(T)
    for j in 1:nr
        correction = 0.0
        for k in 1:ns
            q = max(theta[k], 0.0)
            energy = q*(m.coverage_energy[j,k,1] + q*(m.coverage_energy[j,k,2] +
                     q*(m.coverage_energy[j,k,3] + q*m.coverage_energy[j,k,4])))
            correction += log(10.0)*m.coverage_a[j,k]*q - energy/(R*T)
            if m.coverage_m[j,k] != 0
                correction += m.coverage_m[j,k]*log(max(q, 1e-300))
            end
        end
        kf = m.arrhenius[j,1] * exp(m.arrhenius[j,2]*logT - m.arrhenius[j,3]/(R*T) + correction)
        if m.sticking_species[j] > 0
            m.motz_wise[j] && (kf /= 1 - 0.5*kf)
            kf *= m.site_density^(-m.sticking_order[j]) * sqrt(T) * m.sticking_factor[j]
        end
        delta_mu = 0.0
        for k in 1:nt
            delta_mu += m.stoichiometry[k,j] * w.standard_mu_RT[k]
        end
        kr = m.reversible[j] ? kf * exp(clamp(delta_mu, -700.0, 700.0)) : 0.0
        w.forward_rate_constants[j], w.reverse_rate_constants[j] = kf, kr
        qf, qr = kf, kr
        for k in 1:nt
            m.orders[k,j] != 0 && (qf *= w.concentrations[k]^m.orders[k,j])
            m.reversible[j] && m.products[k,j] != 0 &&
                (qr *= w.concentrations[k]^m.products[k,j])
        end
        w.forward_rates[j], w.reverse_rates[j] = qf, qr
        w.net_rates[j] = qf - qr
    end
    mul!(w.production_rates, m.stoichiometry, w.net_rates)
    for k in 1:ns
        w.coverage_rates[k] = m.site_sizes[k]*w.production_rates[k]/m.site_density
    end
    return w
end
surface_rates!(w::SurfaceWorkspace, s::IdealSurface, theta=s.coverages) =
    surface_rates!(w, s.mechanism, s.temperature, s.pressure, s.mole_fractions, theta;
                   gas_temperature=s.gas_temperature)

"Calculate an independent copy of gas/surface/bulk production and coverage rates."
function surface_rates(s::IdealSurface, theta=s.coverages)
    w = surface_rates!(SurfaceWorkspace(s.mechanism), s, theta)
    ns, ng = s.mechanism.n_surface, s.mechanism.n_gas
    return (surface=w.production_rates[1:ns], gas=w.production_rates[ns+1:ns+ng],
        bulk=w.production_rates[ns+ng+1:end], coverages=w.coverage_rates,
        forward=w.forward_rates, reverse=w.reverse_rates, net=w.net_rates,
        forward_rate_constants=w.forward_rate_constants,
        reverse_rate_constants=w.reverse_rate_constants,
        elemental_rates=s.mechanism.elemental_matrix*w.production_rates,
        heat_release=-dot(w.enthalpy, w.production_rates))
end

struct SurfaceRHS
    surface::IdealSurface
    workspace::SurfaceWorkspace
    trial::Vector{Float64}
    plus::Vector{Float64}
    minus::Vector{Float64}
end
surface_rhs(s::IdealSurface) = SurfaceRHS(s, SurfaceWorkspace(s.mechanism),
    copy(s.coverages), similar(s.coverages), similar(s.coverages))
function (rhs::SurfaceRHS)(du, u, p, t)
    surface_rates!(rhs.workspace, rhs.surface, u)
    copyto!(du, rhs.workspace.coverage_rates)
    return nothing
end

"Second-order finite-difference coverage Jacobian, including rate dependencies."
function surface_jacobian!(J, u, rhs::SurfaceRHS, t=0.0)
    n = length(u)
    size(J) == (n, n) || throw(DimensionMismatch("surface Jacobian"))
    copyto!(rhs.trial, u)
    for j in 1:n
        h = cbrt(eps(Float64))*max(abs(u[j]), 1e-7)
        rhs.trial[j] = u[j] + h
        rhs(rhs.plus, rhs.trial, nothing, t)
        if u[j] > h
            rhs.trial[j] = u[j] - h
            rhs(rhs.minus, rhs.trial, nothing, t)
            for i in 1:n
                J[i,j] = (rhs.plus[i] - rhs.minus[i])/(2h)
            end
        else
            rhs.trial[j] = u[j] + 2h
            rhs(rhs.minus, rhs.trial, nothing, t)
            surface_rates!(rhs.workspace, rhs.surface, u)
            for i in 1:n
                J[i,j] = (-3rhs.workspace.coverage_rates[i] + 4rhs.plus[i] - rhs.minus[i])/(2h)
            end
        end
        rhs.trial[j] = u[j]
    end
    return J
end

function surface_problem(s::IdealSurface, tspan; initial_coverages=s.coverages)
    length(tspan) == 2 && all(isfinite, tspan) && tspan[2] > tspan[1] ||
        throw(ArgumentError("finite increasing integration interval required"))
    rhs = surface_rhs(s)
    theta = _surface_restart(s, initial_coverages)
    return (f=rhs, jac=(J,u,p,t)->surface_jacobian!(J,u,rhs,t),
        tgrad=(du,u,p,t)->(fill!(du, 0.0); nothing), u0=theta,
        tspan=(float(tspan[1]), float(tspan[2])), p=nothing,
        isoutofdomain=(u,p,t)->any(x->!isfinite(x) || x < -1e-13, u))
end

function _surface_restart(s, theta)
    if theta isa AbstractVector
        all(x -> isfinite(x) && x >= -1e-13, theta) ||
            throw(ArgumentError("invalid restart coverages"))
        theta = max.(theta, 0.0)
    end
    return _surface_composition(s.mechanism.species_names[1:s.mechanism.n_surface], theta)
end

"Integrate surface coverages with a caller-owned stiff ODE integrator."
function solve_surface(s::IdealSurface, tspan; integrator, initial_coverages=s.coverages, kwargs...)
    return integrator(surface_problem(s, tspan; initial_coverages); kwargs...)
end

"""
    steady_coverages(surface; integrator, interval=1, max_time=1e5,
                     steady_tolerance=1e-6, kwargs...)

Integrate fixed-gas coverage equations until the maximum absolute coverage
derivative is below `steady_tolerance` in s⁻¹. Returns the coverages, residual,
elapsed physical time, and final integration solution. Temperatures and
pressure are exactly those supplied to `IdealSurface` throughout the solve.
At a restart, negative roundoff coverages within 1e-13 are set to zero and
the site fractions are normalized before continuing.
"""
function steady_coverages(s::IdealSurface; integrator, interval=1.0, max_time=1e5,
                          steady_tolerance=1e-6, initial_coverages=s.coverages, kwargs...)
    interval > 0 && max_time > 0 && steady_tolerance > 0 ||
        throw(ArgumentError("positive steady integration controls required"))
    u = _surface_restart(s, initial_coverages)
    time = 0.0
    rhs = surface_rhs(s)
    derivative = similar(u)
    solution = nothing
    while time < max_time
        stop = min(time + interval, max_time)
        solution = solve_surface(s, (time, stop); integrator, initial_coverages=u, kwargs...)
        u = _surface_restart(s, solution.u[end])
        time = stop
        rhs(derivative, u, nothing, time)
        residual = maximum(abs, derivative)
        residual <= steady_tolerance && return (coverages=u, residual, time, solution)
        interval = min(2interval, max_time-time)
    end
    error("surface did not reach the requested stationary residual by t=$max_time")
end

export SurfaceMechanism, IdealSurface, SurfaceWorkspace, SurfaceRHS
export surface_rates!, surface_rates, surface_rhs, surface_jacobian!, surface_problem
export solve_surface, steady_coverages

"""
    ReactorSurface(node, mechanism; area, coverages=mechanism.initial_coverages)

Attach a fixed-area ideal surface to an isothermal `WellStirredReactor` named
`node`. The surface follows that reactor's temperature, pressure, and gas
composition. Persistent bulk products are tracked as signed deposited amounts.
"""
struct ReactorSurface
    node::Symbol
    mechanism::SurfaceMechanism
    area::Float64
    coverages::Vector{Float64}
end
function ReactorSurface(node::Symbol, m::SurfaceMechanism; area,
                        coverages=m.initial_coverages)
    isfinite(area) && area > 0 || throw(ArgumentError("positive finite surface area required"))
    theta = _surface_composition(m.species_names[1:m.n_surface], coverages)
    return ReactorSurface(node, m, area, theta)
end

"""
    CatalyticNetwork(network; surfaces)

Extend a native reactor network with gas/surface species exchange. Reactors
carrying surfaces must be isothermal; other network vessels retain their own
energy settings. The state contains the network state followed by each surface's
coverages and signed bulk deposition in kmol. No solid depletion or surface
energy equation is assumed. `catalytic_diagnostics` reports combined gas, surface,
and deposited-bulk mass and elemental inventories and their rates.
"""
struct CatalyticNetwork{N,S}
    network::N
    surfaces::S
    node_indices::Vector{Int}
    offsets::Vector{Int}
    gas_maps::Vector{Vector{Int}}
    initial_state::Vector{Float64}
end
function CatalyticNetwork(network::ReactorNetwork; surfaces)
    surfaces = Tuple(surfaces)
    isempty(surfaces) && throw(ArgumentError("at least one reactor surface required"))
    indices, offsets, maps = Int[], Int[], Vector{Int}[]
    u = network_state(network)
    for surface in surfaces
        i = findfirst(==(surface.node), keys(network.nodes))
        isnothing(i) && throw(ArgumentError("unknown surface node $(surface.node)"))
        node = network.nodes[i]
        node isa WellStirredReactor && node.initial.energy === :isothermal ||
            throw(ArgumentError("reactors with attached surfaces must be isothermal"))
        m, gas = surface.mechanism, node.initial.gas
        names = m.species_names[m.n_surface+1:m.n_surface+m.n_gas]
        length(names) == gas.n_species && Set(names) == Set(gas.species_names) ||
            throw(ArgumentError("reactor gas must match the surface's companion gas species"))
        mapping = [findfirst(==(name), gas.species_names)::Int for name in names]
        for (k,j) in enumerate(mapping)
            isapprox(gas.MW[j],m.molecular_weights[m.n_surface+k];rtol=1e-10) ||
                throw(ArgumentError("surface and gas molecular weights differ"))
            for name in union(gas.elements,m.element_names)
                a,b = findfirst(==(name),gas.elements),findfirst(==(name),m.element_names)
                gas_atoms = isnothing(a) ? 0.0 : gas.ele_matrix[a,j]
                surface_atoms = isnothing(b) ? 0.0 : m.elemental_matrix[b,m.n_surface+k]
                gas_atoms == surface_atoms || throw(ArgumentError("surface and gas elemental compositions differ"))
            end
        end
        push!(indices,i); push!(offsets,length(u)+1); push!(maps,mapping)
        append!(u,surface.coverages)
        append!(u,zeros(length(m.bulk_molar_volumes)))
    end
    return CatalyticNetwork(network,surfaces,indices,offsets,maps,u)
end
catalytic_state(n::CatalyticNetwork) = copy(n.initial_state)

struct CatalyticRHS{N,R}
    system::N
    network_rhs::R
    workspaces::Vector{SurfaceWorkspace}
    mole_fractions::Vector{Vector{Float64}}
    trial::Vector{Float64}
    plus::Vector{Float64}
    minus::Vector{Float64}
    base::Vector{Float64}
end
function catalytic_rhs(n::CatalyticNetwork)
    return CatalyticRHS(n,network_rhs(n.network),
        [SurfaceWorkspace(s.mechanism) for s in n.surfaces],
        [zeros(s.mechanism.n_gas) for s in n.surfaces],catalytic_state(n),
        zeros(length(n.initial_state)),zeros(length(n.initial_state)),zeros(length(n.initial_state)))
end
function (rhs::CatalyticRHS)(du,u,p,t)
    n = rhs.system
    length(u) == length(du) == length(n.initial_state) || throw(DimensionMismatch("catalytic state"))
    nn = length(n.network.initial_state)
    fill!(du,0.0)
    rhs.network_rhs(@view(du[1:nn]),@view(u[1:nn]),p,t)
    for (a,surface) in enumerate(n.surfaces)
        m, i, offset = surface.mechanism,n.node_indices[a],n.offsets[a]
        state = rhs.network_rhs.states[i]
        gas = n.network.nodes[i].initial.gas
        X = rhs.mole_fractions[a]
        for (k,j) in enumerate(n.gas_maps[a])
            X[k] = rhs.network_rhs.workspaces[i].X[j]
        end
        w = surface_rates!(rhs.workspaces[a],m,state.temperature,state.pressure,X,
                           @view(u[offset:offset+m.n_surface-1]))
        gasoffset = n.network.offsets[i]
        for (k,j) in enumerate(n.gas_maps[a])
            du[gasoffset+j-1] += surface.area*w.production_rates[m.n_surface+k]*gas.MW[j]
        end
        for k in 1:m.n_surface
            du[offset+k-1] = w.coverage_rates[k]
        end
        for k in eachindex(m.bulk_molar_volumes)
            du[offset+m.n_surface+k-1] = surface.area*w.production_rates[m.n_surface+m.n_gas+k]
        end
    end
    return nothing
end

function catalytic_jacobian!(J,u,rhs::CatalyticRHS,t=0.0)
    n = rhs.system
    length(u) == size(J,1) == size(J,2) || throw(DimensionMismatch("catalytic Jacobian"))
    copyto!(rhs.trial,u)
    rhs(rhs.base,u,nothing,t)
    nn = length(n.network.initial_state)
    for j in eachindex(u)
        scale = 1e-7
        if j <= nn
            for (i,node) in enumerate(n.network.nodes)
                offset = n.network.offsets[i]
                offset == 0 && continue
                if offset <= j < offset+node.initial.gas.n_species
                    scale = 1e-7*rhs.network_rhs.states[i].mass
                    break
                end
            end
        end
        h = cbrt(eps(Float64))*max(abs(u[j]),scale)
        rhs.trial[j] = u[j]+h
        rhs(rhs.plus,rhs.trial,nothing,t)
        if u[j] > h
            rhs.trial[j] = u[j]-h
            rhs(rhs.minus,rhs.trial,nothing,t)
            for i in eachindex(u)
                J[i,j] = (rhs.plus[i]-rhs.minus[i])/(2h)
            end
        else
            rhs.trial[j] = u[j]+2h
            rhs(rhs.minus,rhs.trial,nothing,t)
            for i in eachindex(u)
                J[i,j] = (-3rhs.base[i]+4rhs.plus[i]-rhs.minus[i])/(2h)
            end
        end
        rhs.trial[j] = u[j]
    end
    return J
end

function catalytic_problem(n::CatalyticNetwork,tspan;initial_state=catalytic_state(n))
    length(tspan) == 2 && all(isfinite,tspan) && tspan[2] > tspan[1] ||
        throw(ArgumentError("finite increasing integration interval required"))
    length(initial_state) == length(n.initial_state) || throw(DimensionMismatch("catalytic state"))
    rhs = catalytic_rhs(n)
    nn = length(n.network.initial_state)
    function isoutofdomain(u,p,t)
        all(isfinite,u) || return true
        network_isoutofdomain(n.network,@view(u[1:nn])) && return true
        for (a,s) in enumerate(n.surfaces)
            any(<(-1e-13),@view(u[n.offsets[a]:n.offsets[a]+s.mechanism.n_surface-1])) && return true
        end
        return false
    end
    function tgrad(du,u,p,t)
        h = sqrt(eps(Float64))*max(abs(t),1.0)
        rhs(rhs.plus,u,p,t+h)
        rhs(rhs.minus,u,p,t)
        @. du = (rhs.plus-rhs.minus)/h
        return nothing
    end
    return (f=rhs,jac=(J,u,p,t)->catalytic_jacobian!(J,u,rhs,t),tgrad=tgrad,
        u0=Float64.(initial_state),tspan=(float(tspan[1]),float(tspan[2])),p=nothing,isoutofdomain=isoutofdomain)
end
function solve_catalytic(n::CatalyticNetwork,tspan;integrator,initial_state=catalytic_state(n),kwargs...)
    return integrator(catalytic_problem(n,tspan;initial_state);kwargs...)
end

"Combined inventory and conservation rates, including surface storage and bulk deposition."
function catalytic_diagnostics(rhs::CatalyticRHS,u,t=0.0)
    n = rhs.system
    du = similar(u)
    rhs(du,u,nothing,t)
    names = sort!(unique(vcat([String.(node.initial.gas.elements) for node in n.network.nodes]...,
                              [s.mechanism.element_names for s in n.surfaces]...)))
    inventory, rates = zeros(length(names)),zeros(length(names))
    mass, mass_rate, energy, energy_rate = 0.0,0.0,0.0,0.0
    function add_elements!(elements,coefficient,amount,rate)
        for (e,name) in enumerate(elements)
            j = findfirst(==(String(name)),names)
            inventory[j] += coefficient[e]*amount
            rates[j] += coefficient[e]*rate
        end
    end
    for (i,node) in enumerate(n.network.nodes)
        offset = n.network.offsets[i]
        offset == 0 && continue
        gas,state,w = node.initial.gas,rhs.network_rhs.states[i],rhs.network_rhs.workspaces[i]
        mass += state.mass
        energy += state.mass*state.internal_energy
        energy_rate += state.mass*state.cv*du[offset+gas.n_species]
        for k in 1:gas.n_species
            mk, dmk = u[offset+k-1],du[offset+k-1]
            mass_rate += dmk
            energy_rate += (w.h_mole[k]-R*state.temperature)/gas.MW[k]*dmk
            add_elements!(gas.elements,@view(gas.ele_matrix[:,k]),mk/gas.MW[k],dmk/gas.MW[k])
        end
    end
    for (a,s) in enumerate(n.surfaces)
        m,offset,w = s.mechanism,n.offsets[a],rhs.workspaces[a]
        P = rhs.network_rhs.states[n.node_indices[a]].pressure
        for k in 1:m.n_surface
            amount = s.area*m.site_density*u[offset+k-1]/m.site_sizes[k]
            rate = s.area*m.site_density*du[offset+k-1]/m.site_sizes[k]
            mass += amount*m.molecular_weights[k]
            mass_rate += rate*m.molecular_weights[k]
            energy += amount*w.enthalpy[k]
            energy_rate += rate*w.enthalpy[k]
            add_elements!(m.element_names,@view(m.elemental_matrix[:,k]),amount,rate)
        end
        for b in eachindex(m.bulk_molar_volumes)
            k = m.n_surface+m.n_gas+b
            amount,rate = u[offset+m.n_surface+b-1],du[offset+m.n_surface+b-1]
            uk = w.enthalpy[k]-P*m.bulk_molar_volumes[b]
            mass += amount*m.molecular_weights[k]
            mass_rate += rate*m.molecular_weights[k]
            energy += amount*uk
            energy_rate += rate*uk
            add_elements!(m.element_names,@view(m.elemental_matrix[:,k]),amount,rate)
        end
    end
    # Flow enthalpy and fixed-wall energy already include internal cancellation.
    boundary_power = sum(rhs.network_rhs.energy_flux)
    return (mass=mass,mass_rate=mass_rate,element_names=names,element_inventory=inventory,
        element_rates=rates,internal_energy=energy,internal_energy_rate=energy_rate,
        boundary_power=boundary_power,thermostat_power=energy_rate-boundary_power,
        derivative=du)
end
catalytic_diagnostics(n::CatalyticNetwork,u,t=0.0) = catalytic_diagnostics(catalytic_rhs(n),u,t)

export ReactorSurface, CatalyticNetwork, CatalyticRHS, catalytic_state, catalytic_rhs
export catalytic_jacobian!, catalytic_problem, solve_catalytic, catalytic_diagnostics
