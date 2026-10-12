struct IonizedFlameWorkspace{K,C}
    kinetics::K
    transports::Vector{IonTransportWorkspace}
    X::Matrix{Float64}
    rho::Vector{Float64}
    cp::Matrix{Float64}
    h::Matrix{Float64}
    source::Matrix{Float64}
    flux::Matrix{Float64}
    conductivity::Vector{Float64}
    C::Vector{Float64}
    entropy::Vector{Float64}
    xmid::Vector{Float64}
    ymid::Vector{Float64}
    gradient::Vector{Float64}
    band::Matrix{Float64}
    perturbed::Matrix{Float64}
    residual_perturbed::Matrix{Float64}
    steps::Vector{Float64}
    rate_caches::Vector{_KineticsTemperatureCache}
    excess::Vector{Int}
    conservative::C
end

function FlameWorkspace(f::IonizedFlame)
    n, N = f.gas.n_species, length(f.grid)
    B = n + 3
    return IonizedFlameWorkspace(
        KineticsWorkspace(f.gas.reaction),
        [IonTransportWorkspace(f.ion_data) for _ in 1:N-1],
        zeros(n, N),
        zeros(N),
        zeros(n, N),
        zeros(n, N),
        zeros(n, N),
        zeros(n, N - 1),
        zeros(N - 1),
        zeros(n),
        zeros(n),
        zeros(n),
        zeros(n),
        zeros(n),
        zeros(6B - 2, B * N),
        zeros(B, N),
        zeros(B, N),
        zeros(N),
        [_KineticsTemperatureCache(f.gas.reaction) for _ in 1:N],zeros(Int,2),
        f.discretization==:conservative ? IonizedConservativeWorkspace(f.ion_data,N) : nothing)
end

function _ion_flame_properties!(w::IonizedFlameWorkspace, f::IonizedFlame, u;
        update_transport=true, nodes=eachindex(f.grid))
    gas, n, N = f.gas, f.gas.n_species, length(f.grid)
    B = n + 3

    size(u) == (B, N) || throw(DimensionMismatch("ionized flame state must have size (n+3, N)"))
    length(f.ion_data.species_names) == n || throw(DimensionMismatch("ion transport species count mismatch"))
    length(w.transports) == N - 1 || throw(DimensionMismatch("one ion transport workspace per face is required"))
    size(w.X) == (n, N) || throw(DimensionMismatch("workspace X must have size (n, N)"))
    length(w.rho) == N || throw(DimensionMismatch("workspace rho must have length N"))
    size(w.cp) == (n, N) || throw(DimensionMismatch("workspace cp must have size (n, N)"))
    size(w.h) == (n, N) || throw(DimensionMismatch("workspace h must have size (n, N)"))
    size(w.source) == (n, N) || throw(DimensionMismatch("workspace source must have size (n, N)"))
    size(w.flux) == (n, N - 1) || throw(DimensionMismatch("workspace flux must have size (n, N-1)"))
    length(w.conductivity) == N - 1 || throw(DimensionMismatch("workspace conductivity must have length N-1"))
    length(w.C) == n || throw(DimensionMismatch("workspace C must have length n"))
    length(w.entropy) == n || throw(DimensionMismatch("workspace entropy must have length n"))
    length(w.xmid) == n || throw(DimensionMismatch("workspace xmid must have length n"))
    length(w.ymid) == n || throw(DimensionMismatch("workspace ymid must have length n"))
    length(w.gradient) == n || throw(DimensionMismatch("workspace gradient must have length n"))
    size(w.band) == (6B - 2, B * N) || throw(DimensionMismatch("workspace band has invalid size"))
    size(w.perturbed) == (B, N) || throw(DimensionMismatch("workspace perturbed must have size (n+3, N)"))
    size(w.residual_perturbed) == (B, N) || throw(DimensionMismatch("workspace residual_perturbed must have size (n+3, N)"))
    length(w.steps) == N || throw(DimensionMismatch("workspace steps must have length N"))
    length(w.rate_caches) == N || throw(DimensionMismatch("one kinetics temperature cache per node is required"))
    (isfinite(f.pressure) && f.pressure > 0) || throw(ArgumentError("flame pressure must be positive and finite"))

    MW = gas.MW
    if update_transport
        w.excess[1]=argmax(@view(u[2:n+1,1]))
        w.excess[2]=argmax(@view(u[2:n+1,N]))
    end

    @inbounds for j in nodes
        T = 1000 * u[1, j]
        inverseMW = 0.0
        for k in 1:n
            y = u[k + 1, j]
            isfinite(y) || throw(ArgumentError("species state must be finite"))
            inverseMW += y / MW[k]
        end
        (isfinite(T) && T > 0) || throw(DomainError(T, "node temperature must be positive and finite"))
        (isfinite(inverseMW) && inverseMW > 0) ||
            throw(DomainError(inverseMW, "node inverse mean molecular weight must be positive and finite"))

        meanMW = 1 / inverseMW
        w.rho[j] = f.pressure * meanMW / (R * T)
        for k in 1:n
            w.X[k, j] = (u[k + 1, j] / MW[k]) * meanMW
            w.C[k] = f.pressure / (R * T) * w.X[k, j]
        end

        x = @view w.X[:, j]
        h = @view w.h[:, j]
        cp = @view w.cp[:, j]
        cal_h_RT!(h, gas, T, f.pressure, x)
        cal_cp_R!(cp, gas, T, f.pressure, x)
        cal_s0_R!(w.entropy, gas, T, f.pressure, x)
        h .*= R * T
        cp .*= R
        w.entropy .*= R
        wdot!(@view(w.source[:, j]), gas.reaction, T, w.C, w.entropy, h, w.kinetics;
            temperature_cache=w.rate_caches[j])
    end

    frozen = !f.field_enabled
    @inbounds for j in 1:N-1
        dz = f.grid[j + 1] - f.grid[j]
        (isfinite(dz) && dz > 0) || throw(ArgumentError("flame grid spacing must be positive and finite"))

        Tmid = 500 * (u[1, j] + u[1, j + 1])
        inverseMW = 0.0
        for k in 1:n
            ymid = 0.5 * (u[k + 1, j] + u[k + 1, j + 1])
            w.ymid[k] = ymid
            isfinite(ymid) || throw(ArgumentError("midpoint species state must be finite"))
            inverseMW += ymid / MW[k]
        end
        (isfinite(Tmid) && Tmid > 0) ||
            throw(DomainError(Tmid, "face temperature must be positive and finite"))
        (isfinite(inverseMW) && inverseMW > 0) ||
            throw(DomainError(inverseMW, "face inverse mean molecular weight must be positive and finite"))

        meanMW = 1 / inverseMW
        for k in 1:n
            w.xmid[k] = (w.ymid[k] / MW[k]) * meanMW
        end

        if update_transport
            _, lambda, _ = ionized_transport!(w.transports[j], f.ion_data, f.pressure, Tmid, w.xmid;
                mean_molecular_weight=meanMW)
            w.conductivity[j] = lambda
            if w.conservative !== nothing
                rhomid=f.pressure*meanMW/(R*Tmid)
                for k in 1:n
                    w.conservative.density_diffusion[k,j]=rhomid*w.transports[j].diffusion[k]
                end
            end
        end

        for k in 1:n
            w.gradient[k] = (w.X[k, j + 1] - w.X[k, j]) / dz
        end
        ionized_flux!(@view(w.flux[:, j]), w.transports[j], f.ion_data, w.gradient,
            @view(u[2:n+1, j]), @view(u[2:n+1, j + 1]);
            density=w.rho[j], electric_field=1000 * u[n + 2, j], frozen=frozen)
    end
    if update_transport && w.conservative !== nothing
        ionized_transport!(w.conservative.outlet_transport,f.ion_data,f.pressure,
            1000*u[1,N],@view(w.X[:,N]);
            mean_molecular_weight=inv(sum(u[k+1,N]/MW[k] for k in 1:n)))
    end
    return w
end
