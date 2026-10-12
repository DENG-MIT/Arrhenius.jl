# Mass/enthalpy conservative numerical fluxes for ionized flames.
#
# Reuses the shared `_conservative_flame_fluxes!` helper (with the mass-flux
# dispatch hook reading rho*v at ionized nodes) and adds the charged
# zero-gradient electric drift contribution at the outlet. State rows are
# T/1000, n signed raw Y, E/1000, physical velocity; all algebra stays signed.

struct IonizedConservativeWorkspace
    density_diffusion::Matrix{Float64}          # (n, N-1) rho_mid*D per species face
    species_flux::Matrix{Float64}               # (n, N) kg/m2/s species fluxes
    enthalpy_flux::Vector{Float64}              # (N,)
    enthalpy::Vector{Float64}                   # (N,) mixture mass-specific enthalpy [J/kg]
    cp::Vector{Float64}                         # (N,)
    previous_enthalpy::Vector{Float64}          # (N,)
    previous_species_enthalpy::Vector{Float64}  # (n,) species h/RT scratch
    outlet_transport::IonTransportWorkspace
    outlet_flux::Vector{Float64}                # (n,) charged drift outlet flux
    zero_gradient::Vector{Float64}              # (n,) zero dY/dz at outlet
end

function IonizedConservativeWorkspace(data::IonTransportData, N::Integer)
    n = length(data.species_names)
    return IonizedConservativeWorkspace(
        zeros(n, N - 1),
        zeros(n, N),
        zeros(N),
        zeros(N),
        zeros(N),
        zeros(N),
        zeros(n),
        IonTransportWorkspace(data),
        zeros(n),
        zeros(n),
    )
end

function _ion_conservative_fluxes!(c, f, u, w)
    n = f.gas.n_species
    N = length(f.grid)
    _conservative_flame_fluxes!(c, f, u, w)
    # c.outlet_transport is initialized by the property updater at the
    # endpoint T and current X (frozen during Jacobian assembly); zero-gradient
    # electric drift replaces the diffusive part of the natural outflow.
    Yout = @view(u[2:n+1, N])
    ionized_flux!(c.outlet_flux, c.outlet_transport, f.ion_data,
        c.zero_gradient, Yout, Yout;
        density=w.rho[N],
        electric_field=1000 * u[n+2, N],
        frozen=!f.field_enabled)
    MW = f.gas.MW
    @inbounds for k in 1:n
        c.species_flux[k, N] += c.outlet_flux[k]
        c.enthalpy_flux[N] += w.h[k, N] * c.outlet_flux[k] / MW[k]
    end
    return c
end

function _ion_conservative_cell!(res, f, u, w, j; previous=nothing, dt=Inf,
        previous_enthalpy=nothing)
    c = w.conservative
    n = f.gas.n_species
    N = length(f.grid)
    z = f.grid
    MW = f.gas.MW
    cell = (z[min(j + 1, N)] - z[j - 1]) / 2
    factor = _flame_timescale / w.rho[j]
    @inbounds for k in 1:n
        res[k+1, j] = factor * (MW[k] * w.source[k, j] -
            (c.species_flux[k, j] - c.species_flux[k, j-1]) / cell)
    end
    res[1, j] = -factor / (1000 * c.cp[j]) *
        (c.enthalpy_flux[j] - c.enthalpy_flux[j-1]) / cell
    if previous !== nothing
        accumulation = _flame_timescale / dt
        @inbounds for k in 1:n
            res[k+1, j] -= accumulation * (u[k+1, j] - previous[k+1, j])
        end
        res[1, j] -= accumulation / (1000 * c.cp[j]) *
            (c.enthalpy[j] - previous_enthalpy[j])
    end
    return res
end