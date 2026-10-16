# native_dump_cf.jl — native-side counterflow dump tool mirroring ct_dump.py
# (test/cantera_parity) for the four axisymmetric counterflow kinds
# (diffusion / premixed / twin / stagnation). Diagnosis only: no src changes;
# the per-term arrays below are computed in this file by mirroring
# `counterflow_residual!` (src/CounterflowFlames.jl) line by line from the same
# property arrays, and are self-checked against the production residual.
#
# Usage (from the validation harness directory (NH3_HARNESS_DIR) so mechanism sidecars resolve):
#   julia --project=<worktree> native_dump_cf.jl spec.json outdir state CT_NPZ
#       CT_NPZ: ct_state.npz (initial guess) or ct_solved.npz (converged state)
#       written by ct_dump.py. Loads grid, T, u, V, Lambda, Y into the native
#       counterflow state on the same grid (mdot = rho*u, V and Lambda in the
#       native 1/100 / 1/10000 scalings) and dumps:
#           outdir/native_cf.npz:
#             grid, state (n+4, N native rows), residual (raw, native scaling),
#             dependent_species, mirror_error (max |mirror - production|),
#             terms/{mdot,lambda} (1, N),
#             terms/V (4, N):  shear, -Lambda, -mdot*dVdz, -rho*V^2
#             terms/T (4, N):  conduction, -chemical, -enthalpy_flux, -mdot*cp*dT
#             terms/Y (n, 3, N): MW*source, -(flux[j]-flux[j-1])/cell, -mdot*dY/dz
#           outdir/native_cf_meta.json: row names + point tags
#       Row scalings vs the physical Cantera residual follow the production
#       code: continuity ~ -dz_j, Lambda interior ~ 1/Lscale, momentum ~
#       ts/Vscale/rho factors, energy ~ -ts/1000, species ~ ts; the compare
#       script estimates them empirically per (component, region).
using Arrhenius, YAML, NPZ, LinearAlgebra
BLAS.set_num_threads(1)

# One Julia session loops over every case directory given as an argument and
# over each Cantera dump (ct_state.npz / ct_solved.npz) present in it.
const TAGS = ("ct_state", "ct_solved")

function build_native(spec, gas, grid)
    kind = spec["kind"]
    P = spec["P_atm"] * one_atm
    T_in = Float64(spec["T"])
    if kind == "diffusion"
        return CounterflowDiffusionFlame(gas; fuel=spec["fuel"],
            oxidizer=spec["oxidizer"], mdot_fuel=spec["mdot_fuel"],
            mdot_oxidizer=spec["mdot_oxidizer"], T_fuel=T_in,
            T_oxidizer=Float64(get(spec,"T_oxidizer",T_in)), P, width=grid[end]-grid[1], grid)
    elseif kind == "premixed"
        return CounterflowPremixedFlame(gas; reactants=spec["X"],
            mdot_reactants=spec["mdot"], mdot_products=spec["mdot_products"],
            T_reactants=T_in, P, width=grid[end]-grid[1], grid)
    elseif kind == "twin"
        return CounterflowTwinPremixedFlame(gas; reactants=spec["X"],
            mdot=spec["mdot"], T=T_in, P, width=grid[end]-grid[1], grid)
    elseif kind == "stagnation"
        return ImpingingJet(gas; reactants=spec["X"], mdot=spec["mdot"],
            T_inlet=T_in, T_surface=Float64(spec["T_surface"]), P,
            width=grid[end]-grid[1], grid)
    else
        error("unsupported kind $kind")
    end
end

"Load the Cantera state (grid, T, u, V, Lambda, Y) into native rows."
function load_ct_state!(f, ct, names, gas, P)
    z = vec(ct["grid"]); T = vec(ct["values/T"]); u = vec(ct["values/velocity"])
    V = vec(ct["values/spreadRate"]); Lam = vec(ct["values/Lambda"])
    names == gas.species_names || error("species order mismatch: $names vs $(gas.species_names)")
    n = gas.n_species
    N = length(z)
    Y = reduce(hcat, [vec(ct["values/$s"]) for s in names])
    meanMW = 1 ./ (Y * (1 ./ gas.MW))
    rho = P .* meanMW ./ (Arrhenius.R .* T)
    state = zeros(n + 4, N)
    state[1, :] = T ./ 1000
    state[2:n+1, :] = permutedims(Y)
    state[n+2, :] = V ./ Arrhenius._counterflow_Vscale
    state[n+3, :] = Lam ./ Arrhenius._counterflow_Lscale
    state[end, :] = rho .* u
    f.grid = z
    f.state = state
    return f
end

# Per-term mirror of counterflow_residual! (CounterflowFlames.jl:316-422) with
# the AbstractPremixedCounterflow boundary overrides (PremixedCounterflowFlames.jl)
# applied afterwards, exactly as production does. Records each additive term
# separately; the assembled mirror must reproduce the production residual to
# roundoff (reported as mirror_error).
function residual_terms(f, u, w)
    fbase = f isa CounterflowDiffusionFlame ? f : getfield(f, :flow)
    p = w.properties
    Arrhenius._flame_properties!(p, fbase, u; update_transport=true,
        nodes=eachindex(fbase.grid), freeze_source=false,
        viscosity=fbase.transport_model == :mixture_averaged ? w.viscosity : nothing)
    Arrhenius._counterflow_radiation!(w.radiative_heat_loss, w.planck_absorption,
        fbase, u, p.X)
    n, N = fbase.gas.n_species, length(fbase.grid)
    z, MW = fbase.grid, fbase.gas.MW
    ts = Arrhenius._flame_timescale
    Vs = Arrhenius._counterflow_Vscale
    Ls = Arrhenius._counterflow_Lscale
    iv, il, im = size(u, 1) - 2, size(u, 1) - 1, size(u, 1)
    r_mdot = zeros(1, N)
    r_lambda = zeros(1, N)
    t_V = zeros(4, N)   # shear, -Lambda, -mdot*dVdz, -rho*V^2
    t_T = zeros(4, N)   # conduction, -chemical, -enthalpy_flux, -mdot*cp*dT
    t_Y = zeros(n, 3, N) # MW*source, -(flux diff)/cell, -mdot*dY/dz
    mirror = zeros(size(u))
    @inbounds for j in 1:N
        if j < N
            r_mdot[1, j] = u[im, j+1] - u[im, j] + (z[j+1] - z[j]) *
                Vs * (p.rho[j+1] * u[iv, j+1] + p.rho[j] * u[iv, j])
        else
            r_mdot[1, j] = u[im, j] + fbase.oxidizer_mass_flux
        end
        r_lambda[1, j] = j == 1 ? u[im, j] - fbase.fuel_mass_flux :
            u[il, j] - u[il, j-1]
        mirror[im, j] = r_mdot[1, j]
        mirror[il, j] = r_lambda[1, j]
        if j == 1 || j == N
            mirror[iv, j] = u[iv, j]
            mirror[1, j] = u[1, j] - (j == 1 ? fbase.fuel_temperature :
                fbase.oxidizer_temperature) / 1000
            for k in 1:n
                mirror[k+1, j] = j == 1 ?
                    fbase.fuel_mass_flux * fbase.fuel_Y[k] - u[im, j] * u[k+1, j] - p.flux[k, 1] :
                    p.flux[k, N-1] + u[im, j] * u[k+1, j] +
                    fbase.oxidizer_mass_flux * fbase.oxidizer_Y[k]
            end
        else
            left = z[j] - z[j-1]; right = z[j+1] - z[j]; cell = .5 * (left + right)
            jm, jp = u[im, j] > 0 ? (j-1, j) : (j, j+1)
            dz = z[jp] - z[jm]
            V = Vs * u[iv, j]
            dV = Vs * (u[iv, jp] - u[iv, jm]) / dz
            shear = Vs * (w.viscosity[j] * (u[iv, j+1] - u[iv, j]) / right -
                w.viscosity[j-1] * (u[iv, j] - u[iv, j-1]) / left) / cell
            t_V[1, j] = shear
            t_V[2, j] = -Ls * u[il, j]
            t_V[3, j] = -u[im, j] * dV
            t_V[4, j] = -p.rho[j] * V^2
            mirror[iv, j] = ts / (p.rho[j] * Vs) * sum(t_V[:, j])
            dT = 1000 * (u[1, jp] - u[1, jm]) / dz
            cpmean = 0.; chemical = 0.; enthalpyflux = 0.
            for k in 1:n
                cpmean += u[k+1, j] * p.cp[k, j] / MW[k]
                chemical += p.h[k, j] * p.source[k, j]
                enthalpyflux += .5 * (p.flux[k, j-1] + p.flux[k, j]) *
                    (p.h[k, jp] - p.h[k, jm]) / (dz * MW[k])
                t_Y[k, 1, j] = MW[k] * p.source[k, j]
                t_Y[k, 2, j] = -(p.flux[k, j] - p.flux[k, j-1]) / cell
                t_Y[k, 3, j] = -u[im, j] * (u[k+1, jp] - u[k+1, jm]) / dz
                mirror[k+1, j] = ts / p.rho[j] * sum(t_Y[k, :, j])
            end
            conduction = 1000 * (p.conductivity[j] * (u[1, j+1] - u[1, j]) / right -
                p.conductivity[j-1] * (u[1, j] - u[1, j-1]) / left) / cell
            t_T[1, j] = conduction
            t_T[2, j] = -chemical
            t_T[3, j] = -enthalpyflux
            t_T[4, j] = -u[im, j] * cpmean * dT
            mirror[1, j] = ts / (1000 * p.rho[j] * cpmean) * sum(t_T[:, j])
        end
        mirror[fbase.dependent_species+1, j] = sum(@view(u[2:n+1, j])) - 1
    end
    # Boundary overrides of the premixed-counterflow wrappers, applied after
    # the base residual exactly as production orders them.
    if (f isa CounterflowTwinPremixedFlame || f isa ImpingingJet) && !f.soret_enabled
        for k in 1:n
            k == fbase.dependent_species && continue
            mirror[k+1, end] = u[k+1, end] - u[k+1, end-1]
        end
    end
    if f isa CounterflowTwinPremixedFlame
        mirror[end-2, end] = u[end-2, end] - u[end-2, end-1]
        if isempty(f.fixed_temperature)
            mirror[1, end] = u[1, end] - u[1, end-1]
        end
    end
    return (mdot=r_mdot, lambda=r_lambda, V=t_V, T=t_T, Y=t_Y, mirror=mirror)
end

function dump_state(outdir, spec, gas, ctfile)
    ct = npzread(joinpath(outdir, ctfile))
    meta = YAML.load_file(joinpath(outdir, replace(ctfile, ".npz" => "_meta.json")))
    grid = vec(ct["grid"])
    f = build_native(spec, gas, grid)
    load_ct_state!(f, ct, String.(meta["species_names"]), gas,
        spec["P_atm"] * one_atm)
    w = Arrhenius.CounterflowWorkspace(f)
    r = Arrhenius.counterflow_residual!(similar(f.state), f, f.state, w)
    terms = residual_terms(f, f.state, Arrhenius.CounterflowWorkspace(f))
    mirror_error = maximum(abs.(terms.mirror .- r))
    p = w.properties
    npzwrite(joinpath(outdir, "native_cf_" * ctfile), Dict(
        "grid" => f.grid, "state" => f.state, "residual" => r,
        "mirror" => terms.mirror,
        "dependent_species" => [f.dependent_species],
        "mirror_error" => [mirror_error],
        "terms/mdot" => terms.mdot, "terms/lambda" => terms.lambda,
        "terms/V" => terms.V, "terms/T" => terms.T, "terms/Y" => terms.Y,
        "props/rho" => p.rho, "props/cp" => p.cp, "props/h" => p.h,
        "props/X" => p.X, "props/source" => p.source, "props/flux" => p.flux,
        "props/conductivity" => p.conductivity, "props/viscosity" => w.viscosity))
    if mirror_error > 1e-10
        d = abs.(terms.mirror .- r)
        for j in axes(d, 2), k in axes(d, 1)
            d[k, j] > 1e-10 && println("  MIRROR-DIFF row k=$k j=$j ",
                "prod=$(r[k,j]) mirror=$(terms.mirror[k,j])")
        end
    end
    open(joinpath(outdir, "native_cf_" * replace(ctfile, ".npz" => "_meta.json")), "w") do io
        names = join(("\"" * s * "\"" for s in
            ["T"; gas.species_names; "V"; "Lambda"; "mdot"]), ", ")
        write(io, "{\"row_names\": [$names], \"kind\": \"$(spec["kind"])\"}")
    end
    println("native_dump_cf $ctfile: N=$(length(f.grid)) dependent=",
        f.dependent_species, " (", gas.species_names[f.dependent_species], ")",
        " mirror_error=$mirror_error")
    return mirror_error
end

for dir in ARGS
    spec = YAML.load_file(joinpath(dir, "spec.json"))
    gas = CreateSolution(spec["mechanism"] * ".yaml")
    for tag in TAGS
        isfile(joinpath(dir, tag * ".npz")) || continue
        dump_state(dir, spec, gas, tag * ".npz")
    end
end
