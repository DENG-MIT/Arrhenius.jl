# xresid_native.jl — native half of the xresid cross-residual tool (xresid.py).
#
# Installs physical-unit states onto the native flames, re-extracts them (round
# trip), evaluates the production residual, and converts every row to Cantera's
# physical units with the exact production scale factors. One Julia process per
# invocation handles every job; xresid.py orchestrates the Cantera side.
#
# Usage: julia --project=<worktree> xresid_native.jl JOBS.json
#   JOBS.json (YAML/JSON) = {"spec": <spec.json>, "discretization": "finite_difference",
#     "extract":   [{"id","state","out","out_json"}...],   # install + re-extract + round trip
#     "solve":     [{"id","out","out_json"[,"initial_points"]},  # from-scratch solve -> state
#     "residual":  [{"id","state","out_npz","out_json"}...]}
# Canonical state npz (physical units; what xresid.py writes and reads back):
#   grid (N,), T (N,) [K], Y (N,n) points-major, velocity (N,) [m/s] (u, not rho*u),
#   spreadRate/Lambda (N,) for counterflow kinds, species_names_utf8.
#   ct-style dumps (values/<name> + (n,N) Y) also load.
#
# Native row -> physical scale factors (src refs in audit.md), exact from
# counterflow_residual! (src/CounterflowFlames.jl:348-455), the twin wall
# override (src/PremixedCounterflowFlames.jl:222-241) and flame_residual!
# (src/Flames.jl:571-637); ts=_flame_timescale, Vs=_counterflow_Vscale,
# Ls=_counterflow_Lscale. The interior residual rows already divide their
# bracket by rho (species), rho*cp (T) or rho*Vs (V) internally, so converting
# to Cantera's row only removes the time (and state) scaling - no property
# factor appears (verified on the h2o2 burner mini: with a rho*cp factor the
# native/CT ratio is exactly -rho*cp for T and rho for species):
#   counterflow:
#     T interior  x 1000/ts [K/s]
#     T boundary  x 1000 [K]
#     Y_k interior x 1/ts [1/s], boundary x 1 [kg/m2/s]
#     V interior  x Vs/ts, boundary x Vs [m/s]
#     Lambda interior x Ls [Pa], left x (-1) [kg/m2/s]  (CT: mdot - rho*u)
#     mdot interior x (-1/dz_j) [kg/m3/s], right x 1 [kg/m2/s]
#     dependent-species row (sum Y - 1) x 1 everywhere (CT: interior PDE,
#       boundary 1 - sum Y on its excess species)
#   free/burner/fixed:
#     T interior x 1000/ts, boundary/imposed x 1000
#     Y_k interior x 1/ts, j=1 flux row x 1, j=N zero-grad x 1
#     mdot: burner j=1 x 1 and interior x 1 (undivided, as in CT unstrained);
#           free j=1/interior x (-1/dz), anchor x 1000 (T-Tfix [K]), j=N x 1
# Physical residual rows and per-(row,eqtype) term scales go to out_npz /
# out_json; xresid.py applies the same scales to Cantera's rows for the common
# relative norm  rel(row) = max_j |r_row(j)| / S_row,  S_row = max over points
# and additive terms of |term| on the scales state (the converged solution).
using Arrhenius, YAML, NPZ, LinearAlgebra
BLAS.set_num_threads(1)

const VS = Arrhenius._counterflow_Vscale
const LS = Arrhenius._counterflow_Lscale
const TS = Arrhenius._flame_timescale
const KINDS_CF = ("diffusion", "premixed", "twin", "stagnation")
# Equation-type codes written to the npz (type_codes matrix).
const EQTYPES = ["energy", "species", "momentum", "continuity", "lambda",
    "lambda_left", "T_bc", "flux_bc", "zero_grad", "closure", "mdot_bc", "T_anchor",
    "zg_T", "zg_Y", "zg_V"]

parse_X(text) = Dict(String(n) => parse(Float64, v)
    for (n, v) in split.(split(text, ","), ":"))

is_counterflow(spec) = spec["kind"] in KINDS_CF
utf8(s) = Vector{UInt8}(join(s, "\n"))

function build(spec, gas, grid; discretization=:finite_difference)
    kind = spec["kind"]
    P = Float64(spec["P_atm"]) * one_atm
    T_in = Float64(spec["T"])
    transport = Symbol(spec["transport"])
    width = grid[end] - grid[1]
    if kind == "free"
        return FreeFlame(gas; T=T_in, P, X=parse_X(spec["X"]), width,
            transport_model=transport, discretization, grid)
    elseif kind in ("burner", "fixed")
        f = BurnerFlame(gas; mdot=spec["mdot"], T=T_in, P, X=parse_X(spec["X"]),
            width, transport_model=transport, discretization, grid)
        kind == "fixed" && let z = Float64.(spec["profile_z"])
            # val_final_default specs store ABSOLUTE z (0..width); a
            # self-authored spec may store normalized z. Normalize either way
            # (identity for already-normalized input), matching xresid.py.
            set_temperature_profile!(f, (z .- z[1]) ./ (z[end] - z[1]),
                Float64.(spec["profile_T"]); relative=true)
        end
        return f
    elseif kind == "diffusion"
        return CounterflowDiffusionFlame(gas; fuel=spec["fuel"],
            oxidizer=spec["oxidizer"], mdot_fuel=spec["mdot_fuel"],
            mdot_oxidizer=spec["mdot_oxidizer"], T_fuel=T_in,
            T_oxidizer=Float64(get(spec,"T_oxidizer",T_in)), P, width, grid)
    elseif kind == "premixed"
        return CounterflowPremixedFlame(gas; reactants=spec["X"],
            mdot_reactants=spec["mdot"], mdot_products=spec["mdot_products"],
            T_reactants=T_in, P, width, grid)
    elseif kind == "twin"
        return CounterflowTwinPremixedFlame(gas; reactants=spec["X"],
            mdot=spec["mdot"], T=T_in, P, width, grid)
    elseif kind == "stagnation"
        return ImpingingJet(gas; reactants=spec["X"], mdot=spec["mdot"],
            T_inlet=T_in, T_surface=Float64(spec["T_surface"]), P, width, grid)
    end
    error("unsupported kind $kind")
end

"Load a canonical or ct-style state npz into physical arrays."
function canonical(state, gas)
    z = vec(state["grid"])
    T = haskey(state, "T") ? vec(state["T"]) : vec(state["values/T"])
    u = haskey(state, "velocity") ? vec(state["velocity"]) :
        vec(state["values/velocity"])
    n, N = gas.n_species, length(z)
    Y = if haskey(state, "Y")
        raw = state["Y"]
        size(raw) == (n, N) && size(raw) != (N, n) ? permutedims(raw) : raw
    else
        reduce(hcat, [vec(state["values/" * s]) for s in gas.species_names])
    end
    size(Y) == (N, n) || error("Y has shape $(size(Y)), expected ($N, $n)")
    out = Dict{String,Any}("grid" => z, "T" => T, "velocity" => u, "Y" => Y)
    for key in ("spreadRate", "Lambda")
        haskey(state, key) ? out[key] = vec(state[key]) :
            haskey(state, "values/" * key) && (out[key] = vec(state["values/" * key]))
    end
    return out
end

"Install physical arrays onto the native state (the inverse of extract)."
function install!(f, can, gas, spec)
    z, T, u, Y = can["grid"], can["T"], can["velocity"], can["Y"]
    n, N = gas.n_species, length(z)
    meanMW = 1 ./ (Y * (1 ./ gas.MW))
    rho = f.pressure .* meanMW ./ (Arrhenius.R .* T)
    state = zeros(eltype(T), is_counterflow(spec) ? n + 4 : n + 2, N)
    state[1, :] .= T ./ 1000
    state[2:n+1, :] .= permutedims(Y)
    if is_counterflow(spec)
        haskey(can, "spreadRate") || error("counterflow state needs spreadRate")
        haskey(can, "Lambda") || error("counterflow state needs Lambda")
        state[n+2, :] .= can["spreadRate"] ./ VS
        state[n+3, :] .= can["Lambda"] ./ LS
    end
    state[end, :] .= rho .* u
    f.grid = z
    f.state = state
    if spec["kind"] == "free"
        # Re-select the fixed-temperature anchor on the installed profile, then
        # pin its temperature from the state (mirrors the ct-state loader).
        f.anchor = clamp(argmin(abs.(state[1, :] .- (f.fixed_temperature / 1000))), 2, N - 1)
        f.fixed_temperature = 1000 * state[1, f.anchor]
    end
    return f
end

"Physical-unit extraction of the installed native state (inverse of install!)."
function extract(f, gas, spec)
    n = gas.n_species
    out = Dict{String,Any}("grid" => copy(f.grid), "T" => temperature(f),
        "Y" => copy(permutedims(f.state[2:n+1, :])), "velocity" => velocity(f),
        "species_names_utf8" => utf8(gas.species_names))
    is_counterflow(spec) || return out
    out["spreadRate"] = spread_rate(f)
    out["Lambda"] = pressure_curvature(f)
    return out
end

"Componentwise round-trip error: max|a-b| / max(max|a|, max|b|)."
roundtrip(a, b) = maximum(abs.(a .- b)) / max(maximum(abs.(a)), maximum(abs.(b)), 1e-300)

function roundtrips(f, can, spec)
    back = extract(f, f.gas, spec)
    rt = Dict{String,Float64}()
    for key in ("T", "velocity", "spreadRate", "Lambda")
        haskey(back, key) && (rt[key] = roundtrip(Float64.(can[key]), Float64.(back[key])))
    end
    rt["Y"] = roundtrip(Float64.(can["Y"]), Float64.(back["Y"]))
    return rt
end

"Per-row, per-point equation-type classification (codes index EQTYPES)."
function row_eqtypes(f, spec)
    n, N = f.gas.n_species, length(f.grid)
    B = size(f.state, 1)
    code(t) = findfirst(==(t), EQTYPES)
    types = fill(code("closure"), B, N)   # closure rows overwritten below
    for j in 1:N
        boundary = j == 1 || j == N
        if is_counterflow(spec)
            types[1, j] = boundary ? code("T_bc") : code("energy")
            excess=argmax(@view f.state[2:n+1,j])
            for k in 1:n
                boundary && k==excess && continue
                types[k+1,j]=boundary ? code("flux_bc") : code("species")
            end
            types[n+2, j] = boundary ? code("zg_V") : code("momentum")
            types[n+3, j] = j == 1 ? code("lambda_left") : code("lambda")
            types[end, j] = j == N ? code("mdot_bc") : code("continuity")
        else
            types[1, j] = boundary || spec["kind"] == "fixed" ? code("T_bc") : code("energy")
            for k in 1:n
                boundary && k == f.dependent_species && continue
                types[k+1, j] = j == 1 ? code("flux_bc") : j == N ? code("zero_grad") : code("species")
            end
            if spec["kind"] in ("burner", "fixed")
                types[end, j] = j == 1 ? code("mdot_bc") : code("continuity")
            elseif j == f.anchor
                types[end, j] = code("T_anchor")
            else
                # j == 1 (anchor >= 2) and j == N included: the free-flame
                # mass-flux row is continuity below the anchor, zero-gradient
                # at the outlet (both undivided forms see physical_residual!).
                types[end, j] = j == N ? code("zero_grad") : code("continuity")
            end
        end
    end
    if spec["kind"] == "twin"
        # Wall overrides (PremixedCounterflowFlames.jl:222-241): zero-gradient
        # Y/V/T rows replace the flux forms at the symmetry plane.
        types[1, N] = isempty(f.fixed_temperature) ? code("zg_T") : code("T_bc")
        excess=argmax(@view f.state[2:n+1,N])
        for k in 1:n
            k == excess || (types[k+1,N] = code("zg_Y"))
        end
        types[n+2, N] = code("zg_V")
    end
    return types
end

"Convert the native residual to Cantera's physical units, row by row."
function physical_residual!(R, r, f, spec, w)
    n, N = f.gas.n_species, length(f.grid)
    z = f.grid
    for j in 1:N
        boundary = j == 1 || j == N
        if is_counterflow(spec)
            # Interior rows are the Cantera row times the time scale (and the
            # state scaling for T/V): counterflow_residual! divides its bracket
            # by rho (species), rho*cp (T), rho*Vs (V) internally, so no extra
            # property factor may appear here. Verified: native/CT = -rho*cp
            # (T) and rho (species) with the old factors on the burner mini.
            R[1, j] = boundary ? 1000r[1, j] : 1000 * r[1, j] / TS
            for k in 1:n
                R[k+1, j] = boundary ? r[k+1, j] : r[k+1, j] / TS
            end
            R[n+2, j] = boundary ? VS * r[n+2, j] : VS * r[n+2, j] / TS
            R[n+3, j] = j == 1 ? -r[n+3, j] : LS * r[n+3, j]
            R[end, j] = j == N ? r[end, j] : -r[end, j] / (z[j+1] - z[j])
        else
            R[1, j] = boundary || spec["kind"] == "fixed" ? 1000r[1, j] :
                1000 * r[1, j] / TS
            for k in 1:n
                R[k+1, j] = boundary ? r[k+1, j] : r[k+1, j] / TS
            end
            if spec["kind"] in ("burner", "fixed")
                R[end, j] = r[end, j]
            elseif j == f.anchor
                R[end, j] = 1000r[end, j]
            elseif j == N
                R[end, j] = r[end, j]
            else
                R[end, j] = -r[end, j] / (j < f.anchor ? z[j+1] - z[j] : z[j] - z[j-1])
            end
        end
    end
    return R
end

# Term scales: each row's additive physical terms on the installed state,
# reduced to S(row,eqtype) = max over points and terms of |term| in the same
# physical units as the row. The counterflow terms mirror counterflow_residual!
# (cfscheme-verified decomposition); the free/burner terms follow
# flame_residual! (src/Flames.jl:601-620).
function term_scales(f, spec, w)
    n, N = f.gas.n_species, length(f.grid)
    z, MW = f.grid, f.gas.MW
    # CounterflowWorkspace nests the property arrays under .properties;
    # FlameWorkspace (free/burner/fixed) carries them directly.
    p = is_counterflow(spec) ? w.properties : w
    u = f.state
    scales = Dict{String,Float64}()
    note!(k, v) = scales[k] = max(get(scales, k, 0.0), v)
    scales["closure"] = 1.0
    if is_counterflow(spec)
        iv, il, im = size(u, 1) - 2, size(u, 1) - 1, size(u, 1)
        for j in 1:N
            if j < N
                # terms of the row (not their difference): d(rho*u)/dz + rho*V
                dz = z[j+1] - z[j]
                note!("continuity", max(abs(u[im, j+1]) / dz, abs(u[im, j]) / dz,
                    abs(VS * (p.rho[j+1] * u[iv, j+1] + p.rho[j] * u[iv, j]))))
            end
            note!("lambda", j == 1 ? abs(LS * u[il, j]) :
                max(abs(LS * u[il, j]), abs(LS * u[il, j-1])))
            if j == 1 || j == N
                flux = j == 1 ? f.fuel_mass_flux : f.oxidizer_mass_flux
                yin = j == 1 ? f.fuel_Y : f.oxidizer_Y
                for k in 1:n
                    note!("flux_bc/$k", max(abs(flux * yin[k]), abs(u[im, j] * u[k+1, j]),
                        abs(p.flux[k, j == 1 ? 1 : N-1])))
                end
                note!("T_bc", max(abs(1000u[1, j]), abs(j == 1 ? f.fuel_temperature :
                    f.oxidizer_temperature)))
                note!("zg_V", abs(VS * u[iv, j]))
                j == 1 && note!("lambda_left", max(abs(u[im, 1]), abs(f.fuel_mass_flux)))
                j == N && note!("mdot_bc", max(abs(u[im, N]), abs(f.oxidizer_mass_flux)))
            else
                left = z[j]-z[j-1]; right = z[j+1]-z[j]; cell = .5(left + right)
                jm, jp = u[im, j] > 0 ? (j-1, j) : (j, j+1)
                dz = z[jp]-z[jm]
                V = VS * u[iv, j]
                shear = VS * (w.viscosity[j] * (u[iv, j+1]-u[iv, j]) / right -
                    w.viscosity[j-1] * (u[iv, j]-u[iv, j-1]) / left) / cell
                note!("momentum", max(abs(shear), abs(LS * u[il, j]),
                    abs(u[im, j] * VS * (u[iv, jp]-u[iv, jm]) / dz),
                    abs(p.rho[j] * V^2)) / p.rho[j])
                dT = 1000 * (u[1, jp]-u[1, jm]) / dz
                cpmean = sum(u[k+1, j] * p.cp[k, j] / MW[k] for k in 1:n)
                chemical = sum(p.h[k, j] * p.source[k, j] for k in 1:n)
                enthalpyflux = sum(.5 * (p.flux[k, j-1] + p.flux[k, j]) *
                    (p.h[k, jp]-p.h[k, jm]) / (dz * MW[k]) for k in 1:n)
                conduction = 1000 * (p.conductivity[j] * (u[1, j+1]-u[1, j]) / right -
                    p.conductivity[j-1] * (u[1, j]-u[1, j-1]) / left) / cell
                note!("energy", max(abs(conduction), abs(chemical), abs(enthalpyflux),
                    abs(u[im, j] * cpmean * dT), abs(w.radiative_heat_loss[j])) /
                    (p.rho[j] * cpmean))
                for k in 1:n
                    note!("species/$k", max(abs(MW[k] * p.source[k, j]),
                        abs((p.flux[k, j]-p.flux[k, j-1]) / cell),
                        abs(u[im, j] * (u[k+1, jp]-u[k+1, jm]) / dz)) / p.rho[j])
                end
            end
        end
        if spec["kind"] == "twin"
            note!("zg_T", max(abs(1000u[1, N]), abs(1000u[1, N-1])))
            for k in 1:n
                note!("zg_Y/$k", max(abs(u[k+1, N]), abs(u[k+1, N-1])))
            end
        end
    else
        for j in 1:N
            if j == 1
                for k in 1:n
                    note!("flux_bc/$k", max(abs(u[end, 1] * f.inlet_Y[k]),
                        abs(u[end, 1] * u[k+1, 1]), abs(p.flux[k, 1])))
                end
                note!("T_bc", max(abs(1000u[1, 1]), abs(f.inlet_temperature)))
                if spec["kind"] in ("burner", "fixed")
                    note!("mdot_bc", abs(spec["mdot"]))
                elseif N > 1
                    note!("continuity", max(abs(u[end, 2]), abs(u[end, 1])) / (z[2] - z[1]))
                end
            elseif j == N
                for k in 1:n
                    note!("zero_grad", max(abs(u[k+1, N]), abs(u[k+1, N-1])))
                end
                note!("T_bc", max(abs(1000u[1, N]), abs(1000u[1, N-1])))
                spec["kind"] == "free" &&
                    note!("zero_grad", max(abs(u[end, N]), abs(u[end, N-1])))
            else
                left = z[j]-z[j-1]; right = z[j+1]-z[j]; cell = .5(left + right)
                mdot = u[end, j]
                dT = 1000 * (u[1, j]-u[1, j-1]) / left
                cpmean = sum(u[k+1, j] * p.cp[k, j] / MW[k] for k in 1:n)
                chemical = sum(p.h[k, j] * p.source[k, j] for k in 1:n)
                enthalpyflux = sum(.5 * (p.flux[k, j-1] + p.flux[k, j]) *
                    (p.h[k, j]-p.h[k, j-1]) / (left * MW[k]) for k in 1:n)
                conduction = 1000 * (p.conductivity[j] * (u[1, j+1]-u[1, j]) / right -
                    p.conductivity[j-1] * (u[1, j]-u[1, j-1]) / left) / cell
                note!("energy", max(abs(conduction), abs(chemical), abs(enthalpyflux),
                    abs(mdot * cpmean * dT)) / (p.rho[j] * cpmean))
                for k in 1:n
                    note!("species/$k", max(abs(MW[k] * p.source[k, j]),
                        abs((p.flux[k, j]-p.flux[k, j-1]) / cell),
                        abs(mdot * (u[k+1, j]-u[k+1, j-1]) / left)) / p.rho[j])
                end
                if spec["kind"] in ("burner", "fixed")
                    note!("continuity", abs(u[end, j] - u[end, j-1]))
                elseif j != f.anchor
                    note!("continuity", abs((j < f.anchor ? u[end, j+1] - u[end, j] :
                        u[end, j] - u[end, j-1]) /
                        (j < f.anchor ? z[j+1] - z[j] : z[j] - z[j-1])))
                end
            end
        end
        if spec["kind"] == "free"
            scales["T_anchor"] = max(abs(1000u[1, f.anchor]), abs(f.fixed_temperature))
        end
        if spec["kind"] == "fixed"
            note!("T_bc", maximum(abs.(f.imposed_temperature)))
        end
    end
    return scales
end

jesc(s) = replace(string(s), "\\" => "\\\\", "\"" => "\\\"")

function write_json(path, id, rt, f, summary, scales)
    anchor = f isa FreeFlame ? f.anchor : 0
    io = open(path, "w")
    print(io, "{\"id\": \"", jesc(id), "\", \"anchor\": ", anchor, ", \"roundtrip\": {")
    for (i, kv) in enumerate(sort(collect(rt)))
        print(io, i == 1 ? "" : ", ", "\"", jesc(kv[1]), "\": ", kv[2])
    end
    print(io, "}, \"rows\": [")
    for (i, s) in enumerate(summary)
        print(io, i == 1 ? "" : ", ", "{\"row\": \"", jesc(s.row),
            "\", \"eqtype\": \"", jesc(s.eqtype), "\", \"max\": ", s.max,
            ", \"j\": ", s.j, ", \"z\": ", s.z, "}")
    end
    print(io, "], \"scales\": {")
    for (i, kv) in enumerate(sort(collect(scales)))
        print(io, i == 1 ? "" : ", ", "\"", jesc(kv[1]), "\": ", kv[2])
    end
    println(io, "}, \"dependent_species\": ", f.dependent_species, "}")
    close(io)
end

"Residual job: physical residual + types + scales to npz, summary to json."
function residual_job(f, spec, job, rt)
    n, N = f.gas.n_species, length(f.grid)
    if is_counterflow(spec)
        w = Arrhenius.CounterflowWorkspace(f)
        r = Arrhenius.counterflow_residual!(similar(f.state), f, f.state, w)
    else
        w = Arrhenius.FlameWorkspace(f)
        r = Arrhenius.flame_residual!(similar(f.state), f, f.state, w)
    end
    R = physical_residual!(similar(r), r, f, spec, w)
    types = row_eqtypes(f, spec)
    scales = term_scales(f, spec, w)
    rows = is_counterflow(spec) ? vcat(["T"], f.gas.species_names, ["V", "Lambda", "mdot"]) :
        vcat(["T"], f.gas.species_names, ["mdot"])
    npzwrite(job["out_npz"], Dict{String,Any}("grid" => f.grid, "physical" => R, "raw" => r,
        "type_codes" => Int16.(types), "type_names" => utf8(EQTYPES),
        "row_names" => utf8(rows),
        "species_names_utf8" => utf8(f.gas.species_names),
        "dependent_species" => [f.dependent_species],
        [(k => [v]) for (k, v) in scales]...))
    # JSON summary: per (row, eqtype) max |r| physical + location, and scales.
    summary = []
    for (i, row) in enumerate(rows)
        bytype = Dict{Int,Vector{Int}}()
        for j in 1:N
            push!(get!(bytype, types[i, j], Int[]), j)
        end
        for (t, js) in sort(collect(bytype))
            vals = abs.(R[i, js])
            m = findmax(vals)
            push!(summary, (row=row, eqtype=EQTYPES[t], max=vals[m[2]], j=js[m[2]],
                z=f.grid[js[m[2]]]))
        end
    end
    write_json(job["out_json"], job["id"], rt, f, summary, scales)
    return summary, scales
end

function main(jobs)
    spec = YAML.load_file(jobs["spec"])
    haskey(spec, "kind") || (spec = spec["spec"])  # reproducer jsons nest it
    disc = Symbol(get(jobs, "discretization", "finite_difference"))
    gas = CreateSolution(spec["mechanism"] * ".yaml")
    for job in get(jobs, "extract", Any[])
        can = canonical(npzread(job["state"]), gas)
        f = build(spec, gas, can["grid"]; discretization=disc)
        install!(f, can, gas, spec)
        rt = roundtrips(f, can, spec)
        npzwrite(job["out"], extract(f, gas, spec))
        write_json(job["out_json"], job["id"], rt, f, [], Dict{String,Float64}())
        println("xresid_native extract $(job["id"]): N=$(length(f.grid)) ",
            "anchor=$(f isa FreeFlame ? f.anchor : 0) roundtrip=",
            join(["$k=$(round(v, sigdigits=3))" for (k, v) in sort(collect(rt))], " "))
    end
    for job in get(jobs, "solve", Any[])
        # From-scratch native solve (the extract/residual jobs reinstall
        # states onto exact grids instead). Refine criteria come from the
        # spec, matching the validation harness.
        n0 = get(job, "initial_points", 16)
        grid = collect(range(0.0, Float64(spec["width"]); length=n0))
        f = build(spec, gas, grid; discretization=disc)
        # Clamp prune exactly like Cantera's Refiner::setCriteria requires
        # (ct_dump.py and the CT side do the same).
        slope, curve = Float64(get(spec, "slope", 0.06)), Float64(get(spec, "curve", 0.12))
        prune = min(Float64(get(spec, "prune", 0.0)), slope, curve)
        solve!(f; refine_grid=true,
            ratio=Float64(get(spec, "ratio", 3.0)), slope, curve,
            prune, loglevel=0)
        npzwrite(job["out"], extract(f, gas, spec))
        write_json(job["out_json"], job["id"], Dict{String,Float64}(), f,
            Any[], Dict{String,Float64}())
        println("xresid_native solve $(job["id"]): N=$(length(f.grid)) ",
            "converged=", f.converged)
    end
    for job in get(jobs, "residual", Any[])
        can = canonical(npzread(job["state"]), gas)
        f = build(spec, gas, can["grid"]; discretization=disc)
        install!(f, can, gas, spec)
        rt = roundtrips(f, can, spec)
        summary, scales = residual_job(f, spec, job, rt)
        worst = sort(summary; by=s -> -s.max)[1:min(end, 5)]
        println("xresid_native residual $(job["id"]): N=$(length(f.grid)) ",
            "anchor=$(f isa FreeFlame ? f.anchor : 0) rt_Y=$(round(rt["Y"], sigdigits=3)) worst=",
            join(["$(s.row)/$(s.eqtype)=$(round(s.max, sigdigits=3))@j$(s.j)" for s in worst], " "))
    end
end

jobs = YAML.load_file(ARGS[1])  # jobs files are YAML-compatible JSON
if abspath(PROGRAM_FILE) == @__FILE__
    main(jobs)
end
