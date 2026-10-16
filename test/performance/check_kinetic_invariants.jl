# Recheck state closure and ideal-gas flame invariants on saved native profiles.
# Usage: julia --project test/performance/check_kinetic_invariants.jl MANIFEST OUT_DIR
using Arrhenius, YAML, NPZ, LinearAlgebra
BLAS.set_num_threads(1)

function check_case(id, spec_path, out_dir)
    spec = YAML.load_file(spec_path)
    gas = Arrhenius.CreateSolution(spec["mechanism"] * ".yaml")
    data = NPZ.npzread(joinpath(out_dir, id * ".npz"))
    T = Float64.(data["T"])
    Y = Float64.(data["Y"])
    velocity = Float64.(data["velocity"])
    P = Float64(spec["P_atm"] * Arrhenius.one_atm)
    n, N = size(Y)
    n == gas.n_species || error("$id: profile species dimension mismatch")
    length(T) == N == length(velocity) || error("$id: profile node dimension mismatch")
    sum_y_error = maximum(abs, vec(sum(Y; dims=1) .- 1.0); init=0.0)
    w = Arrhenius.KineticsWorkspace(gas.reaction)
    source = zeros(n)
    C = zeros(n)
    h = zeros(n); cp = zeros(n); s = zeros(n)
    max_element_relative = 0.0
    max_element_absolute = 0.0
    mass_flux_relative = NaN
    if spec["kind"] in ("free", "burner", "fixed")
        mass_flux = Vector{Float64}(undef, N)
        for j in 1:N
            inverse_mw = sum(@view(Y[:,j]) ./ gas.MW)
            rho = P / (Arrhenius.R * T[j] * inverse_mw)
            mass_flux[j] = rho * velocity[j]
        end
        mean_flux = sum(mass_flux) / N
        mass_flux_relative = maximum(abs, mass_flux .- mean_flux) / max(abs(mean_flux), 1e-300)
    end
    for j in 1:N
        X = @view(Y[:,j]) ./ gas.MW
        X ./= sum(X)
        C .= X .* (P / (Arrhenius.R * T[j]))
        Arrhenius.cal_h_cp_s0!(h,cp,s,gas,T[j],P,X)
        h .*= Arrhenius.R * T[j]
        s .*= Arrhenius.R
        Arrhenius.wdot!(source,gas.reaction,T[j],C,s,h,w;pressure=P)
        for e in axes(gas.ele_matrix,1)
            balance = 0.0
            scale = 0.0
            @inbounds for k in 1:n
                term = gas.ele_matrix[e,k] * source[k]
                balance += term
                scale += abs(term)
            end
            max_element_absolute = max(max_element_absolute,abs(balance))
            max_element_relative = max(max_element_relative,abs(balance)/max(scale,1e-300))
        end
    end
    println("INVARIANT case=$id nodes=$N sumY_max=$sum_y_error element_source_abs_max=$max_element_absolute element_source_relative_max=$max_element_relative mass_flux_relative_max=$mass_flux_relative")
    return sum_y_error,max_element_absolute,max_element_relative,mass_flux_relative
end

manifest, out_dir = ARGS
sum_y_max = 0.0; element_abs_max = 0.0; element_relative_max = 0.0
mass_flux_relative_max = 0.0; mass_flux_cases = 0; cases = 0
for line in eachline(manifest)
    isempty(strip(line)) && continue
    id,spec_path = split(line, ',', limit=2)
    sy,ea,er,mf = check_case(id,spec_path,out_dir)
    sum_y_max = max(sum_y_max,sy)
    element_abs_max = max(element_abs_max,ea)
    element_relative_max = max(element_relative_max,er)
    if isfinite(mf)
        mass_flux_relative_max = max(mass_flux_relative_max,mf)
        mass_flux_cases += 1
    end
    cases += 1
end
println("INVARIANT_SUM cases=$cases sumY_max=$sum_y_max element_source_abs_max=$element_abs_max element_source_relative_max=$element_relative_max mass_flux_relative_max=$mass_flux_relative_max mass_flux_cases=$mass_flux_cases")
