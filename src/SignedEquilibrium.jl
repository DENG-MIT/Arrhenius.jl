# Stable positive/negative conservation balances for charged ideal-gas equilibrium.

function _signed_equilibrium_system(gas, X)
    all(isfinite, X) || throw(ArgumentError("mole fractions must be finite"))
    E = gas.ele_matrix
    all(isfinite, E) || throw(ArgumentError("element matrix must be finite"))
    b_all = E * X
    all(isfinite, b_all) || throw(ArgumentError("elemental targets must be finite"))
    nrows, nsp = size(E)
    kept = trues(nsp)
    dropped = falses(nrows)
    # Iterate species removal to a fixed point: a zero-target row whose
    # remaining atom counts are one-sided forbids every species carrying a
    # nonzero count on that row. Removing those species can make further rows
    # one-sided (e.g. a lone charge carrier once its partners vanish), so the
    # sweep repeats until stable. Two-sided signed rows with b == 0 are never
    # used as a removal criterion.
    changed = true
    while changed
        changed = false
        @inbounds for e in 1:nrows
            dropped[e] && continue
            target = b_all[e]
            haspos = false
            hasneg = false
            for k in 1:nsp
                kept[k] || continue
                c = E[e, k]
                c > 0 && (haspos = true)
                c < 0 && (hasneg = true)
                haspos && hasneg && break
            end
            if target == 0
                if !haspos && !hasneg
                    dropped[e] = true
                    changed = true
                elseif !haspos || !hasneg
                    for k in 1:nsp
                        if kept[k] && E[e, k] != 0
                            kept[k] = false
                            changed = true
                        end
                    end
                end
            else
                feasible = target > 0 ? haspos : hasneg
                feasible || throw(ArgumentError(
                    "composition target is infeasible for a conserved row"))
            end
        end
    end
    species = findall(kept)
    isempty(species) && throw(ArgumentError(
        "no species are compatible with the zero-abundance balances"))
    rows = findall(.!dropped)
    any(e -> b_all[e] != 0, rows) || throw(ArgumentError(
        "composition must contain a conserved element with nonzero abundance"))
    A = Float64.(E[rows, species])
    @inbounds for k in axes(A, 2)
        any(!=(0), view(A, :, k)) || throw(ArgumentError(
            "species without conserved elements are unsupported"))
    end
    # Select independent original rows, preserving row/target correspondence.
    # An SVD basis could mix rows and destroy the signed two-sided structure.
    selector = qr(transpose(A), ColumnNorm())
    diagonal = abs.(diag(selector.R))
    rankA = count(>(maximum(diagonal) * max(size(A)...) * eps(Float64)), diagonal)
    independent = selector.p[1:rankA]
    A = A[independent, :]
    b = Float64.(b_all[rows[independent]])
    ne, ns = size(A)
    logb = Vector{Float64}(undef, ne)
    @inbounds for i in 1:ne
        logb[i] = b[i] == 0 ? -Inf : log(abs(b[i]))
    end
    logabs_atoms = Matrix{Float64}(undef, ne, ns)
    positive_indices = Vector{Vector{Int}}(undef, ne)
    negative_indices = Vector{Vector{Int}}(undef, ne)
    @inbounds for i in 1:ne
        positive_indices[i] = findall(>(0), view(A, i, :))
        negative_indices[i] = findall(<(0), view(A, i, :))
        for k in 1:ns
            a = A[i, k]
            logabs_atoms[i, k] = a == 0 ? -Inf : log(abs(a))
        end
    end
    return (gas=gas, X=X, species=species, A=A, At=Matrix(transpose(A)),
            b=b, state=zeros(ne+1),
            g=zeros(ns), initialized=Ref(false), logb=logb,
            mole=zeros(ns), potentials=zeros(ns), abar=zeros(ne),
            residual=zeros(ne+1),
            jacobian=zeros(ne+1, ne+1), factor=zeros(ne+1, ne+1),
            trial=zeros(ne+1), step=zeros(ne+1),
            signed_balance=true,
            side_plus=zeros(ne), side_minus=zeros(ne),
            positive_indices=positive_indices,
            negative_indices=negative_indices,
            logabs_atoms=logabs_atoms)
end

# src/SignedEquilibrium.jl — signed (charged) element-balance equilibrium
# evaluator. Included after Equilibrium.jl. Expects a system NamedTuple with
# the shared equilibrium fields plus signed_balance=true, side_plus,
# side_minus, positive_indices, negative_indices and logabs_atoms.

function _signed_logaddexp(a::Float64, b::Float64)
    a == -Inf && return b
    b == -Inf && return a
    hi = max(a, b)
    lo = min(a, b)
    return hi + log1p(exp(lo - hi))
end

function _signed_side_logsumexp(v, logabs_atoms, i, indices)
    isempty(indices) && return -Inf
    vmax = -Inf
    @inbounds for k in indices
        vmax = max(vmax, v[k] + logabs_atoms[i, k])
    end
    vmax == -Inf && return -Inf
    total = 0.0
    @inbounds for k in indices
        total += exp(v[k] + logabs_atoms[i, k] - vmax)
    end
    return vmax + log(total)
end

function _signed_side_moments!(side, A, v, logabs_atoms, i, indices, q)
    fill!(side, 0.0)
    q == -Inf && return side
    ne = size(A, 1)
    @inbounds for k in indices
        w = exp(v[k] + logabs_atoms[i, k] - q)
        for j in 1:ne
            side[j] += w * A[j, k]
        end
    end
    return side
end

function _signed_equilibrium_evaluate(system, state, logpressure, constant_volume; jacobian=false)
    A, At, g, b = system.A, system.At, system.g, system.b
    ne = size(A, 1)
    v, mole, abar, f = system.potentials, system.mole, system.abar, system.residual
    logabs = system.logabs_atoms
    pos, neg = system.positive_indices, system.negative_indices
    mul!(v, At, view(state, 1:ne))
    logN = state[end]
    offset = logpressure + (constant_volume ? logN : 0.0)
    @. v = v - g - offset
    vmax = maximum(v)
    total = 0.0
    @inbounds for k in eachindex(v)
        total += exp(v[k] - vmax)
    end
    z = vmax + log(total)
    @. mole = exp(v - z)
    mul!(abar, A, mole)
    J = system.jacobian
    side_plus, side_minus = system.side_plus, system.side_minus
    @inbounds for i in 1:ne
        qplus = _signed_side_logsumexp(v, logabs, i, pos[i])
        qminus = _signed_side_logsumexp(v, logabs, i, neg[i])
        argplus = logN + qplus - z
        argminus = logN + qminus - z
        lres = _signed_logaddexp(argplus, log(max(-b[i], 0.0)))
        rres = _signed_logaddexp(argminus, log(max(b[i], 0.0)))
        f[i] = lres - rres
        if jacobian
            _signed_side_moments!(side_plus, A, v, logabs, i, pos[i], qplus)
            _signed_side_moments!(side_minus, A, v, logabs, i, neg[i], qminus)
            wl = argplus == -Inf ? 0.0 : exp(argplus - lres)
            wr = argminus == -Inf ? 0.0 : exp(argminus - rres)
            for j in 1:ne
                J[i, j] = wl * (side_plus[j] - abar[j]) - wr * (side_minus[j] - abar[j])
            end
            J[i, end] = wl - wr
        end
    end
    f[end] = z
    jacobian || return f, mole
    @inbounds for j in 1:ne
        J[end, j] = abar[j]
    end
    J[end, end] = constant_volume ? -1.0 : 0.0
    return f, mole, J
end
