"""
    flame_speed_sensitivities(f::FreeFlame)

Normalized laminar flame speed sensitivities ``s_i = (k_i/S_u)\\,dS_u/dk_i`` to the
rate constant of every reaction at a converged conservative free flame. Construct
with `discretization=:conservative` and solve before calling. One
adjoint solve with the steady Newton Jacobian gives ``λ`` with ``Jᵀλ = ∂S_u/∂u``;
scaling reaction ``i`` scales its net rate of progress, which enters only the
species balances of each node, so ``s_i = -λᵀ ∂F/∂\\ln k_i / S_u``. The adjoint
uses a finite-difference Jacobian with transport coefficients updated in every
column, as Cantera's adjoint sensitivities do.
"""
function flame_speed_sensitivities(f::FreeFlame)
    f.converged || throw(ArgumentError("solve the flame before computing sensitivities"))
    _conservative_flame(f) || throw(ArgumentError("sensitivities require the conservative discretization"))
    gas = f.gas; reaction = gas.reaction; n = gas.n_species
    u = f.state; B,N = size(u); bw = 2B-1
    w = FlameWorkspace(f)
    r = flame_residual!(similar(u),f,u,w)
    band = _flame_jacobian(f,u,w,r;analytic=false,update_transport=true)
    _,pivots = LinearAlgebra.LAPACK.gbtrf!(bw,bw,length(u),band)
    # S_u depends only on the inlet node through its density and mass flux.
    su = flame_speed(f); gradient = zeros(B,N)
    for k in 1:B
        step = 1e-7*max(abs(u[k,1]),k == 1 ? .1 : 1e-5)
        saved = u[k,1]
        u[k,1] = saved+step; gradient[k,1] = (flame_speed(f)-su)/step
        u[k,1] = saved
    end
    adjoint = vec(gradient)
    LinearAlgebra.LAPACK.gbtrs!('T',bw,bw,length(u),band,pivots,adjoint)
    adjoint = reshape(adjoint,B,N)
    sensitivities = zeros(reaction.n_reactions)
    progress = zeros(reaction.n_reactions); entropy = zeros(n); concentration = zeros(n)
    kinetics = KineticsWorkspace(reaction)
    rows,coefficients = rowvals(reaction.vk),nonzeros(reaction.vk)
    @inbounds for j in 2:N
        T = 1000u[1,j]; factor = _flame_timescale/w.rho[j]
        for k in 1:n
            concentration[k] = f.pressure/(R*T)*w.X[k,j]
        end
        cal_s0_R!(entropy,gas,T,f.pressure,@view(w.X[:,j])); entropy .*= R
        copyto!(progress,wdot!(similar(concentration),reaction,T,concentration,entropy,@view(w.h[:,j]),kinetics;
            get_qdot=true,pressure=f.pressure))
        for i in eachindex(progress)
            q = progress[i]; q == 0 && continue
            projection = 0.0
            for slot in nzrange(reaction.vk,i)
                k = rows[slot]; k == f.dependent_species && continue
                projection += adjoint[k+1,j]*factor*gas.MW[k]*coefficients[slot]
            end
            sensitivities[i] -= projection*q
        end
    end
    return sensitivities ./ su
end
export flame_speed_sensitivities
