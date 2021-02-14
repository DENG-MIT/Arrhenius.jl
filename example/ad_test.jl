using ForwardDiff

function f(u)
    T = u[end]
    Y = @view(u[1:end-1])
    return T
end

u0 = rand(40)
@time f(u0)
@time f(u0)

@time grad = ForwardDiff.gradient(f, u0)
@time grad = ForwardDiff.gradient(f, u0)

function dudt(u)
    T = u[end]
    Y = u[1:end-1]
    # mean_MW = 1. / dot(Y, 1 ./ gas.MW)
    # ρ_mass = P / R / T * mean_MW
    # X = Y2X(gas, Y, mean_MW)
    # C = Y2C(gas, Y, ρ_mass)
    # cp_mole, cp_mass = get_cp(gas, T, X, mean_MW)
    # h_mole = get_H(gas, T, Y, X)
    # S0 = get_S(gas, T, P, X)
    # wdot = wdot_func(gas.reaction, T, C, S0, h_mole)
    # Ydot = wdot / ρ_mass .* gas.MW
    # Tdot = -dot(h_mole, wdot) / ρ_mass / cp_mass
    return T
end
u0 = ones(41)
@time dudt(u0)
@time dudt(u0)
@time du0 = ForwardDiff.gradient(dudt, u0)
@time du0 = ForwardDiff.gradient(dudt, u0)

@time du0 = ForwardDiff.gradient(x -> dudt(x), u0)
@time du0 = ForwardDiff.gradient(x -> dudt(x), u0)
