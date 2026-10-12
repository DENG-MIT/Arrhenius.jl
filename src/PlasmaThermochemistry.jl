"""
    set_reduced_electric_field!(state, reduced_field)

Set E/N in V·m² using the current total number density. The absolute electric
field is retained across subsequent gas-state changes. This does not update the
electron distribution or its cached mobility; call update_eedf! afterward.
"""
function set_reduced_electric_field!(s::PlasmaState, reduced_field)
    s.mechanism.thermal === nothing &&
        throw(ArgumentError("an electric field requires a Boltzmann plasma mechanism"))
    EN = Float64(reduced_field)
    isfinite(EN) && EN >= 0 || throw(ArgumentError("E/N must be finite and nonnegative"))
    N = s.density * sum(s.mass_fractions ./ s.mechanism.MW) * _EEDF_AVOGADRO_KMOL
    s.electric_field = EN * N
    return s
end

"""
    update_eedf!(state; options=TwoTermOptions())

Solve the Boltzmann distribution at the current gas state and electric field,
continuing from the previous distribution when appropriate. Cache the result
and mobility without changing the independently stored electron temperature.
"""
function update_eedf!(s::PlasmaState; options=TwoTermOptions())
    data = s.mechanism.thermal
    data === nothing && throw(ArgumentError("EEDF updates require a Boltzmann plasma mechanism"))
    p = plasma_properties(s)
    m = s.mechanism
    N = s.density * sum(s.mass_fractions ./ m.MW) * _EEDF_AVOGADRO_KMOL
    state = _accepted_eedf_state(data.eedf_model; T=p.T, P=p.P,
        mole_fractions=Dict(zip(m.species_names,p.X)),
        molecular_weights=Dict(zip(m.species_names,m.MW)),
        reduced_field=s.electric_field/N, number_density=N)
    result = _solve_accepted_eedf(data.eedf_model,state;options,initial=s.eedf)
    result.converged || throw(ErrorException("Boltzmann EEDF did not converge"))
    s.eedf = result
    return s
end

function _plasma_species_thermo!(h,cp,s0,m::PlasmaMechanism,T,Te)
    data = m.thermal
    data === nothing && throw(ArgumentError("this plasma mechanism has no gas thermodynamics"))
    thermo = data.thermo
    length(h) == length(cp) == length(s0) == m.n_species ||
        throw(DimensionMismatch("one thermochemistry entry per plasma species required"))
    @inbounds for i in eachindex(h)
        Ti = i == m.electron_index ? Te : T
        if thermo.extra === nothing
            a = _nasa7_coefficients(thermo,i,Ti)
            cp_R = a[i,1]+Ti*(a[i,2]+Ti*(a[i,3]+Ti*(a[i,4]+Ti*a[i,5])))
            h_RT = a[i,1]+Ti*(a[i,2]/2+Ti*(a[i,3]/3+Ti*(a[i,4]/4+Ti*a[i,5]/5)))+a[i,6]/Ti
            s_R = a[i,1]*log(Ti)+Ti*(a[i,2]+Ti*(a[i,3]/2+Ti*(a[i,4]/3+Ti*a[i,5]/4)))+a[i,7]
        else
            cp_R,h_RT,s_R = _extended_thermo(thermo,i,Ti)
        end
        h[i],cp[i],s0[i] = R*Ti*h_RT,R*cp_R,R*s_R
    end
    return nothing
end

function _plasma_species_thermo(m::PlasmaMechanism,T,Te)
    h,cp,s0 = zeros(m.n_species),zeros(m.n_species),zeros(m.n_species)
    _plasma_species_thermo!(h,cp,s0,m,T,Te)
    return h,cp,s0
end


"Reusable storage for signed thermal-plasma source evaluations."
mutable struct _PlasmaThermoRateWorkspace{T}
    X::Vector{T}
    C::Vector{T}
    h::Vector{T}
    cp::Vector{T}
    s0::Vector{T}
    wdot::Vector{T}
    collision_rates::Vector{T}
    collision_integrand::Vector{T}
    kinetics::KineticsWorkspace{T}
end

function _PlasmaThermoRateWorkspace(m::PlasmaMechanism,::Type{T}=Float64) where {T}
    data=m.thermal
    data===nothing && throw(ArgumentError("thermal plasma workspace requires gas thermodynamics"))
    return _PlasmaThermoRateWorkspace(zeros(T,m.n_species),zeros(T,m.n_species),
        zeros(T,m.n_species),zeros(T,m.n_species),zeros(T,m.n_species),
        zeros(T,m.n_species),zeros(T,m.n_reactions),zeros(T,length(m.energy_levels)),
        KineticsWorkspace(data.reaction,T))
end

function _plasma_collision_rates!(rates,integrand,m::PlasmaMechanism,eedf::EEDFResult)
    length(rates)==m.n_reactions || throw(DimensionMismatch("one collision-rate entry per reaction required"))
    length(integrand)==length(m.energy_levels)==length(eedf.edge_eedf) ||
        throw(DimensionMismatch("collision quadrature grid length"))
    fill!(rates,zero(eltype(rates)))
    f=eedf.edge_eedf
    @inbounds for j in eachindex(rates)
        m.rate_types[j]==0x03 || continue
        for i in eachindex(f)
            energy=m.energy_levels[i]
            sigma=_linear_interp_hold(energy,m.collision_energy[j],m.cross_sections[j])
            integrand[i]=energy*f[i]*sigma
        end
        rates[j]=_EEDF_GAMMA*_EEDF_AVOGADRO_KMOL*_eedf_simpson(integrand,m.energy_levels)
    end
    return rates
end

function _plasma_thermal_sources!(workspace::_PlasmaThermoRateWorkspace,m::PlasmaMechanism,
        T,Te,P,rho,Y)
    data=m.thermal
    data===nothing && throw(ArgumentError("thermal plasma rates require gas thermodynamics"))
    length(Y)==m.n_species || throw(DimensionMismatch("one mass fraction per plasma species required"))
    _plasma_species_thermo!(workspace.h,workspace.cp,workspace.s0,m,T,Te)
    inverse_mw=zero(T)
    @inbounds for i in eachindex(Y)
        amount=Y[i]/m.MW[i]
        inverse_mw+=amount
        workspace.X[i]=amount
        workspace.C[i]=rho*amount
    end
    @inbounds for i in eachindex(Y)
        workspace.X[i]/=inverse_mw
    end
    e=m.electron_index
    workspace.s0[e]=workspace.s0[e]*(Te/T)+R*(1-Te/T)*log(P/one_atm)
    _rate_factors!(data.reaction,T,workspace.C,workspace.s0,workspace.h,
        workspace.kinetics,data.reverse_plan;pressure=P)
    kf=workspace.kinetics.kf
    @inbounds for j in eachindex(kf)
        kind=m.rate_types[j]
        if kind==0x02
            A,b,Eg,Ee,bg,inverse_T=view(m.rate_parameters,:,j)
            kf[j]=A*exp(bg*log(T)+b*log(Te)-Eg/(R*T)+
                Ee*(Te-T)/(R*Te*T)-T*inverse_T)
        elseif kind==0x03
            kf[j]=workspace.collision_rates[j]
        end
    end
    for (i,j) in enumerate(data.chebyshev_indices)
        kf[j]=_plasma_chebyshev_rate(data.chebyshev_coefficients[i],
            data.chebyshev_temperature_ranges[i],data.chebyshev_pressure_ranges[i],T,P)
    end
    @inbounds for j in data.reaction.index_three_body
        if m.rate_types[j]!=0x01
            kf[j]*=dot(view(data.reaction.efficiencies_coeffs,:,j),workspace.C)
        end
    end
    _mass_action!(workspace.kinetics,data.reaction,workspace.C)
    mul!(workspace.wdot,data.reaction.vk,workspace.kinetics.rates_of_progress)
    return workspace
end

"""
    plasma_thermodynamics(state)

Return partial molar enthalpies (J/kmol), heat capacities (J/kmol/K), reference
entropies (J/kmol/K), and mixture mass-specific enthalpy and heat capacity.
Heavy species use the gas temperature; electrons use the stored electron
temperature. The mixture heat capacity includes each species' heat capacity.
"""
function plasma_thermodynamics(s::PlasmaState)
    m = s.mechanism
    h,cp,s0 = _plasma_species_thermo(m,s.temperature,
        _plasma_electron_temperature(s.mean_electron_energy))
    amounts = s.mass_fractions ./ m.MW
    return (partial_molar_enthalpies=h,partial_molar_heat_capacities=cp,
        reference_molar_entropies=s0,h_mass=dot(amounts,h),cp_mass=dot(amounts,cp))
end

"""
    set_plasma_enthalpy!(state, enthalpy; pressure, rtol=1e-12, maxiter=500)

Set mass-specific enthalpy in J/kg and pressure in Pa by solving for gas
temperature. Retain composition, electron energy, electric field, EEDF and
mobility. A failed inversion leaves the state unchanged.
"""
function set_plasma_enthalpy!(s::PlasmaState,enthalpy;
        pressure=plasma_properties(s).P,rtol=1e-12,maxiter=500)
    target = Float64(enthalpy)
    isfinite(target) || throw(ArgumentError("enthalpy must be finite"))
    P = _plasma_positive(pressure,"pressure")
    tolerance = _plasma_positive(rtol,"enthalpy tolerance")
    maxiter isa Integer && maxiter > 0 || throw(ArgumentError("maxiter must be a positive integer"))
    T = s.temperature
    Te = _plasma_electron_temperature(s.mean_electron_energy)
    amounts = s.mass_fractions ./ s.mechanism.MW
    for _ in 1:maxiter
        h,cp,_ = _plasma_species_thermo(s.mechanism,T,Te)
        error = target-dot(amounts,h)
        if abs(error) <= tolerance*max(abs(target),R*T*sum(amounts),1.0)
            return set_plasma_state!(s;temperature=T,pressure=P)
        end
        slope = dot(amounts,cp)
        isfinite(slope) && slope > 0 || throw(DomainError(slope,"nonpositive plasma heat capacity during enthalpy inversion"))
        step = clamp(error/slope,-100.0,100.0)
        T += max(step,-0.5T)
        isfinite(T) || throw(DomainError(T,"nonfinite temperature during enthalpy inversion"))
    end
    throw(ErrorException("plasma enthalpy inversion did not converge in $maxiter iterations"))
end

function _plasma_chebyshev_rate(a,trange,prange,T,P)
    nt,np = size(a)
    # A single pressure coefficient needs no pressure reduction, so equal
    # pressure bounds are valid in that case.
    xr = (2/T-inv(trange[1])-inv(trange[2]))/(inv(trange[2])-inv(trange[1]))
    yr = np == 1 ? 0.0 :
        (2log10(P)-log10(prange[1])-log10(prange[2]))/(log10(prange[2])-log10(prange[1]))
    tprev,tcur = 1.0,xr
    value = 0.0
    @inbounds for i in 1:nt
        ti = i == 1 ? 1.0 : tcur
        pprev,pcur = 1.0,yr
        for j in 1:np
            pj = j == 1 ? 1.0 : pcur
            value += a[i,j]*ti*pj
            j >= 2 && ((pprev,pcur) = (pcur,2yr*pcur-pprev))
        end
        i >= 2 && ((tprev,tcur) = (tcur,2xr*tcur-tprev))
    end
    return 10.0^value
end

function _plasma_thermal_rates(s::PlasmaState)
    m,data = s.mechanism,s.mechanism.thermal
    s.eedf === nothing && throw(ArgumentError("call update_eedf! before evaluating Boltzmann plasma rates"))
    T = s.temperature
    Te = _plasma_electron_temperature(s.mean_electron_energy)
    C = s.mass_fractions .* s.density ./ m.MW
    h,_,s0 = _plasma_species_thermo(m,T,Te)
    P = plasma_properties(s).P
    # Express electron chemical potential in the common gas-RT convention
    # used by the shared equilibrium-constant kernel.
    e = m.electron_index
    s0[e] = s0[e]*(Te/T)+R*(1-Te/T)*log(P/one_atm)
    workspace = KineticsWorkspace(data.reaction)
    _rate_factors!(data.reaction,T,C,s0,h,workspace,data.reverse_plan;pressure=P)
    kf = workspace.kf
    f = s.eedf.edge_eedf
    integrand = similar(f)
    @inbounds for j in eachindex(kf)
        kind = m.rate_types[j]
        if kind == 0x02
            A,b,Eg,Ee,bg,inverse_T = view(m.rate_parameters,:,j)
            kf[j] = A*exp(bg*log(T)+b*log(Te)-Eg/(R*T)+
                Ee*(Te-T)/(R*Te*T)-T*inverse_T)
        elseif kind == 0x03
            for i in eachindex(f)
                energy = m.energy_levels[i]
                sigma = _linear_interp_hold(energy,m.collision_energy[j],m.cross_sections[j])
                integrand[i] = energy*f[i]*sigma
            end
            kf[j] = _EEDF_GAMMA*_EEDF_AVOGADRO_KMOL*_eedf_simpson(integrand,m.energy_levels)
        end
    end
    for (i,j) in enumerate(data.chebyshev_indices)
        kf[j] = _plasma_chebyshev_rate(data.chebyshev_coefficients[i],
            data.chebyshev_temperature_ranges[i],data.chebyshev_pressure_ranges[i],T,P)
    end
    all(isfinite,kf) || throw(DomainError(kf,"nonfinite plasma rate constant"))
    reported_kf = copy(kf)
    # The internal factor includes M. Report raw coefficients even when M=0.
    for j in data.reaction.index_three_body
        if m.rate_types[j] == 0x01
            A,b,Ea = view(data.reaction.Arrhenius_coeffs,j,:)
            reported_kf[j] = A*exp(b*log(T)-Ea*4184/(R*T))
        else
            # Specialized factors replaced the dummy thermal row above.
            kf[j] *= dot(view(data.reaction.efficiencies_coeffs,:,j),C)
        end
    end
    _mass_action!(workspace,data.reaction,C)
    q = workspace.rates_of_progress
    wdot = data.reaction.vk*q
    return (forward_rate_constants=reported_kf,net_rates_of_progress=q,
        net_production_rates=wdot,concentrations=C,dYdt=wdot.*m.MW./s.density,
        electron_energy_distribution=copy(f),electron_energy_levels=copy(m.energy_levels))
end

export set_reduced_electric_field!, update_eedf!, plasma_thermodynamics, set_plasma_enthalpy!
