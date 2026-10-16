module PlasmaEnergyReactorTests
using Arrhenius, Test, YAML, LinearAlgebra
using ..PlasmaThermochemistryTests: fixture

@testset "Closed plasma energy reactor and signed states" begin
    mktempdir() do directory
        path, root = fixture(directory)
        # A negative enthalpy reference is physically valid and is not a
        # temperature/domain constraint.
        for species in root["species"]
            species["thermo"]["data"][1][6] = -10000.0
        end
        YAML.write_file(path, root)
        m = PlasmaMechanism(path)
        s = PlasmaState(m)
        @test_throws ArgumentError PlasmaEnergyReactor(s)
        update_eedf!(s)
        set_reduced_electric_field!(s, 2e-21)
        r = PlasmaEnergyReactor(s; volume=0.25)
        rhs = reactor_rhs(r)
        u = reactor_state(r)
        saved = copy(u)
        props = reactor_properties(rhs, u)
        thermo = plasma_thermodynamics(s)
        @test length(u) == m.n_species + 2
        @test u[1] ≈ s.density * 0.25
        @test u[2] < 0
        @test u[2] ≈ u[1] * thermo.h_mass
        @test u[3:end] == s.mass_fractions
        @test props.Y == s.mass_fractions
        @test props.T ≈ s.temperature rtol=1e-11
        @test props.Te == plasma_properties(s).Te
        @test props.rho ≈ s.density rtol=1e-11
        @test props.volume ≈ 0.25
        @test props.h_mass ≈ u[2]/u[1] rtol=1e-11
        @test u == saved
        du = similar(u)
        rhs(du, u, nothing, 0.0)
        @test du[1] == 0
        @test du[3:end] ≈ plasma_rates(s).dYdt rtol=1e-11
        @test du[2] ≈ plasma_properties(s).joule_heating * 0.25 rtol=1e-11
        @test abs(sum(du[3:end])) <= 1e-12 * sum(abs,du[3:end])
        @test u == saved

        # Constant gas-temperature derivatives exclude electron Cp because
        # electron temperature is an independent stored quantity.
        amount = u[3:end] ./ m.MW
        capacity = dot(amount, thermo.partial_molar_heat_capacities) -
                   amount[m.electron_index] * thermo.partial_molar_heat_capacities[m.electron_index]
        deltaH = 1e-5 * abs(u[2])
        plus = copy(u); minus = copy(u)
        plus[2] += deltaH; minus[2] -= deltaH
        dTdH = (reactor_properties(rhs,plus).T-reactor_properties(rhs,minus).T)/(2deltaH)
        @test dTdH ≈ inv(u[1]*capacity) rtol=1e-8

        # Signed nonlinear solver trial states must not be normalized or
        # clipped; only the Joule source has the physical zero branch.
        e = m.electron_index + 2
        trial = copy(u); trial[e] = -abs(trial[e])
        negative = similar(u); rhs(negative,trial,nothing,0.0)
        @test all(isfinite,negative)
        @test negative[2] == 0
        @test reactor_properties(rhs,trial).Y == trial[3:end]
        @test trial[e] < 0
        trial[e] = 0
        rhs(negative,trial,nothing,0.0)
        @test negative[2] == 0

        # Independent RHS workspaces and externally cached distributions.
        other = reactor_rhs(r)
        other_du = similar(u)
        other(other_du,u,nothing,0.0)
        @test other_du == du
        rhs(negative,trial,nothing,0.0)
        other(other_du,u,nothing,0.0)
        @test other_du == du
        @test s.mass_fractions == saved[3:end]
        @test u == saved

        J = zeros(length(u),length(u))
        reactor_jacobian!(J,u,rhs)
        @test all(isfinite,J)
        @test all(iszero,view(J,1,:))
        @test J[2,1] ≈ du[2]/u[1] rtol=1e-11
        @test J[2,e] ≈ du[2]/u[e] rtol=1e-11
        @test J[2,2] == 0
        @test u == saved

        problem = reactor_problem(r,(0.0,1e-9))
        @test problem.u0 == u && problem.p === nothing
        @test problem.f !== rhs
        time_derivative = fill(NaN,length(u))
        problem.tgrad(time_derivative,u,nothing,0.0)
        @test all(iszero,time_derivative)
        @test_throws ArgumentError reactor_problem(r,(1.0,0.0))
        @test_throws ArgumentError PlasmaEnergyReactor(s;volume=0.0)
        invalid = copy(u); invalid[1] = 0
        @test_throws DomainError rhs(negative,invalid,nothing,0.0)
        @test_throws DimensionMismatch rhs(zeros(2),u,nothing,0.0)
    end
end

@testset "Signed accepted plasma targets retain strict standalone EEDF inputs" begin
    mktempdir() do directory
        path, root = fixture(directory)
        push!(root["collisions"], Dict("equation"=>"H + e => H + e",
            "type"=>"electron-collision-plasma", "energy-levels"=>[0.,10.],
            "cross-sections"=>[1e-20,1e-20]))
        YAML.write_file(path,root)
        m = PlasmaMechanism(path)
        s = PlasmaState(m)
        # Accepted solver states can contain tiny signed trace targets. The
        # standalone input contract remains strict; the plasma bridge retains
        # their signed weights as in the source collision operator.
        target = findfirst(==("H2"),m.species_names)
        s.mass_fractions[target] = -1e-40
        set_plasma_state!(s;pressure=101325.)
        saved = copy(s.mass_fractions)
        p = plasma_properties(s)
        fractions = Dict(zip(m.species_names,p.X))
        weights = Dict(zip(m.species_names,m.MW))
        @test_throws ArgumentError EEDFState(m.thermal.eedf_model;
            T=p.T,P=p.P,mole_fractions=fractions,molecular_weights=weights,reduced_field=0.)
        unchecked = EEDFState(p.T,p.P,fractions,weights,0.,0.)
        @test_throws ArgumentError solve_eedf(m.thermal.eedf_model,unchecked)
        update_eedf!(s)
        @test s.eedf.converged
        @test s.mass_fractions == saved
        @test s.mass_fractions[target] < 0
        @test plasma_properties(s).Te == p.Te
        r = PlasmaEnergyReactor(s)
        rhs = reactor_rhs(r)
        u = reactor_state(r)
        @test update_eedf!(rhs,u;reduced_field=2e-21) === rhs
        before = reactor_properties(rhs,u)
        @test update_eedf!(rhs,u) === rhs
        after = reactor_properties(rhs,u)
        @test after.electric_field == before.electric_field
        @test after.Y == saved
        @test after.Te == p.Te
    end
end
end
