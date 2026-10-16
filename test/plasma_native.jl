module PlasmaNativeTests
using Arrhenius, Test, YAML, LinearAlgebra

function fixture(directory)
    species = [Dict("name"=>name, "composition"=>composition) for (name,composition) in
        [("O",Dict("O"=>1)), ("O2",Dict("O"=>2)), ("O2-",Dict("O"=>2,"E"=>1)),
         ("O-",Dict("O"=>1,"E"=>1)), ("O2+",Dict("O"=>2,"E"=>-1))]]
    YAML.write_file(joinpath(directory,"ions.yaml"),Dict("species"=>species))
    phase = Dict("name"=>"test-plasma", "thermo"=>"plasma", "elements"=>["O","E"],
        "species"=>[Dict("species"=>["e"]),Dict("ions.yaml/species"=>["O","O2","O2-","O-","O2+"])],
        "reactions"=>[Dict("reactions"=>"all"),Dict("collisions"=>"all")],
        "state"=>Dict("T"=>300.0,"P"=>"0.01 atm","X"=>Dict("O2"=>1.,"e"=>.005,"O2+"=>.005)),
        "electron-energy-distribution"=>Dict("type"=>"isotropic","shape-factor"=>1.0,
            "mean-electron-energy"=>"2 eV","energy-levels"=>[0.,1.,2.,3.,4.]))
    root = Dict("units"=>Dict("length"=>"cm","quantity"=>"molec","activation-energy"=>"K"),
        "phases"=>[phase],"species"=>[Dict("name"=>"e","composition"=>Dict("E"=>1))],
        "reactions"=>[
            Dict("equation"=>"O2+ + e => O + O","type"=>"two-temperature-plasma",
                 "rate-constant"=>Dict("A"=>6e-5,"b"=>-1.,"Ea-gas"=>0.,"Ea-electron"=>0.,"b-gas"=>.2,"T-inv"=>500.)),
            Dict("equation"=>"E + O2 + O2 => O2- + O2","type"=>"two-temperature-plasma",
                 "rate-constant"=>Dict("A"=>4.2e-27,"b"=>-1.,"Ea-gas"=>600.,"Ea-electron"=>700.)),
            Dict("equation"=>"O2- + O2 => O2 + e + O2",
                 "rate-constant"=>Dict("A"=>1.559e-11,"b"=>.5,"Ea"=>5590.)),
            Dict("equation"=>"O2- + O2+ + M => O2 + O2 + M","type"=>"three-body",
                 "rate-constant"=>Dict("A"=>3.118e-19,"b"=>-2.5,"Ea"=>0.),
                 "efficiencies"=>Dict("O2"=>2.0))],
        "collisions"=>[Dict("equation"=>"O2 + e => O- + O","type"=>"electron-collision-plasma",
            "energy-levels"=>[0.,4.],"cross-sections"=>[1e-20,1e-20])])
    path = joinpath(directory,"plasma.yaml");YAML.write_file(path,root)
    return path, YAML.load_file(path)
end

@testset "native plasma import and state" begin
    mktempdir() do directory
        path,root = fixture(directory)
        m = PlasmaMechanism(path)
        @test m.species_names == ["e","O","O2","O2-","O-","O2+"]
        @test m.n_species == 6 && m.n_reactions == 5
        electron_mass = 9.1093837015e-31 * 6.02214076e26
        @test m.MW ≈ [electron_mass,15.999,31.998,31.998+electron_mass,15.999+electron_mass,31.998-electron_mass] rtol=1e-14
        @test m.elemental_matrix * m.stoichiometry == zeros(2,5)
        @test m.thirdbody == [false,true,true,true,false]
        @test m.reactants[3,2] == 1 && m.products[3,2] == 0
        @test m.reactants[3,3] == 0 && m.products[3,3] == 1
        @test m.efficiencies[:,2] == [0,0,1,0,0,0]
        @test m.efficiencies[:,4] == [1,1,2,1,1,1]
        conversion = 6.02214076e26 / 1e6
        @test m.rate_parameters[1,1] ≈ 6e-5 * conversion rtol=1e-14
        @test m.rate_parameters[1,2] ≈ 4.2e-27 * conversion^2 rtol=1e-14
        @test m.rate_parameters[1,3] ≈ 1.559e-11 * conversion rtol=1e-14
        @test m.rate_parameters[3,2] ≈ 600R
        @test m.rate_parameters[4,2] ≈ 700R
        @test m.initial_temperature == 300.0
        @test m.initial_pressure ≈ 1013.25
        @test m.initial_mean_electron_energy == 2.0
        @test m.initial_mole_fractions ≈ [.005,0,1,0,0,.005]/1.01
        s = PlasmaState(m)
        before = plasma_properties(s)
        @test before.X ≈ m.initial_mole_fractions
        @test before.P ≈ 1013.25 rtol=1e-14
        set_mean_electron_energy!(s,10.)
        after = plasma_properties(s)
        @test after.rho == before.rho
        @test after.Y == before.Y && after.T == before.T
        @test after.Te ≈ 5before.Te
        @test after.P / before.P ≈ (after.T + after.X[1]*(after.Te-after.T))/(before.T+before.X[1]*(before.Te-before.T)) rtol=1e-14
        set_plasma_state!(s; pressure=1013.25)
        @test plasma_properties(s).P ≈ 1013.25 rtol=1e-14
        @test s.density < after.rho
        @test s.mean_electron_energy == 10.
        saved = plasma_properties(s)
        @test_throws ArgumentError set_mean_electron_energy!(s,0.)
        @test_throws ArgumentError set_plasma_state!(s;temperature=-1.)
        @test_throws ArgumentError set_plasma_state!(s;mass_fractions=[1.,0,0,0,-1,0])
        @test plasma_properties(s) == saved
        @test_throws ArgumentError PlasmaState(m;mole_fractions=Dict("missing"=>1.))
        @test_throws ArgumentError PlasmaState(m;mole_fractions=ones(6),mass_fractions=ones(6))
        @test PlasmaMechanism(path;atomic_weights=Dict("O"=>18.0)).MW[3] == 36.0
        @test_throws ArgumentError PlasmaMechanism(path;phase="missing")
        root["reactions"][1]["rate-constant"]["Ea-gas"]="2 kJ/mol"
        root["reactions"][1]["rate-constant"]["Ea-electron"]="3 cal/mol"
        YAML.write_file(path,root)
        explicit_units=PlasmaMechanism(path)
        @test explicit_units.rate_parameters[3,1] == 2e6
        @test explicit_units.rate_parameters[4,1] ≈ 3*4184.0 rtol=1e-14
        root["phases"][1]["electron-energy-distribution"]["energy-levels"]=[0.,2.,1.]
        YAML.write_file(path,root)
        @test_throws ArgumentError PlasmaMechanism(path)
        path,root=fixture(directory);root["reactions"][1]["equation"]="O2+ + e <=> O + O";YAML.write_file(path,root)
        @test_throws ArgumentError PlasmaMechanism(path)
        path,root=fixture(directory);root["reactions"][1]["orders"]=Dict("e"=>-.5);YAML.write_file(path,root)
        @test_throws ArgumentError PlasmaMechanism(path)
    end
end

@testset "native plasma rates and species-only reactor" begin
    mktempdir() do directory
        path,_ = fixture(directory);m=PlasmaMechanism(path);s=PlasmaState(m)
        set_mean_electron_energy!(s,10.)
        values=plasma_rates(s);properties=plasma_properties(s)
        # Maxwellian shape1 has an independent closed-form distribution.
        e=m.energy_levels;expected=(2/sqrt(pi))*(1.5/10)^1.5 .* exp.(-1.5 .*e./10)
        @test values.electron_energy_distribution ≈ expected rtol=1e-13
        integrand=e.*expected.*1e-20
        integral=(integrand[1]+4integrand[2]+2integrand[3]+4integrand[4]+integrand[5])/3
        expected_collision=sqrt(2*1.602176634e-19/9.1093837015e-31)*6.02214076e26*integral
        @test values.forward_rate_constants[5] ≈ expected_collision rtol=1e-13
        C=values.concentrations;kf=values.forward_rate_constants
        @test kf[1] ≈ (6e-5*6.02214076e20)*300.0^.2/properties.Te*exp(-300.0/500.0) rtol=1e-13
        @test values.net_rates_of_progress ≈ [kf[1]*C[6]*C[1],kf[2]*C[1]*C[3]^2,
            kf[3]*C[4]*C[3],kf[4]*C[4]*C[6]*(sum(C)+C[3]),kf[5]*C[3]*C[1]] rtol=1e-13
        r=PlasmaReactor(s);u=reactor_state(r);rhs=reactor_rhs(r);du=similar(u);rhs(du,u,nothing,0.)
        @test length(u)==m.n_species
        @test du ≈ values.dYdt rtol=2e-13
        @test abs(sum(du)) <= 2e-13*sum(abs,du)
        @test norm(m.elemental_matrix*(du./m.MW),Inf) <= 2e-13*norm(du./m.MW,1)
        snapshot=reactor_state(r);set_plasma_state!(s;temperature=900.,mole_fractions=Dict("O2"=>1.))
        @test reactor_state(r)==snapshot
        @test r.temperature==300.
        @test reactor_properties(r).P ≈ properties.P rtol=1e-14
        problem=reactor_problem(r,(0.,1e-6));J=zeros(6,6)
        problem.jac(J,u,nothing,0.)
        @test problem.jac.calls[]==1
        @test all(isfinite,J)
        for signed in (false,true)
            y=copy(u)
            if signed;y[2]=-1e-14;y[5]=-1e-15;y[3]+=1.1e-14;end
            problem.jac(J,y,nothing,0.)
            for j in eachindex(y)
                # Higher-precision differencing avoids cancellation when a zero
                # species perturbs a large production rate by a tiny amount.
                h=big"1e-5"*max(abs(BigFloat(y[j])),big"1e-7")
                yp=BigFloat.(y);ym=copy(yp);yp[j]+=h;ym[j]-=h
                dp=similar(yp);dm=similar(ym);rhs(dp,yp);rhs(dm,ym)
                numerical=(dp.-dm)./(2h)
                @test norm(J[:,j]-numerical,Inf)<=1e-5*max(norm(J[:,j],Inf),1.)
            end
        end
        other=reactor_problem(r,(0.,1e-6))
        @test other.f.C !== problem.f.C && other.f.q !== problem.f.q
        @test other.jac.config !== problem.jac.config
        @test_throws DimensionMismatch rhs(zeros(5),zeros(5))
        @test_throws ArgumentError reactor_problem(r,(0.,0.))
        @test solve_reactor(r,(0.,1e-6);integrator=(p;kwargs...)->length(p.u0))==6
    end
end
end
