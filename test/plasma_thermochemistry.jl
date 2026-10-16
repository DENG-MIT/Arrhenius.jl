module PlasmaThermochemistryTests
using Arrhenius, Test, YAML, LinearAlgebra

function fixture(directory)
    specs = Any[]
    for (name,composition,cp) in [
            ("H2",Dict("H"=>2),3.5),("H",Dict("H"=>1),2.5),
            ("H+",Dict("H"=>1,"E"=>-1),2.5),("e",Dict("E"=>1),2.5)]
        push!(specs,Dict("name"=>name,"composition"=>composition,
            "thermo"=>Dict("model"=>"NASA7","temperature-ranges"=>[200.,6000.],
                "data"=>[[cp,0.,0.,0.,0.,-100.,1.]])))
    end
    YAML.write_file(joinpath(directory,"gas.yaml"),Dict(
        "units"=>Dict("length"=>"cm","quantity"=>"mol","activation-energy"=>"cal/mol"),
        "reactions"=>[Dict("equation"=>"H + H <=> H2",
            "rate-constant"=>Dict("A"=>1e3,"b"=>0.,"Ea"=>1.))]))
    phase=Dict("name"=>"gas-plasma","thermo"=>"plasma","elements"=>["H","E"],
        "species"=>["H2","H","H+","e"],
        "reactions"=>[Dict("gas.yaml/reactions"=>"all"),Dict("reactions"=>"all"),Dict("collisions"=>"all")],
        "electron-energy-distribution"=>Dict("type"=>"Boltzmann-two-term","energy-levels"=>collect(0.:.25:10.)),
        "state"=>Dict("T"=>300.,"P"=>101325.,"X"=>Dict("H2"=>.7,"H"=>.298,"H+"=>.001,"e"=>.001)))
    root=Dict("units"=>Dict("length"=>"m","quantity"=>"kmol","activation-energy"=>"J/kmol"),
        "phases"=>[phase],"species"=>specs,"reactions"=>Any[
        Dict("equation"=>"H + H + M => H2 + M","type"=>"three-body",
             "rate-constant"=>Dict("A"=>2e8,"b"=>0.,"Ea"=>0.)),
        Dict("equation"=>"H + H (+ M) => H2 (+ M)","type"=>"falloff",
             "high-P-rate-constant"=>Dict("A"=>3e8,"b"=>0.,"Ea"=>0.),
             "low-P-rate-constant"=>Dict("A"=>5e10,"b"=>0.,"Ea"=>0.),
             "Troe"=>Dict("A"=>.5,"T1"=>1000.,"T2"=>10000.,"T3"=>100.)),
        Dict("equation"=>"H2 => H + H","type"=>"Chebyshev",
             "temperature-range"=>[200.,6000.],"pressure-range"=>["1 atm","10 atm"],
             "data"=>[[1.,.1],[.2,.05]]),
        Dict("equation"=>"H+ + e + M => H + M","type"=>"two-temperature-plasma",
             "rate-constant"=>Dict("A"=>1e8,"b"=>-1.5,"Ea-gas"=>0.,"Ea-electron"=>0.)),
        Dict("equation"=>"H2 => H + H","type"=>"Chebyshev",
             "temperature-range"=>[200.,6000.],"pressure-range"=>["1 atm","1 atm"],
             "data"=>[[1.],[.2]])],
        "collisions"=>[Dict("equation"=>"H2 + e => H2 + e","type"=>"electron-collision-plasma",
            "energy-levels"=>[0.,10.],"cross-sections"=>[1e-20,1e-20])])
    path=joinpath(directory,"plasma.yaml");YAML.write_file(path,root)
    return path,root
end

@testset "Boltzmann thermochemistry and cached state" begin
    mktempdir() do directory
        path,root=fixture(directory)
        m=PlasmaMechanism(path);s=PlasmaState(m);p0=plasma_properties(s)
        @test m.n_species==4 && m.n_reactions==7
        @test m.thermal.reaction.index_three_body==[2,5]
        @test m.thermal.reaction.index_falloff==[3]
        @test m.thermal.reaction.Arrhenius_coeffs[1,:] ≈ [1.,0.,1.]
        @test p0.Te ≈ .001 rtol=1e-15
        @test_throws ArgumentError plasma_rates(s)
        @test_throws ArgumentError PlasmaReactor(s)
        update_eedf!(s); r0=plasma_rates(s);t0=plasma_thermodynamics(s)
        @test r0.forward_rate_constants[1] ≈ exp(-4184/(Arrhenius.R*300))
        @test r0.forward_rate_constants[2] == 2e8
        @test r0.forward_rate_constants[5] ≈ 1e8*p0.Te^(-1.5) rtol=1e-14
        @test r0.net_rates_of_progress[5] ≈ r0.forward_rate_constants[5]*
            r0.concentrations[3]*r0.concentrations[4]*sum(r0.concentrations) rtol=1e-14
        @test abs(sum(r0.dYdt)) <= 1e-12*sum(abs,r0.dYdt)
        saved=s.eedf
        set_reduced_electric_field!(s,2e-21);field=s.electric_field
        @test s.eedf === saved
        @test plasma_rates(s).forward_rate_constants == r0.forward_rate_constants
        set_plasma_state!(s;pressure=2p0.P)
        rp=plasma_rates(s)
        @test rp.forward_rate_constants[[2,5,6]] == r0.forward_rate_constants[[2,5,6]]
        @test rp.net_rates_of_progress[5] ≈ 8r0.net_rates_of_progress[5] rtol=1e-14
        @test plasma_properties(s).reduced_electric_field ≈ 1e-21 rtol=1e-14
        @test s.eedf === saved && s.electric_field == field
        set_plasma_state!(s;temperature=600.)
        th=plasma_thermodynamics(s)
        @test th.partial_molar_enthalpies[4] == t0.partial_molar_enthalpies[4]
        @test th.partial_molar_enthalpies[3]-t0.partial_molar_enthalpies[3] ≈ 2.5Arrhenius.R*300
        set_plasma_enthalpy!(s,t0.h_mass;pressure=p0.P)
        @test s.temperature ≈ 300 rtol=1e-11
        @test s.eedf === saved && s.electric_field == field
        @test plasma_properties(s).Te == p0.Te
        @test plasma_thermodynamics(s).h_mass ≈ t0.h_mass rtol=1e-11
        r=plasma_rates(s);r.electron_energy_distribution .= 0
        @test any(>(0),s.eedf.edge_eedf)
        other=PlasmaState(m)
        @test other.eedf === nothing && other.electric_field==0
        @test other.mass_fractions !== s.mass_fractions
        before=(s.temperature,s.density,copy(s.mass_fractions),s.eedf,s.electric_field)
        @test_throws ErrorException set_plasma_enthalpy!(s,t0.h_mass+1e6;maxiter=1)
        @test (s.temperature,s.density,s.mass_fractions,s.eedf,s.electric_field)==before
        @test_throws ArgumentError set_reduced_electric_field!(s,-1)
        # At both pressure and temperature endpoints, the 2x2 expansion
        # becomes a signed sum of its four coefficients.
        cheb=Arrhenius._plasma_chebyshev_rate
        @test cheb([1. .1;.2 .05],(200.,6000.),(1e5,1e6),200.,1e5) ≈ 10.0^(.75)
        @test cheb([1. .1;.2 .05],(200.,6000.),(1e5,1e6),6000.,1e6) ≈ 10.0^(1.35)
        root["reactions"][3]["pressure-range"]=["1 atm","1 atm"]
        YAML.write_file(path,root)
        @test_throws ArgumentError PlasmaMechanism(path)
    end
end
end
