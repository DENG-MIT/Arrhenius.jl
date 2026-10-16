using Arrhenius, Test, ForwardDiff

function thermo_fixture(models)
    return Dict{String,Any}("phases"=>[Dict("species"=>["species_$i" for i in eachindex(models)])],
        "species"=>[Dict("name"=>"species_$i","thermo"=>model) for (i,model) in enumerate(models)])
end

@testset "species thermo regions and reference states" begin
    nasa9 = Dict("model"=>"NASA9","temperature-ranges"=>[200.,1000.,6000.,20000.],
        "data"=>[[0.,0.,cp,0.,0.,0.,0.,0.,0.] for cp in (2.5,3.5,4.5)])
    thermo = IdealGasThermo(thermo_fixture([nasa9]))
    for (T,expected) in ((999.999,2.5),(1000.,3.5),(1000.001,3.5),
                         (5999.999,3.5),(6000.,4.5),(6000.001,4.5))
        result = species_thermo(thermo,T)
        @test result.cp_R == [expected]
        @test result.h_RT == [expected]
        @test result.s_R ≈ [expected*log(T)]
    end
    shomate = Dict("model"=>"Shomate","temperature-ranges"=>[200.,1000.,6000.],
        "data"=>[[cp*R/1000,0.,0.,0.,0.,0.,0.] for cp in (2.5,3.5)])
    thermo = IdealGasThermo(thermo_fixture([shomate]))
    @test species_thermo(thermo,1000.).cp_R ≈ [2.5]
    @test species_thermo(thermo,1000.001).cp_R ≈ [3.5]

    nasa7 = Dict("model"=>"NASA7","temperature-ranges"=>[200.,1000.,6000.],
        "data"=>[[cp,0.,0.,0.,0.,0.,0.] for cp in (2.5,3.5)])
    thermo = IdealGasThermo(thermo_fixture([nasa7]))
    @test isnothing(thermo.extra)
    @test species_thermo(thermo,1000.).cp_R == [2.5]
    @test species_thermo(thermo,1000.001).cp_R == [3.5]
    original = species_thermo(thermo,1500.;P=3one_atm)
    reference = deepcopy(nasa7)
    reference["reference-pressure"] = "1 bar"
    converted = species_thermo(IdealGasThermo(thermo_fixture([reference])),1500.;P=3one_atm)
    @test converted.cp_R == original.cp_R
    @test converted.h_RT == original.h_RT
    @test converted.s_R-original.s_R ≈ [log(1e5/one_atm)] atol=1e-14
    # The original four-field constructor stays a NASA7-only representation.
    copied = IdealGasThermo(copy(thermo.nasa_low),copy(thermo.nasa_high),copy(thermo.Trange),thermo.isTcommon)
    @test isnothing(copied.extra)
    @test species_thermo(copied,1500.) == species_thermo(thermo,1500.)
    @test isnothing(Arrhenius._convert_precision(copied,Float32).extra)
    @test IdealGasThermo{Float32}(thermo.nasa_low,thermo.nasa_high,thermo.Trange,thermo.isTcommon) isa IdealGasThermo{Float32}

    constant = Dict("model"=>"constant-cp","T0"=>"1000 K","h0"=>"9.22 kcal/mol",
        "s0"=>"-3.02 cal/mol/K","cp0"=>"5.95 cal/(mol*K)")
    thermo = IdealGasThermo(thermo_fixture([constant]))
    result = species_thermo(thermo,1000.)
    @test result.cp_R ≈ [5.95*4184/R]
    @test result.h_RT ≈ [9.22*4184/R]
    @test result.s_R ≈ [-3.02*4184/R]
    result = species_thermo(thermo,2000.)
    @test result.h_RT ≈ [(9.22*4184000+5.95*4184*1000)/(R*2000)]
    @test result.s_R ≈ [(-3.02+5.95*log(2))*4184/R]
    @test thermo.Trange[1,1] == 0
    @test thermo.Trange[1,end] == Inf
    defaults = IdealGasThermo(thermo_fixture([Dict("model"=>"constant-cp")]))
    @test species_thermo(defaults,700.).cp_R == [0.]
    @test species_thermo(defaults,700.).h_RT == [0.]
    @test species_thermo(defaults,700.).s_R == [0.]

    # Root and nested unit defaults must agree with explicit unit strings.
    implicit = Dict("model"=>"constant-cp","T0"=>1000.,"h0"=>9220.,"s0"=>-3.02,"cp0"=>5.95)
    fixture = thermo_fixture([implicit])
    fixture["units"] = Dict("energy"=>"cal","quantity"=>"mol")
    @test species_thermo(IdealGasThermo(fixture),2000.) == result
    typed = Arrhenius._convert_precision(thermo,Float32)
    @test eltype(typed.extra.coefficients[1]) == Float32
    @test typed.extra.coefficients[1] !== thermo.extra.coefficients[1]
    @test species_thermo(typed,1000.f0;P=Float32(one_atm)).cp_R ≈ species_thermo(thermo,1000.).cp_R rtol=1e-6
    @test ForwardDiff.derivative(T -> species_thermo(thermo,T).h_RT[1]*R*T,1500.) ≈ 5.95*4184
end

@testset "species thermo input rejection" begin
    good = Dict("model"=>"NASA9","temperature-ranges"=>[200.,1000.],
        "data"=>[[0.,0.,3.5,0.,0.,0.,0.,0.,0.]])
    for (key,value) in (("model","unimplemented"),("temperature-ranges",[1000.,200.]),
                        ("temperature-ranges",[200.,200.]),("temperature-ranges",[0.,1000.]),
                        ("temperature-ranges",[200.,Inf]),("temperature-ranges",[200.,500.,1000.]),
                        ("data",[zeros(7)]),("data",[fill(NaN,9)]),
                        ("reference-pressure",-1.),("reference-pressure","one atm"))
        candidate = deepcopy(good)
        candidate[key] = value
        @test_throws ArgumentError IdealGasThermo(thermo_fixture([candidate]))
    end
    for (key,value) in (("T0",0.),("T0","100 C"),("h0",NaN),("cp0","1 bad/mol/K"),
                        ("cp0","1 J/kg/K"),("T-min",-1.),("T-max",-1.))
        candidate = Dict{String,Any}("model"=>"constant-cp",key=>value)
        @test_throws ArgumentError IdealGasThermo(thermo_fixture([candidate]))
    end
    thermo = IdealGasThermo(thermo_fixture([good]))
    @test_throws ArgumentError species_thermo(thermo,0.)
    @test_throws ArgumentError species_thermo(thermo,Inf)
    @test_throws ArgumentError species_thermo(thermo,300.;P=0.)
end
