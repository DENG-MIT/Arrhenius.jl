using Arrhenius, Test, YAML

function audit_species(name,comp; diameter=3.0,well_depth=100.0)
    Dict("name"=>name,"composition"=>comp,
         "thermo"=>Dict("model"=>"NASA7","temperature-ranges"=>[200.,3000.],"data"=>[[3.5,0.,0.,0.,0.,0.,1.]]),
         "transport"=>Dict("model"=>"gas","geometry"=>length(comp)==1 ? "atom" : "linear","well-depth"=>well_depth,"diameter"=>diameter,"polarizability"=>1.0))
end

@testset "gas mechanism audit and NASA7 continuity refit" begin
    mktempdir() do dir
        path=joinpath(dir,"discontinuous.yaml")
        YAML.write_file(path,Dict(
            "phases"=>[Dict("name"=>"gas","thermo"=>"ideal-gas","elements"=>["H"],"species"=>["H2"],"reactions"=>"all")],
            "species"=>[Dict("name"=>"H2","composition"=>Dict("H"=>2),"note"=>"annotation retained only in input","thermo"=>Dict("model"=>"NASA7","temperature-ranges"=>[200.,1000.,3000.],"data"=>[[3.5,0.,0.,0.,0.,0.,1.],[4.5,0.,0.,0.,0.,0.,2.]]))],
            "reactions"=>Any[]))
        before=audit_mechanism(path)
        @test before.species==1 && before.reactions==0
        @test count(x->hasproperty(x,:check)&&x.check==:nasa_continuity,before.issues)==1
        output=joinpath(dir,"refit.yaml")
        refused,_=refit_nasa7_continuity(path,output)
        @test refused.species_changed==0 && refused.species_refused==1
        @test only(refused.changes).accepted==false
        @test count(x->hasproperty(x,:check)&&x.check==:nasa_continuity,audit_mechanism(output).issues)==1
        report,_=refit_nasa7_continuity(path,output;force=true)
        @test report.species_changed==1 && report.species_refused==0
        @test all(isfinite(getproperty(only(report.changes),k)) for k in (:max_abs_delta_cp_R,:max_abs_delta_h_RT,:max_abs_delta_s_R,:max_abs_delta_g_RT))
        @test isfile(output)
        after=audit_mechanism(output)
        @test count(x->hasproperty(x,:check)&&x.check==:nasa_continuity,after.issues)==0
        parsed=YAML.load_file(output)
        @test length(parsed["species"])==1
        @test !haskey(parsed["species"][1],"note")
    end
    mktempdir() do dir
        path=joinpath(dir,"collision.yaml")
        YAML.write_file(path,Dict(
          "phases"=>[Dict("name"=>"gas","thermo"=>"ideal-gas","elements"=>["H"],"species"=>["H","H2"],"reactions"=>"all")],
          "species"=>[audit_species("H",Dict("H"=>1);diameter=2.0),audit_species("H2",Dict("H"=>2);diameter=30.0)],
          "reactions"=>[Dict("equation"=>"2 H => H2","rate-constant"=>Dict("A"=>1e50,"b"=>0.,"Ea"=>0.))]))
        result=audit_mechanism(path)
        @test count(x->x.check==:collision_limit_screen,result.issues)==1
        @test count(x->x.check==:transport_sanity,result.issues)>=1
        @test !any(x->x.level==:error,result.issues)
    end
end


@testset "PLOG log-space range and underflow audit" begin
    mktempdir() do dir
        path=joinpath(dir,"plog-underflow.yaml")
        species(name,n)=Dict("name"=>name,"composition"=>Dict("H"=>n),
            "thermo"=>Dict("model"=>"NASA7","temperature-ranges"=>[200.,3000.],"data"=>[[3.5,0.,0.,0.,0.,0.,1.]]))
        pentry(P)=Dict("P"=>P,"A"=>1e10,"b"=>0.,"Ea"=>1e6)
        rxn(eq)=Dict("equation"=>eq,"type"=>"pressure-dependent-Arrhenius",
            "rate-constants"=>[pentry("0.1 atm"),pentry("10 atm")])
        rxn2=rxn("H2 + H => H + H + H")
        rxn2["rate-constants"]=[pentry("0.1 atm"),merge(pentry("0.1 atm"),Dict("A"=>-1e9)),pentry("10 atm")]
        YAML.write_file(path,Dict("units"=>Dict("activation-energy"=>"cal/mol"),
            "phases"=>[Dict("name"=>"gas","thermo"=>"ideal-gas","elements"=>["H"],"species"=>["H2","H"],"reactions"=>"all")],
            "species"=>[species("H2",2),species("H",1)],
            "reactions"=>[rxn("H2 => H + H"),rxn2]))
        report=audit_mechanism(path)
        under=filter(x->x.check==:plog_underflow_risk,report.issues)
        ranges=filter(x->x.check==:plog_valid_range,report.issues)
        @test length(ranges)==2
        @test length(under)>=4
        @test minimum(x.temperature for x in under)==200.0
        @test all(x.log_rate < log(floatmin(Float64)) for x in under)
        @test !any(x->x.check==:plog_nonpositive_rate,report.issues)
        @test all(x.pressure_range==(10132.5,1013250.0) for x in ranges)
        output=joinpath(dir,"cheb.yaml")
        fit,_=refit_plog_to_chebyshev(path,output,1;temperature_range=(200.,3000.),
            temperature_order=2,pressure_order=1,samples=12,max_log10_error=1e-8)
        @test fit.accepted
        @test fit.max_log10_error < 1e-8
        fitted=YAML.load_file(output)["reactions"][1]
        @test fitted["type"]=="Chebyshev"
        @test_throws ArgumentError CreateSolution(output)
        @test YAML.load_file(path)["reactions"][1]["type"]=="pressure-dependent-Arrhenius"
        @test length(fitted["data"])==2
        @test all(length(row)==1 for row in fitted["data"])
        refused,_=refit_plog_to_chebyshev(path,joinpath(dir,"refused.yaml"),1;
            temperature_range=(200.,3000.),temperature_order=2,pressure_order=1,
            samples=12,max_log10_error=0.0)
        @test !refused.accepted
        @test YAML.load_file(joinpath(dir,"refused.yaml"))["reactions"][1]["type"]=="pressure-dependent-Arrhenius"
    end
end
