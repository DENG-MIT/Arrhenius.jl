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
        report,_=refit_nasa7_continuity(path,output)
        @test report.species_changed==1
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
