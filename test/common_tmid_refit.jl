using Arrhenius, Test, YAML

@testset "common NASA7 midpoint refit" begin
    mktempdir() do dir
        path=joinpath(dir,"source.yaml"); output=joinpath(dir,"common.yaml")
        species(name,elem,mid,low,high)=Dict("name"=>name,"composition"=>Dict(elem=>1),
            "thermo"=>Dict("model"=>"NASA7","temperature-ranges"=>[250.0,mid,4000.0],"data"=>[low,high]))
        single=Dict("name"=>"C","composition"=>Dict("C"=>1),"thermo"=>Dict(
            "model"=>"NASA7","temperature-ranges"=>[300.0,5000.0],"data"=>[[3.5,0.,0.,0.,0.,0.,1.]]))
        YAML.write_file(path,Dict("phases"=>[Dict("name"=>"gas","thermo"=>"ideal-gas",
            "elements"=>["H","O","C"],"species"=>["H","O","C"],"reactions"=>"all")],
            "species"=>[species("H","H",900.0,[3.5,0.,0.,0.,0.,0.,1.],[4.2,0.,0.,0.,0.,0.,2.]),
                        species("O","O",1600.0,[3.3,0.,0.,0.,0.,0.,1.],[4.4,0.,0.,0.,0.,0.,2.]),single],
            "reactions"=>Any[]))
        refused,_=refit_nasa7_common_tmid(path,output;Tmid=1200.0,samples=32)
        @test refused.species_changed==1
        @test refused.species_refused==2
        @test count(c->!c.accepted,refused.changes)==2
        @test isempty(refused.skipped_species)
        original=YAML.load_file(path); refused_yaml=YAML.load_file(output)
        @test refused_yaml["species"][1]["thermo"]["temperature-ranges"] == original["species"][1]["thermo"]["temperature-ranges"]
        report,_=refit_nasa7_common_tmid(path,output;Tmid=1200.0,samples=32,force=true)
        @test report.species_changed==3
        @test report.species_refused==0
        @test all(c.accepted for c in report.changes)
        @test report.shared_interval==(300.0,4000.0)
        out=YAML.load_file(output)
        for sp in out["species"]
            th=sp["thermo"]
            @test th["temperature-ranges"][2]==1200.0
            @test length(th["data"])==2
        end
        after=audit_mechanism(output)
        @test !any(x->hasproperty(x,:check)&&x.check==:nasa_continuity,after.issues)
        @test all(isfinite(getproperty(c,k)) for c in report.changes for k in
            (:max_abs_delta_cp_R,:max_abs_delta_h_RT,:max_abs_delta_s_R,:max_abs_delta_g_RT))
        @test_throws ArgumentError refit_nasa7_common_tmid(path,output;Tmid=300.0)
    end
end