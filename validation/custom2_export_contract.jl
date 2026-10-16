using Test
include(joinpath(@__DIR__,"custom2_case.jl"))

function custom2_export_contract(output)
    @testset "custom2 real exporter packing without chemistry" begin
        @test isdefined(Main,:native_network_bdf)
        @test isdefined(Main,:run_moving_wall_example)
        cli=Module(:Custom2CLILoadingContract)
        Core.eval(cli,:(include(path)=Base.include($cli,path)))
        Base.include(cli,joinpath(@__DIR__,"..","example","reactors","custom2.jl"))
        @test isdefined(cli,:run_moving_wall_example)
        @test isdefined(cli,:native_network_bdf)
        times=[0.,.5];state=zeros(16,2)
        species=["species_$i" for i in 1:10];elements=["H","O"]
        scalars=NamedTuple{Tuple(Symbol.(CUSTOM2_SCALARS))}(ntuple(_->1.,length(CUSTOM2_SCALARS)))
        row=merge(scalars,(Y=fill(.1,10),X=fill(.1,10),species_mass=fill(.1,10),elements=[1.,2.]))
        rows=[row,row]
        data=custom2_pack(times,state,species,elements,ones(10),rows)
        @test data["state_full"]===state
        @test data["time"]==times
        @test String(copy(data["species_names_utf8"]))==join(species,'\n')
        for key in CUSTOM2_SCALARS;@test size(data[key])==(2,);end
        for key in ("Y","X","species_mass");@test size(data[key])==(2,10);end
        @test data["elements"]==[1. 2.;1. 2.]
        @test_throws DimensionMismatch custom2_pack(times,state[:,1:1],species,elements,ones(10),rows)
        @test_throws DimensionMismatch custom2_pack(times,zeros(15,2),species,elements,ones(10),rows)
        @test_throws DimensionMismatch custom2_pack(times,state,species,elements,ones(9),rows)
        custom2_write(output,data,Dict{String,Any}("completed"=>false,"fixture"=>true))
        saved=NPZ.npzread(joinpath(output,"native.npz"))
        @test Set(keys(saved))==Set(keys(data))
        for key in keys(data);@test saved[key]==data[key];end
        meta=TOML.parsefile(joinpath(output,"native.toml"))
        @test meta["native_sha256"]==bytes2hex(sha256(read(joinpath(output,"native.npz"))))
    end
end

length(ARGS)==1 || error("supply a new fixture directory")
custom2_export_contract(abspath(ARGS[1]))
