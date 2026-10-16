using Arrhenius, NPZ, Random, Test, LinearAlgebra, YAML

function native_reader_ordered_bits(x::Float64)
    bits = reinterpret(UInt64,x)
    (bits >> 63) == 1 ? ~bits : bits | (UInt64(1) << 63)
end
function native_reader_max_ulp(a,b)
    size(a)==size(b) || return typemax(Int)
    isempty(a) && return 0
    maximum(abs(Int128(native_reader_ordered_bits(Float64(x))) -
                Int128(native_reader_ordered_bits(Float64(y)))) for (x,y) in zip(a,b))
end

@testset "native gas YAML reader" begin
    root = normpath(joinpath(@__DIR__,".."))
    mechanisms = (("2S_CH4_CM2.yaml",6,2),("JP10skeletal.yaml",40,232),
        ("LiDryer.yaml",9,24),("gri30.yaml",53,325),("h2o2.yaml",9,28),
        ("nDodecane_Reitz.yaml",100,553),("nc7_ver3.1_mech.yaml",631,4846))
    for (file,ns,nr) in mechanisms
        gas=CreateSolution(joinpath(root,"mechanism",file))
        @test (gas.n_species,gas.n_reactions)==(ns,nr)
    end

    @test Arrhenius._GAS_ATOMIC_WEIGHTS["E"] == Arrhenius._ION_ELECTRON_MW
    @test Arrhenius._GAS_ATOMIC_WEIGHTS["Be"] == 9.0121831

    fixture=joinpath(@__DIR__,"data","native_gas_features.yaml")
    native=CreateSolution(fixture)
    sidecar=CreateSolution(fixture;kinetics_reader=:sidecar)
    raw=Arrhenius.npzread(fixture*".npz")
    parsed=Arrhenius._native_gas_arrays(YAML.load_file(fixture),YAML.load_file(fixture)["phases"][1],
        native.species_names;mechanism=fixture)
    exact_keys=("reactant_stoich_coeffs","product_stoich_coeffs","reactant_orders",
        "is_reversible","efficiencies_coeffs","index_three_body","index_falloff",
        "index_falloff_Troe","Plog_reaction_indices","Plog_collider_indices",
        "Plog_group_offsets","Plog_rate_offsets","Plog_LMR_offsets","Plog_LMR_species")
    for key in exact_keys
        if haskey(parsed,key) && haskey(raw,key)
            @test parsed[key] == raw[key]
        end
    end
    float_keys=("molecular_weights","Arrhenius_coeffs","Arrhenius_A0","Arrhenius_b0",
        "Arrhenius_Ea0","Troe_A","Troe_T1","Troe_T2","Troe_T3","Plog_pressures",
        "Plog_Arrhenius","Plog_LMR_eps","BlowersMasel_coefficients")
    for key in float_keys
        @test native_reader_max_ulp(parsed[key],raw[key]) <= 3
    end

    # Removing the fixture sidecar still permits direct gas-phase construction.
    mktempdir() do dir
        path=joinpath(dir,"native.yaml")
        cp(fixture,path)
        gas=CreateSolution(path)
        @test gas.n_species==native.n_species
        @test gas.n_reactions==native.n_reactions

        # A stale optional sidecar must not prevent native chemistry loading.
        cp(fixture*".npz",path*".npz")
        open(path,"a") do io; write(io,"\n"); end
        @test CreateSolution(path).n_reactions==native.n_reactions
        @test_throws ArgumentError CreateSolution(path;kinetics_reader=:sidecar)
        @test_throws ArgumentError CreateSolution(path;kinetics_reader=:unknown)
    end

    # These are Cantera's published sidecar values for compositions whose
    # phase element order differs from alphabetical order.
    mw_species=[
        Dict("name"=>"HONO","composition"=>Dict("H"=>1,"N"=>1,"O"=>2)),
        Dict("name"=>"NCO","composition"=>Dict("C"=>1,"N"=>1,"O"=>1)),
    ]
    mw_phase=Dict("elements"=>["H","O","N","C"])
    mw_root=Dict("species"=>mw_species,"elements"=>["C","H","O","N"])
    @test Arrhenius._gas_molecular_weight(mw_root,mw_phase,"HONO") == 47.013
    @test Arrhenius._gas_molecular_weight(mw_root,mw_phase,"NCO") == 42.016999999999996

    Random.seed!(20261016)
    for T in (500.0,1200.0,2600.0)
        Y=rand(native.n_species); Y./=sum(Y)
        X=Y2X(native,Y,1/dot(Y,1 ./ native.MW))
        @test cal_cp_R(native,T,one_atm,X) == cal_cp_R(sidecar,T,one_atm,X)
        @test cal_h_RT(native,T,one_atm,X) == cal_h_RT(sidecar,T,one_atm,X)
        @test cal_s0_R(native,T,one_atm,X) == cal_s0_R(sidecar,T,one_atm,X)
        wn=set_states(native,T,one_atm,Y); ws=set_states(sidecar,T,one_atm,Y)
        scale=max(maximum(abs,ws;init=0.0),floatmin(Float64))
        @test maximum(abs.(wn.-ws);init=0.0)/scale <= 1e-13
    end

    for file in ("gri30.yaml","h2o2.yaml")
        path=joinpath(root,"mechanism",file)
        native=CreateSolution(path); sidecar=CreateSolution(path;kinetics_reader=:sidecar)
        arrays=Arrhenius._native_gas_arrays(YAML.load_file(path),YAML.load_file(path)["phases"][1],
            native.species_names;mechanism=path)
        data=Arrhenius.npzread(path*".npz")
        for key in exact_keys
            if haskey(arrays,key) && haskey(data,key)
                @test arrays[key] == data[key]
            end
        end
        for key in float_keys
            if haskey(arrays,key) && haskey(data,key)
                @test native_reader_max_ulp(arrays[key],data[key]) <= 2
            end
        end
        Y=fill(1/native.n_species,native.n_species)
        wn=set_states(native,1100.0,one_atm,Y); ws=set_states(sidecar,1100.0,one_atm,Y)
        scale=max(maximum(abs,ws;init=0.0),floatmin(Float64))
        @test maximum(abs.(wn.-ws);init=0.0)/scale <= 1e-13
    end
end
@testset "native gas reaction order validation" begin
    order_species=[Dict("name"=>"H2","composition"=>Dict("H"=>2)),
                   Dict("name"=>"H","composition"=>Dict("H"=>1)),
                   Dict("name"=>"AR","composition"=>Dict("Ar"=>1))]
    order_phase=Dict("name"=>"gas","species"=>["H2","H","AR"],"reactions"=>"all")
    order_root=Dict("species"=>order_species,"phases"=>[order_phase],"units"=>Dict())
    reaction(eq,fields=Dict())=merge(Dict("equation"=>eq,"rate-constant"=>Dict("A"=>1.0,"b"=>0.0,"Ea"=>0.0)),fields)
    function parse_order(r)
        root=merge(order_root,Dict("reactions"=>[r]))
        Arrhenius._native_gas_arrays(root,order_phase,["H2","H","AR"])
    end
    @test parse_order(reaction("H2 => H + H"))["reactant_orders"][1,1]==1
    @test parse_order(reaction("H2 => H + H",Dict("orders"=>Dict("H2"=>0.5))))["reactant_orders"][1,1]==0.5
    @test_throws ArgumentError parse_order(reaction("H2 => H + H",Dict("orders"=>Dict("AR"=>1.0))))
    @test parse_order(reaction("H2 => H + H",Dict("orders"=>Dict("AR"=>1.0),"nonreactant-orders"=>true)))["reactant_orders"][3,1]==1
    @test_throws ArgumentError parse_order(reaction("H2 => H + H",Dict("orders"=>Dict("H2"=>-0.5))))
    @test parse_order(reaction("H2 => H + H",Dict("orders"=>Dict("H2"=>-0.5),"negative-orders"=>true)))["reactant_orders"][1,1]==-0.5
    @test_throws ArgumentError parse_order(reaction("H2 <=> H + H",Dict("orders"=>Dict("H2"=>0.5))))
end

@testset "unsupported gas rate forms reject loudly" begin
    order_species=[Dict("name"=>"H2","composition"=>Dict("H"=>2)),
                   Dict("name"=>"H","composition"=>Dict("H"=>1)),
                   Dict("name"=>"AR","composition"=>Dict("Ar"=>1))]
    order_phase=Dict("name"=>"gas","species"=>["H2","H","AR"],"reactions"=>"all")
    root=Dict("species"=>order_species,"phases"=>[order_phase],"units"=>Dict())
    base=Dict("equation"=>"H2 (+M) => H + H (+M)","type"=>"chemically-activated",
        "low-P-rate-constant"=>Dict("A"=>1.0),"high-P-rate-constant"=>Dict("A"=>2.0))
    for fields in (Dict{String,Any}(),Dict("SRI"=>[1.0,0.0,0.0,0.0,0.0]),
                   Dict("Tsang"=>[1.0,0.0,0.0]),Dict("type"=>"Chebyshev"))
        r=merge(base,fields)
        @test_throws ArgumentError Arrhenius._native_gas_arrays(merge(root,Dict("reactions"=>[r])),order_phase,["H2","H","AR"])
    end
    bad_burke=merge(Dict("equation"=>"H2 => H + H","type"=>"linear-Burke"),
        Dict("colliders"=>[Dict("name"=>"AR","type"=>"falloff")]))
    @test_throws ArgumentError Arrhenius._native_gas_arrays(merge(root,Dict("reactions"=>[bad_burke])),order_phase,["H2","H","AR"])
end
