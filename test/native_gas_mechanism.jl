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

@testset "native gas kinetics parity edge cases" begin
    species = [
        Dict("name"=>"H2","composition"=>Dict("H"=>2)),
        Dict("name"=>"H","composition"=>Dict("H"=>1)),
        Dict("name"=>"AR","composition"=>Dict("Ar"=>1)),
    ]
    names = ["H2","H","AR"]
    phase = Dict("name"=>"gas","elements"=>["H","Ar"],"species"=>names,"reactions"=>"all")
    root_units = Dict{String,Any}()
    function edge_arrays(reaction; units=root_units, section_units=nothing, validate=true)
        reactions = section_units === nothing ? Any[reaction] :
            Any[Dict("units"=>section_units), reaction]
        root = Dict("species"=>species,"elements"=>["H","Ar"],"phases"=>[phase],
                    "units"=>units,"reactions"=>reactions)
        Arrhenius._native_gas_arrays(root,phase,names;validate)
    end
    function edge_reaction(equation, fields=Dict{String,Any}())
        merge(Dict{String,Any}("equation"=>equation,
            "rate-constant"=>Dict("A"=>1.0,"b"=>0.0,"Ea"=>0.0)),
            Dict{String,Any}(fields))
    end

    @testset "activation-energy K units" begin
        units=Dict("activation-energy"=>"K")
        parsed=Arrhenius._gas_arrhenius(
            Dict("A"=>1.0e10,"b"=>0.0,"Ea"=>10000.0),1,
            Arrhenius._gas_units(Dict("units"=>units)))
        kf=parsed[1]*exp(-parsed[3]*(4184.0/Arrhenius.R)/1000.0)
        @test isapprox(kf,1.0e10*exp(-10.0);rtol=2e-15)
        @test Arrhenius._gas_quantity("10000 K",1.0,"Ea")==10000.0*Arrhenius.R
    end

    @testset "Troe T2 zero is omitted" begin
        fixture=joinpath(@__DIR__,"data","native_gas_features.yaml")
        yaml=replace(read(fixture,String),
            "Troe: {A: 0.5, T3: 1000, T1: 5000}"=>
            "Troe: {A: 0.5, T3: 100, T1: 1000, T2: 0}")
        mktempdir() do dir
            path=joinpath(dir,"troe_zero.yaml")
            write(path,yaml)
            gas=CreateSolution(path)
            root=YAML.load_file(path)
            reaction=root["reactions"][4]
            units=Arrhenius._gas_units(root,get(reaction,"units",Dict()))
            high=Arrhenius._gas_arrhenius(reaction["high-P-rate-constant"],2,units;name="high-P-rate-constant")
            low=Arrhenius._gas_arrhenius(reaction["low-P-rate-constant"],3,units;name="low-P-rate-constant")
            T=1000.0; P=one_atm
            scale=4184.0/Arrhenius.R/T
            kinf=high[1]*exp(high[2]*log(T)-high[3]*scale)
            k0=low[1]*exp(low[2]*log(T)-low[3]*scale)
            collider=P/(Arrhenius.R*T)
            pr=k0*collider/kinf
            fcent=0.5*exp(-T/100.0)+0.5*exp(-T/1000.0)
            lfc=log10(fcent); lpr=log10(pr)
            c=-0.4-0.67*lfc; n=0.75-1.27*lfc
            f1=(lpr+c)/(n-0.14*(lpr+c))
            expected=kinf*(pr/(1+pr))*10.0^(lfc/(1+f1^2))
            X=fill(1.0/gas.n_species,gas.n_species)
            actual=reaction_rate_constants(gas;T,P,X).forward[4]
            @test isapprox(actual,expected;rtol=2e-13)
            @test fcent < 1.0
        end
    end

    @testset "PLOG pseudo-colliders and pressure units" begin
        units=Dict("length"=>"cm","quantity"=>"mol","pressure"=>"Pa")
        rates=[Dict("P"=>1.0,"A"=>1.0e12,"b"=>0.0,"Ea"=>0.0),
               Dict("P"=>10.0,"A"=>2.0e12,"b"=>0.0,"Ea"=>0.0)]
        pseudo=edge_reaction("H2 (+H) <=> H + H (+H)",
            Dict("type"=>"pressure-dependent-Arrhenius","rate-constants"=>rates))
        parsed=edge_arrays(pseudo;units)
        @test parsed["Plog_collider_indices"]==[0]
        @test parsed["Plog_Arrhenius"][1,1]==1.0e12
        plog_m=edge_reaction("H2 (+M) <=> 2 H (+M)",
            Dict("type"=>"pressure-dependent-Arrhenius","rate-constants"=>rates))
        @test_throws ArgumentError edge_arrays(plog_m;units)
        @test_logs (:warn,r"superfluous M collider.*H2 \(\+M\)") begin
            @test edge_arrays(plog_m;units,validate=false)["Plog_collider_indices"]==[0]
        end
        section=edge_reaction("H2 <=> 2 H",
            Dict("type"=>"pressure-dependent-Arrhenius","rate-constants"=>rates))
        pressure= edge_arrays(section;section_units=Dict("pressure"=>"atm"))
        @test all(isapprox.(pressure["Plog_pressures"],[101325.0,1013250.0];rtol=2e-15))
        local_pressure=edge_arrays(merge(section,Dict("units"=>Dict("pressure"=>"atm"))))
        @test all(isapprox.(local_pressure["Plog_pressures"],[101325.0,1013250.0];rtol=2e-15))
    end

    @testset "load-time rate validation" begin
        negative=edge_reaction("H2 => 2 H",Dict("rate-constant"=>Dict("A"=>-1.0)))
        @test_throws ArgumentError edge_arrays(negative)
        @test_logs (:warn, r"negative A without negative-A") edge_arrays(negative;validate=false)
        flagged=edge_reaction("H2 => 2 H",
            Dict("negative-A"=>true,"rate-constant"=>Dict("A"=>-1.0)))
        @test edge_arrays(flagged)["Arrhenius_coeffs"][1,1]==-1.0

        fixture=joinpath(@__DIR__,"data","native_gas_features.yaml")
        yaml=replace(read(fixture,String),"A: 1.0e10, b: 0.5"=>"A: -1.0e10, b: 0.5")
        mktempdir() do dir
            path=joinpath(dir,"negative_a.yaml")
            write(path,yaml)
            err=try
                CreateSolution(path)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin("H + O <=> OH",sprint(showerror,err))
            @test_logs (:warn, r"H \+ O <=> OH.*negative A") CreateSolution(path;validate=false)
        end

        opposite=edge_reaction("H2 (+M) => 2 H (+M)",
            Dict("type"=>"falloff","negative-A"=>true,
                 "low-P-rate-constant"=>Dict("A"=>1.0),
                 "high-P-rate-constant"=>Dict("A"=>-1.0)))
        @test_throws ArgumentError edge_arrays(opposite)
        @test_logs (:warn,r"same sign.*H2 \(\+M\)") edge_arrays(opposite;validate=false)

        plog_negative=edge_reaction("H2 => 2 H",
            Dict("type"=>"pressure-dependent-Arrhenius",
                 "rate-constants"=>[Dict("P"=>1.0,"A"=>-1.0,"b"=>0.0,"Ea"=>0.0),
                                    Dict("P"=>1.0,"A"=>2.0,"b"=>0.0,"Ea"=>0.0)]))
        @test edge_arrays(plog_negative)["Plog_Arrhenius"][1,1]==-1.0
        nonpositive=edge_reaction("H2 => 2 H",
            Dict("type"=>"pressure-dependent-Arrhenius",
                 "rate-constants"=>[Dict("P"=>1.0,"A"=>-1.0,"b"=>0.0,"Ea"=>0.0)]))
        @test_throws ArgumentError edge_arrays(nonpositive)
        @test_logs (:warn,r"PLOG rate sum must be positive.*H2 => 2 H") edge_arrays(nonpositive;validate=false)

        unknown=edge_reaction("H2 => 2 H",
            Dict("type"=>"unrecognized-rate-type","rate-constant"=>Dict("A"=>1.0)))
        @test_throws ArgumentError edge_arrays(unknown)
        err=try
            edge_arrays(unknown)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("unrecognized-rate-type", sprint(showerror,err))
        @test occursin("H2 => 2 H", sprint(showerror,err))
        # Validation opt-out covers Cantera input checks, not unsupported formulas.
        @test_throws ArgumentError edge_arrays(unknown;validate=false)
    end

    @testset "PLOG logarithmic groups" begin
        @test Arrhenius._positive_floor(0.0)==1.0e-300
        @test Arrhenius._positive_floor(Float32(0))>0

        ea=1000.0*Arrhenius.R*300.0/4184.0
        single=Arrhenius.PlogData([1],[1,3],[1.0e5,1.0e6],[1,2,3],
            [1.0 0.0 ea; 1.0 0.0 0.0])
        rate=Arrhenius._plog_rate(single,1,300.0,sqrt(1.0e11),log(300.0))
        @test isapprox(rate,exp(-500.0);rtol=2e-12)

        multiple=Arrhenius.PlogData([1],[1,3],[1.0e5,1.0e6],[1,3,4],
            [1.0 0.0 ea; 1.0 0.0 ea; 1.0 0.0 0.0])
        grouped=Arrhenius._plog_rate(multiple,1,300.0,sqrt(1.0e11),log(300.0))
        @test isapprox(grouped,1.0e-150;rtol=2e-12)
    end
end
