using Arrhenius, Test, YAML, NPZ, SHA

@testset "selected gas species and sidecar dependencies" begin
    mktempdir() do dir
        path = joinpath(dir,"selected.yaml")
        imported = joinpath(dir,"shared.yaml")
        definitions = [Dict("name"=>name,"composition"=>Dict("Ar"=>1),
            "thermo"=>Dict("model"=>"constant-cp","cp0"=>cp*Arrhenius.R))
            for (name,cp) in (("A",2.5),("B",3.5))]
        YAML.write_file(imported,Dict("species"=>definitions))
        # Receiving-file units must not reinterpret the imported SI thermo.
        root = Dict("units"=>Dict("energy"=>"kJ","quantity"=>"mol","pressure"=>"bar"),
            "phases"=>[Dict("name"=>"selected","thermo"=>"ideal-gas","elements"=>["Ar"],
                "species"=>[Dict("shared.yaml/species"=>["B","A"])],
                "reactions"=>[Dict("shared.yaml/reactions"=>"declared-species")])],
            "reactions"=>[Dict("type"=>"pressure-dependent-Arrhenius"),Dict("type"=>"Blowers-Masel")])
        YAML.write_file(path,root)
        selection = Dict("species_names"=>["B","A"],"n_reactions"=>1,
            "dependencies"=>Dict("shared.yaml"=>bytes2hex(sha256(read(imported)))))
        encoded(x) = collect(codeunits(x))
        payload = Dict{String,Any}("molecular_weights"=>[39.95,39.95],
            "sidecar_format_utf8"=>encoded("arrhenius-sidecar-v4"),
            "source_sha256_utf8"=>encoded(bytes2hex(sha256(read(path)))),
            "phase_selection_utf8"=>encoded(YAML.write(selection)),
            "reactant_stoich_coeffs"=>reshape([1.,0.],2,1),
            "product_stoich_coeffs"=>reshape([0.,1.],2,1),
            "reactant_orders"=>reshape([1.,0.],2,1),
            "efficiencies_coeffs"=>zeros(2,1),"is_reversible"=>[false],
            "Arrhenius_coeffs"=>reshape([1.,0.,0.],1,3))
        save(data=payload) = npzwrite(path*".npz",data)
        save()
        gas = CreateSolution(path)
        @test gas.species_names == ["B","A"]
        @test gas.n_reactions == 1
        @test gas.ele_matrix == ones(1,2)
        cp = zeros(2)
        cal_cp_R!(cp,gas,1000.,one_atm,[.5,.5])
        @test cp ≈ [3.5,2.5] rtol=1e-14
        @test set_states(gas,1000.,one_atm,[.5,.5]) ≈ [-1.,1.] .* (one_atm/Arrhenius.R/1000/2)

        bad=copy(payload); delete!(bad,"phase_selection_utf8"); save(bad)
        @test_throws ArgumentError CreateSolution(path)
        bad=copy(payload); bad["is_reversible"]=[false,false]; save(bad)
        @test_throws ArgumentError CreateSolution(path)
        bad=copy(payload); bad["index_three_body"]=[2]; save(bad)
        @test_throws ArgumentError CreateSolution(path)
        for change in (Dict("species_names"=>["A","B"]),Dict("n_reactions"=>2),
                Dict("dependencies"=>Dict()),Dict("n_reactions"=>true))
            bad=copy(payload)
            bad["phase_selection_utf8"]=encoded(YAML.write(merge(selection,change))); save(bad)
            @test_throws ArgumentError CreateSolution(path)
        end
        bad=copy(payload); bad["sidecar_format_utf8"]=encoded("arrhenius-sidecar-v2"); save(bad)
        @test_throws ArgumentError CreateSolution(path)
        save()
        original = read(imported,String)
        write(imported,replace(original,"\n"=>"\r\n"))
        @test CreateSolution(path).species_names == ["B","A"]
        write(imported,original*"\n# changed dependency\n")
        @test_throws ArgumentError CreateSolution(path)
        write(imported,original)
        searchdir=joinpath(dir,"data"); mkdir(searchdir)
        mv(imported,joinpath(searchdir,"shared.yaml"))
        @test_throws ArgumentError CreateSolution(path)
        @test CreateSolution(path;data_paths=[searchdir]).n_species == 2

        # No reaction arrays are emitted for a selected inert phase.
        inert=copy(payload)
        for key in ("reactant_stoich_coeffs","product_stoich_coeffs","reactant_orders",
                "efficiencies_coeffs","is_reversible","Arrhenius_coeffs")
            delete!(inert,key)
        end
        inert["phase_selection_utf8"]=encoded(YAML.write(merge(selection,Dict("n_reactions"=>0))))
        save(inert)
        gas=CreateSolution(path;data_paths=[searchdir])
        @test gas.n_reactions == 0
        @test set_states(gas,1000.,one_atm,[.5,.5]) == [0.,0.]
    end
end
