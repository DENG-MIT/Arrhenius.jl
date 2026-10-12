module SignedEquilibriumTests
using Arrhenius, Test, NPZ, YAML, SHA, LinearAlgebra

function ion_fixture(dir, delta)
    path=joinpath(dir,"saha.yaml")
    names=["Ar","Ar+","E"]
    comps=[Dict("Ar"=>1),Dict("Ar"=>1,"E"=>-1),Dict("E"=>1)]
    entropy0=2.5-2.5log(1000.)
    specs=[Dict("name"=>names[k],"composition"=>comps[k],"thermo"=>Dict(
        "model"=>"NASA7","temperature-ranges"=>[200.,6000.],
        "data"=>[[2.5,0.,0.,0.,0.,0.,entropy0-(k==3 ? delta : 0.)]])) for k in 1:3]
    YAML.write_file(path,Dict("phases"=>[Dict("name"=>"gas","thermo"=>"ideal-gas",
        "elements"=>["Ar","E"],"species"=>names)],"species"=>specs))
    npzwrite(path*".npz",Dict("molecular_weights"=>[39.95,39.95-Arrhenius._ION_ELECTRON_MW,Arrhenius._ION_ELECTRON_MW],
        "source_sha256_utf8"=>collect(codeunits(bytes2hex(sha256(read(path))))),
        "sidecar_format_utf8"=>collect(codeunits("arrhenius-sidecar-v2"))))
    CreateSolution(path)
end
function rows(g,extra)
    Arrhenius.Solution(g.n_species,g.n_reactions,g.MW,g.species_names,
        [g.elements;"extra"],vcat(g.ele_matrix,extra),g.thermo,g.trans,g.reaction)
end

@testset "signed equilibrium minimal analytic and derivative contracts" begin
    mktempdir() do dir
        # Saha mass action: X_ion*X_e/X_neutral = exp(-delta)/(P/P0).
        # At delta=800 the ion populations are representable although K underflows.
        for delta in (0.,100.,800.), pressure in (one_atm,5one_atm)
            g=ion_fixture(dir,delta)
            eq=equilibrate(g;T=1000.,P=pressure,X="Ar")
            logK=-delta-log(pressure/one_atm)
            rootK=exp(logK/2)
            xe=rootK/(rootK+sqrt(1+exp(logK)))
            @test all(isfinite,eq.X) && all(>=(0),eq.X)
            @test sum(eq.X) ≈ 1 atol=1e-14
            @test eq.X[1] ≈ 1-2xe atol=1e-13
            @test eq.X[2]>0 && eq.X[3]>0
            @test log(eq.X[2]/xe) ≈ 0 atol=5e-10
            @test log(eq.X[3]/xe) ≈ 0 atol=5e-10
        end
        g=ion_fixture(dir,0.)
        for initial in ([1.,0.,0.],[.8,.15,.05],[.8,.05,.15])
            before=copy(initial)
            eq=equilibrate(g;T=1000.,P=one_atm,X=initial)
            @test initial==before
            @test g.ele_matrix*eq.X/dot(g.MW,eq.X) ≈ g.ele_matrix*initial/dot(g.MW,initial) atol=1e-13 rtol=1e-10
            @test eq.X[2]*eq.X[3]/eq.X[1] ≈ 1 rtol=1e-9

            for volume in (false,true)
                sys=Arrhenius._equilibrium_system(g,initial)
                sys.g .= [0.,0.,1600.]
                state=[.4,750.,.2] # one or both charged mole fractions underflow
                # QR may permute the original element rows; put the charge potential
                # on the signed row rather than assuming a particular permutation.
                charge=findfirst(i->any(<(0),@view(sys.A[i,:])),axes(sys.A,1))
                fill!(state,0.2);state[charge]=750.
                f,x,J=Arrhenius._equilibrium_evaluate(sys,state,.3,volume;jacobian=true)
                analytic=copy(J);res=copy(f)
                @test all(isfinite,res) && all(isfinite,analytic)
                finite=zeros(size(J))
                for j in eachindex(state)
                    plus=copy(state);minus=copy(state);h=1e-4
                    plus[j]+=h;minus[j]-=h
                    fp=copy(first(Arrhenius._equilibrium_evaluate(sys,plus,.3,volume)))
                    fm=copy(first(Arrhenius._equilibrium_evaluate(sys,minus,.3,volume)))
                    finite[:,j]=(fp-fm)/(2h)
                end
                @test analytic ≈ finite rtol=2e-7 atol=2e-8
            end
        end
        baseline=equilibrate(g;T=1000.,P=one_atm,X="Ar")
        duplicate=rows(g,2g.ele_matrix[2:2,:])
        @test equilibrate(duplicate;T=1000.,X="Ar").X ≈ baseline.X rtol=1e-10
        # An absent ordinary atom excludes Ar+, which leaves only one charge sign;
        # the remaining electron must then also be excluded to preserve zero charge.
        filtered=rows(g,[0. 1. 0.])
        @test Arrhenius._equilibrium_system(filtered,[1.,0.,0.]).species==[1]
        @test equilibrate(filtered;T=1000.,X="Ar").X==[1.,0.,0.]
        @test equilibrate(g;T=1000.,X="E").X==[0.,0.,1.]
    end
end
end
