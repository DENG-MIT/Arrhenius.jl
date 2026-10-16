module IonTransportTests
using Arrhenius, Test, NPZ, YAML, SHA, LinearAlgebra

function fixture(dir)
    path=joinpath(dir,"ion.yaml")
    names=["Ar","He","Ar+","E"]
    comps=[Dict("Ar"=>1),Dict("He"=>1),Dict("Ar"=>1,"E"=>-1),Dict("E"=>1)]
    specs=[Dict("name"=>name,"composition"=>comp,"thermo"=>Dict("model"=>"NASA7",
        "temperature-ranges"=>[200.,6000.],"data"=>[[2.5,0.,0.,0.,0.,0.,0.]])) for (name,comp) in zip(names,comps)]
    YAML.write_file(path,Dict("phases"=>[Dict("name"=>"gas","thermo"=>"ideal-gas",
        "transport"=>"ionized-gas","elements"=>["Ar","He","E"],"species"=>names)],"species"=>specs))
    mu=zeros(5,4);mu[1,:]=[1e-3,2e-3,1.5e-3,1e-4]
    lam=zeros(5,4);lam[1,:]=[.01,.02,.03,.04]
    binary=zeros(5,16);binary[1,:].=.001
    a=Dict{String,Any}("molecular_weights"=>[39.95,4.002602,39.95-Arrhenius._ION_ELECTRON_MW,Arrhenius._ION_ELECTRON_MW],
        "species_viscosities_poly"=>mu,"thermal_conductivity_poly"=>lam,"binary_diff_coeffs_poly"=>binary,
        "transport_model_utf8"=>collect(codeunits("ionized-gas")),
        "source_sha256_utf8"=>collect(codeunits(bytes2hex(sha256(read(path))))),
        "sidecar_format_utf8"=>collect(codeunits("arrhenius-sidecar-v2")))
    npzwrite(path*".npz",a)
    CreateSolution(path),path,a
end

function as_ionized(g)
    t=g.trans
    trans=Arrhenius.Transport(t.poly_order,t.species_viscosities_poly,t.thermal_conductivity_poly,t.binary_diff_coeffs_poly,:ionized_gas)
    Arrhenius.Solution(g.n_species,g.n_reactions,g.MW,g.species_names,g.elements,g.ele_matrix,g.thermo,trans,g.reaction)
end

@testset "ionized transport minimum contracts" begin
    mktempdir() do dir
        g,path,a=fixture(dir);d=IonTransportData(g);w=IonTransportWorkspace(d)
        @test d.electron==4 && d.ions==[3] && d.neutrals==[1,2]
        @test d.charges==[0.,0.,1.,-1.]
        @test convert_precision(g,Float32).trans.model==:ionized_gas
        @test IonTransportData(convert_precision(g,Float32)).electron==4
        @test_throws ArgumentError mixture_transport!(TransportWorkspace(g),g,one_atm,1000.,fill(.25,4))
        @test FreeFlame(g;T=900.,X=Dict("Ar"=>1.0)) isa IonizedFlame
        bad=copy(a);delete!(bad,"transport_model_utf8");npzwrite(path*".npz",bad)
        @test_throws ArgumentError CreateSolution(path)
        npzwrite(path*".npz",a)

        X=[.7,.28,.01,.01];T=1000.;P=one_atm
        @test_throws DomainError ionized_flux!(zeros(4),w,d,zeros(4),X,X;density=1.)
        mu,lam,sigma=ionized_transport!(w,d,P,T,X)
        @test isfinite(mu) && mu>0 && lam>0 && sigma>0
        @test w.diffusion[4]==.4*Arrhenius._ION_KB*T/Arrhenius._ION_QE
        @test w.mobility[4]==.4 && w.mobility[1:2]==[0.,0.]
        Dp=.001*T*sqrt(T)
        @test w.mobility[3] ≈ Dp/P*Arrhenius._ION_QE/(Arrhenius._ION_KB*T)/.98 rtol=1e-14
        @test sigma ≈ P/(Arrhenius._ION_KB*T)*Arrhenius._ION_QE*(.01*w.mobility[3]+.01*.4) rtol=1e-14
        savedD=copy(w.diffusion);savedM=copy(w.mobility)
        ionized_transport!(w,d,2P,T,X)
        @test w.diffusion[1:3] ≈ savedD[1:3]/2 rtol=1e-14
        @test w.diffusion[4]==savedD[4]
        @test w.mobility[3] ≈ savedM[3]/2 rtol=1e-14
        @test w.mobility[4]==savedM[4]
        for (p,t,x) in ((-P,T,X),(P,0.,X),(P,T,[Inf,0.,0.,0.]),(P,T,zeros(4)))
            @test_throws ArgumentError ionized_transport!(w,d,p,t,x)
        end
        @test_throws DimensionMismatch ionized_transport!(w,d,P,T,[1.])
        @test_throws ArgumentError ionized_transport!(w,d,P,T,X;mean_molecular_weight=0.)
        signed=[.7,.3,-1e-12,1e-12]
        ionized_transport!(w,d,P,T,signed;mean_molecular_weight=20.)
        @test w.X[3]==1e-20 && w.mean_MW==20.
        @test signed[3]==-1e-12
        ionized_transport!(w,d,P,T,X)

        Y=[.6,.38,.019,.001];Yr=[.62,.36,.018,.002];grad=[1.,-2.,.5,.5];rho=.7;E=200.
        pref=d.molecular_weights.*(P/(Arrhenius.R*T)).*w.diffusion
        raw=-pref.*grad
        drift=rho.*((Y.+Yr)./2).*E.*d.charges.*w.mobility
        expected=raw.+drift
        expected[1:2].-=sum(expected).*Y[1:2]/sum(Y[1:2])
        flux=zeros(4)
        @test ionized_flux!(flux,w,d,grad,Y,Yr;density=rho,electric_field=E) === flux
        @test flux ≈ expected rtol=1e-14
        @test abs(sum(flux)) < 1e-13*sum(abs,flux)
        frozen=[raw[1],raw[2],0.,0.]
        frozen[1:2].-=sum(frozen).*Y[1:2]
        ionized_flux!(flux,w,d,grad,Y,Yr;density=rho,electric_field=E,frozen=true)
        @test flux ≈ frozen rtol=1e-14
        @test flux[3:4]==[0.,0.]
        ionized_flux!(flux,w,d,zeros(4),Y,Yr;density=rho,electric_field=E)
        @test flux[3]>0 && flux[4]<0
        plus=copy(flux)
        ionized_flux!(flux,w,d,zeros(4),Y,Yr;density=rho,electric_field=-E)
        @test flux ≈ -plus rtol=1e-14
        for args in ((zeros(3),grad,Y,Yr),(flux,zeros(3),Y,Yr),(flux,grad,zeros(3),Yr))
            @test_throws DimensionMismatch ionized_flux!(args[1],w,d,args[2:4]...;density=rho)
        end
        @test_throws ArgumentError ionized_flux!(flux,w,d,grad,Y,Yr;density=-rho)
        @test_throws ArgumentError ionized_flux!(flux,w,d,grad,Y,Yr;density=rho,electric_field=Inf)
        trial=[.7,.5,-1e-10,2e-10]
        ionized_flux!(flux,w,d,zeros(4),trial,trial;density=rho,electric_field=E)
        expected=rho.*trial.*E.*d.charges.*w.mobility
        expected[1:2].-=sum(expected).*trial[1:2]/(1-sum(trial[3:4]))
        @test flux ≈ expected rtol=1e-14
        @test trial[3]==-1e-10 && sum(trial)>1

        neutral=CreateSolution(joinpath(@__DIR__,"..","mechanism","h2o2.yaml"))
        @test_throws ArgumentError IonTransportData(neutral)
        nd=IonTransportData(as_ionized(neutral));nw=IonTransportWorkspace(nd)
        nx=collect(1.:neutral.n_species);nx./=sum(nx)
        nv,nk=mixture_transport!(TransportWorkspace(neutral),neutral,P,T,nx)
        iv,ik,ic=ionized_transport!(nw,nd,P,T,nx)
        @test iv ≈ nv rtol=1e-13
        @test ik ≈ nk rtol=1e-13
        @test ic==0. && all(iszero,nw.mobility)
    end
end
end
