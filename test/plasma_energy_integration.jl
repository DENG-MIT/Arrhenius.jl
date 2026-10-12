using Arrhenius, Test, LinearAlgebra
using Arrhenius: YAML
include(joinpath(dirname(pathof(Arrhenius)), "..", "example", "reactors", "ode_solver.jl"))

@testset "Energy reactor QNDF and explicit cache restart" begin
    mktempdir() do dir
        names=["H2","H","H+","e"]
        specs=Any[]
        for (name,composition,cp) in [
                ("H2",Dict("H"=>2),3.5),("H",Dict("H"=>1),2.5),
                ("H+",Dict("H"=>1,"E"=>-1),2.5),("e",Dict("E"=>1),2.5)]
            push!(specs,Dict("name"=>name,"composition"=>composition,
                "thermo"=>Dict("model"=>"NASA7","temperature-ranges"=>[200.,6000.],
                    "data"=>[[cp,0.,0.,0.,0.,-10000.,1.]])))
        end
        phase=Dict("name"=>"inert-plasma","thermo"=>"plasma","elements"=>["H","E"],
            "species"=>names,"reactions"=>[Dict("collisions"=>"all")],
            "electron-energy-distribution"=>Dict("type"=>"Boltzmann-two-term",
                "energy-levels"=>collect(0.:.25:10.)),
            "state"=>Dict("T"=>300.,"P"=>101325.,
                "X"=>Dict("H2"=>.998,"H"=>0.,"H+"=>.001,"e"=>.001)))
        root=Dict("units"=>Dict("length"=>"m","quantity"=>"kmol","activation-energy"=>"J/kmol"),
            "phases"=>[phase],"species"=>specs,
            "collisions"=>[Dict("equation"=>"H2 + e => H2 + e",
                "type"=>"electron-collision-plasma",
                "energy-levels"=>[0.,10.],"cross-sections"=>[1e-20,1e-20])])
        path=joinpath(dir,"inert.yaml"); YAML.write_file(path,root)
        m=PlasmaMechanism(path); s=PlasmaState(m)
        update_eedf!(s); set_reduced_electric_field!(s,2e-21)
        r=PlasmaEnergyReactor(s;volume=.25); problem=reactor_problem(r,(0.,1.))
        rhs=problem.f; u0=copy(problem.u0); u0[4]=-1e-40
        p0=reactor_properties(rhs,u0)
        @test u0[2]<0
        @test p0.Y[2]<0
        @test rhs.electric_field>0 && rhs.mobility>0
        n_jac=Ref(0)
        jac=(J,u,p,t)->begin n_jac[]+=1; reactor_jacobian!(J,u,rhs,t); nothing end
        initial=copy(u0)
        roundoff_steps=0
        for chunk in 1:2
            p=reactor_properties(rhs,u0)
            field=rhs.electric_field; mobility=rhs.mobility; distribution=copy(rhs.eedf.center_eedf)
            # Exact inert-mixture solution: fixed Y and total-H derivative.
            # Constants are SI electron charge and Avogadro number per kmol.
            rate=u0[1]*u0[m.electron_index+2]*(1.602176634e-19*6.02214076e26/m.MW[m.electron_index])*mobility*field^2
            capacity=u0[1]*Arrhenius.R*sum((k==m.electron_index ? 0. : (k==1 ? 3.5 : 2.5))*u0[k+2]/m.MW[k] for k in 1:m.n_species)
            dt=5capacity/rate
            times=collect(range(0.,dt;length=5))
            chunk_problem=merge(problem,(u0=copy(u0),tspan=(0.,dt),jac=jac))
            sol=native_bdf(chunk_problem;reltol=1e-12,abstol=1e-19,saveat=times,save_everystep=false)
            roundoff_steps+=sol.destats.naccept+sol.destats.nreject
            # Accumulated arithmetic error grows with integrator updates;
            # the physical conservation ceiling remains independent of steps.
            invariant_rtol=min(1e-10,64eps(Float64)*(sol.destats.naccept+sol.destats.nreject+1))
            @test SciMLBase.successful_retcode(sol)
            @test sol.t==times
            @test n_jac[]>0
            @test maximum(abs(sol.u[i][2]-(u0[2]+rate*times[i])) for i in eachindex(times)) <= 1e-8*abs(u0[2])
            @test all(isapprox(v[1],u0[1];rtol=invariant_rtol,atol=0.) && all(isapprox.(v[3:end],u0[3:end];rtol=invariant_rtol,atol=0.)) for v in sol.u)
            @test all(v[4]<0 && v[2]<0 for v in sol.u)
            last=reactor_properties(rhs,sol.u[end])
            @info "energy integration accuracy" chunk H_error=sol.u[end][2]-(u0[2]+rate*dt) T_error=last.T-(p.T+5) predicted_T_error=(sol.u[end][2]-(u0[2]+rate*dt))/capacity enthalpy_offset_temperature=abs(u0[2])/capacity
            @test last.T≈p.T+5 rtol=1e-9
            @test last.Te==p0.Te
            @test rhs.electric_field==field && rhs.mobility==mobility && rhs.eedf.center_eedf==distribution
            u0=copy(sol.u[end])
            if chunk==1
                update_eedf!(rhs,u0;reduced_field=4e-21)
                @test rhs.electric_field>field
            end
        end
        @test problem.u0==reactor_state(r)
        @test u0[2]>initial[2] && all(isapprox.(u0[3:end],initial[3:end];rtol=min(1e-10,64eps(Float64)*(roundoff_steps+1)),atol=0.))
    end
end
