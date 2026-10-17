using Test, Arrhenius, NPZ

@testset "counterflow pseudo-arclength continuation" begin
    mechanism=joinpath(dirname(pathof(Arrhenius)),"..","mechanism","h2o2.yaml")
    if isfile(mechanism*".npz")
        gas=CreateSolution(mechanism)
        grid=collect(range(0,.018;length=12))
        f=CounterflowDiffusionFlame(gas;fuel="H2:1",oxidizer="O2:1",
            mdot_fuel=.5,mdot_oxidizer=3.,T_oxidizer=500.,width=.018,grid)
        solve!(f;refine_grid=false,max_time_steps=100)
        base=(f.fuel_mass_flux,f.oxidizer_mass_flux)
        curve=pseudo_arclength!(f;parameter=:inlet_flux_scale,steps=2,ds=.01,tolerance=1e-7)
        @test length(curve)==3
        @test all(p->p.residual_inf<1e-7,curve)
        @test all(diff([p.arclength for p in curve]).>0)
        @test curve[2].parameter>curve[1].parameter
        @test f.fuel_mass_flux≈base[1]*curve[end].parameter rtol=1e-10
        @test f.oxidizer_mass_flux≈base[2]*curve[end].parameter rtol=1e-10
        @test all(size(p.state)==size(f.state) for p in curve)
        @test all(p->p.max_strain_rate>0,curve)
        @test !any(p.turning_point for p in curve)
        checkpoint=joinpath(mktempdir(),"counterflow-continuation.npz")
        save_continuation_checkpoint(checkpoint,f,curve[2];parameter=:inlet_flux_scale)
        f2=CounterflowDiffusionFlame(gas;fuel="H2:1",oxidizer="O2:1",
            mdot_fuel=.5,mdot_oxidizer=3.,T_oxidizer=500.,width=.018,grid)
        resumed=pseudo_arclength!(f2;restart=checkpoint,steps=1,tolerance=1e-7)
        rm(dirname(checkpoint);recursive=true,force=true)
        @test resumed[1].state==curve[2].state
        @test resumed[1].arclength==curve[2].arclength
        @test resumed[2].parameter≈curve[3].parameter rtol=1e-12
        @test resumed[2].state≈curve[3].state rtol=1e-9 atol=1e-10

        controlled=two_point_control!(f;max_steps=1)
        @test length(controlled)==1
        @test controlled[1].max_temperature≈maximum(temperature(f))
        @test controlled[1].max_strain_rate>0
        @test controlled[1].solve_seconds>=0
        @test !isnothing(f.control_points)
        @test length(f.grid)==12
        @test controlled[1].state==f.state
    end
end


@testset "Cantera two-point control first fixed-grid step" begin
    mechanism=joinpath(dirname(pathof(Arrhenius)),"..","mechanism","h2o2.yaml")
    fixture=joinpath(dirname(pathof(Arrhenius)),"..","test","data","cantera_h2o2_two_point_12pt.npz")
    if isfile(mechanism*".npz") && isfile(fixture)
        ref=NPZ.npzread(fixture)
        gas=CreateSolution(mechanism)
        z=vec(Float64.(ref["initial_grid_m"])); Y=Matrix{Float64}(ref["initial_Y"])
        T=vec(Float64.(ref["initial_T_K"])); P=Float64(only(ref["pressure_Pa"]))
        mf=Float64(only(ref["initial_fuel_mdot"])); mo=Float64(only(ref["initial_oxidizer_mdot"]))
        f=CounterflowDiffusionFlame(gas;fuel="H2:1",oxidizer="O2:1",mdot_fuel=mf,
            mdot_oxidizer=mo,T_fuel=300.,T_oxidizer=500.,P,width=z[end]-z[1],grid=z)
        rho=P .* (1 ./ (Y*(1 ./ gas.MW))) ./ (Arrhenius.R .* T)
        f.state[1,:] .= T ./ 1000
        f.state[2:gas.n_species+1,:] .= permutedims(Y)
        f.state[gas.n_species+2,:] .= vec(Float64.(ref["initial_spreadRate"])) ./ Arrhenius._counterflow_Vscale
        f.state[gas.n_species+3,:] .= vec(Float64.(ref["initial_Lambda"])) ./ Arrhenius._counterflow_Lscale
        f.state[end,:] .= rho .* vec(Float64.(ref["initial_velocity_m_s"]))
        f.converged=true
        points=two_point_control!(f;max_steps=1,
            temperature_decrement_schedule=vec(Float64.(ref["curve_temperature_decrement_K"][1:1])))
        @test length(points)==1
        @test length(f.grid)==length(z) && f.grid==z
        @test points[1].max_temperature≈Float64(ref["Tmax_curve_K"][1]) rtol=1e-3
        @test points[1].max_strain_rate≈Float64(ref["strain_rate_curve_per_s"][1]) rtol=1e-3
    end
end


@testset "Cantera two-point S-curves on fixed grids" begin
    cases = (
        ("h2o2.yaml", "cantera_h2o2_two_point_12pt.npz", "H2:1", "O2:1"),
        ("gri30.yaml", "cantera_gri30_nh3h2_two_point_24pt.npz", "NH3:0.85,H2:0.15", "O2:0.21,N2:0.79"),
    )
    for (mechanism_name, fixture_name, fuel, oxidizer) in cases
        mechanism=joinpath(dirname(pathof(Arrhenius)),"..","mechanism",mechanism_name)
        fixture=joinpath(dirname(pathof(Arrhenius)),"..","test","data",fixture_name)
        if isfile(mechanism*".npz") && isfile(fixture)
            ref=NPZ.npzread(fixture)
            gas=CreateSolution(mechanism)
            z=vec(Float64.(ref["initial_grid_m"])); Y=Matrix{Float64}(ref["initial_Y"])
            T=vec(Float64.(ref["initial_T_K"])); P=Float64(only(ref["pressure_Pa"]))
            mf=Float64(only(ref["initial_fuel_mdot"])); mo=Float64(only(ref["initial_oxidizer_mdot"]))
            f=CounterflowDiffusionFlame(gas;fuel,oxidizer,mdot_fuel=mf,mdot_oxidizer=mo,
                T_fuel=Float64(only(get(ref,"T_fuel_K",[300.]))),T_oxidizer=Float64(only(get(ref,"T_oxidizer_K",[500.]))),
                P,width=z[end]-z[1],grid=z)
            rho=P .* (1 ./ (Y*(1 ./ gas.MW))) ./ (Arrhenius.R .* T)
            f.state[1,:] .= T ./ 1000
            f.state[2:gas.n_species+1,:] .= permutedims(Y)
            f.state[gas.n_species+2,:] .= vec(Float64.(ref["initial_spreadRate"])) ./ Arrhenius._counterflow_Vscale
            f.state[gas.n_species+3,:] .= vec(Float64.(ref["initial_Lambda"])) ./ Arrhenius._counterflow_Lscale
            f.state[end,:] .= rho .* vec(Float64.(ref["initial_velocity_m_s"]))
            f.converged=true
            schedule=vec(Float64.(ref["curve_temperature_decrement_K"]))
            curve=two_point_control!(f;max_steps=length(schedule),temperature_decrement_schedule=schedule)
            refT=vec(Float64.(ref["Tmax_curve_K"]))
            refA=vec(Float64.(ref["strain_rate_curve_per_s"]))
            @test length(curve)==length(refT)==length(refA)
            @test f.grid==z
            @test maximum(abs(curve[i].max_temperature-refT[i])/refT[i] for i in eachindex(refT))<1e-3
            @test maximum(abs(curve[i].max_strain_rate-refA[i])/refA[i] for i in eachindex(refA))<1e-3
            left_idx=vec(Int.(ref["curve_left_control_index"])) .+ 1
            right_idx=vec(Int.(ref["curve_right_control_index"])) .+ 1
            left_z=vec(Float64.(ref["curve_left_control_z_m"]))
            right_z=vec(Float64.(ref["curve_right_control_z_m"]))
            left_T=vec(Float64.(ref["curve_left_control_T_K"]))
            right_T=vec(Float64.(ref["curve_right_control_T_K"]))
            @test length(left_idx)==length(curve)==length(right_idx)
            @test all(curve[i].left_control_index==left_idx[i] for i in eachindex(curve))
            @test all(curve[i].right_control_index==right_idx[i] for i in eachindex(curve))
            max_dz=0.
            max_dT=0.
            differing_steps=Int[]
            for i in eachindex(curve)
                native=curve[i]
                nzL,nTL,nzR,nTR=native.control_points
                max_dz=max(max_dz,abs(nzL-left_z[i]),abs(nzR-right_z[i]))
                max_dT=max(max_dT,abs(nTL-left_T[i]),abs(nTR-right_T[i]))
                (native.left_control_index==left_idx[i] && native.right_control_index==right_idx[i]) || push!(differing_steps,i)
                @test isapprox(nTL,left_T[i];rtol=1e-3,atol=1e-5)
                @test isapprox(nTR,right_T[i];rtol=1e-3,atol=1e-5)
                println("CONTROL_POINT_STEP ",(mechanism=mechanism_name,step=i,
                    cantera_left=(index=left_idx[i],z=left_z[i],T=left_T[i]),
                    arrhenius_left=(index=native.left_control_index,z=nzL,T=nTL),
                    cantera_right=(index=right_idx[i],z=right_z[i],T=right_T[i]),
                    arrhenius_right=(index=native.right_control_index,z=nzR,T=nTR)))
            end
            println("CONTROL_POINT_SUMMARY ",(mechanism=mechanism_name,steps=length(curve),max_abs_dz=max_dz,max_abs_dT=max_dT,index_differing_steps=differing_steps))
            native_peak=argmax(p.max_strain_rate for p in curve)
            cantera_peak=argmax(refA)
            @test abs(curve[native_peak].max_strain_rate-refA[cantera_peak])/refA[cantera_peak]<1e-3
            @test abs(curve[native_peak].max_temperature-refT[cantera_peak])/refT[cantera_peak]<1e-3
            if mechanism_name=="h2o2.yaml"
                fold_state=copy(f.state); fold_control=f.control_points
                fold_mf=f.fuel_mass_flux; fold_mo=f.oxidizer_mass_flux
                arc_t0=time_ns()
                arc=pseudo_arclength!(f;parameter=:control_temperature,steps=40,ds=.003,
                    min_ds=1e-6,max_ds=.003,tolerance=1e-7,max_corrector_iters=36,state_scale=1000.)
                arc_seconds=(time_ns()-arc_t0)/1e9
                @test length(arc)==41
                @test abs(arc[1].max_strain_rate-refA[cantera_peak])/refA[cantera_peak]<1e-3
                @test abs(arc[1].max_temperature-refT[cantera_peak])/refT[cantera_peak]<1e-3
                @test all(p->p.residual_inf<1e-7,arc)
                @test isapprox(maximum(p.max_strain_rate for p in arc),refA[cantera_peak];rtol=1e-3)
                println("CONTROLLED_ARC_N12 ", (steps=length(arc)-1, points=length(arc), elapsed_s=arc_seconds, peak_strain=maximum(p.max_strain_rate for p in arc), peak_Tmax=arc[argmax(p.max_strain_rate for p in arc)].max_temperature))
                plus=arc; plus_seconds=arc_seconds
                f.state .= fold_state; f.control_points=fold_control
                f.fuel_mass_flux=fold_mf; f.oxidizer_mass_flux=fold_mo; f.converged=true
                minus_t0=time_ns()
                minus=pseudo_arclength!(f;parameter=:control_temperature,steps=1,ds=.0001,
                    min_ds=1e-6,max_ds=.0001,tolerance=1e-7,max_corrector_iters=36,state_scale=1000.,direction=-1)
                minus_seconds=(time_ns()-minus_t0)/1e9
                @test length(plus)==41 && length(minus)==2
                @test all(p->p.residual_inf<1e-7,plus) && all(p->p.residual_inf<1e-7,minus)
                @test isapprox(maximum(p.max_strain_rate for p in plus),refA[cantera_peak];rtol=1e-3)
                @test abs(minus[1].max_strain_rate-refA[cantera_peak])/refA[cantera_peak]<1e-3
                @test abs(minus[1].max_temperature-refT[cantera_peak])/refT[cantera_peak]<1e-3
                println("CONTROLLED_ARC_N12 direction=+1 ", (steps=length(plus)-1, seconds=plus_seconds, end_control_T=plus[end].parameter, peak_strain=maximum(p.max_strain_rate for p in plus)))
                println("CONTROLLED_ARC_N12 direction=-1 ", (steps=length(minus)-1, seconds=minus_seconds, end_control_T=minus[end].parameter, peak_strain=maximum(p.max_strain_rate for p in minus)))
            end
        end
    end
end
