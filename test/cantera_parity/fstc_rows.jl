# Emulation of Cantera example flame_solver_timestep_controls.py counterflow rows
# (H2/O2 counterflow diffusion flame, 20-pt grid over 30 mm, h2o2.yaml,
# finite-difference Jacobian, refine ratio=2 slope=0.06 curve=0.08 prune=0.02).
# Runs every distinct case of the example in ONE Julia session and prints
# ok / pts / Tmax / steps / jacobians per row, mirroring Cantera's table.
# usage: julia --project=<worktree> fstc_rows.jl <mech.yaml> <tag> [--dts]
# With --dts, every accepted transient dt is recorded per row and written to
# <tag>_dts.csv for step-by-step comparison with Cantera's time-step callback.
using Arrhenius, LinearAlgebra, NPZ

# label, pressure [Pa], time_step_regrid, max_time_step_count, strategy, growth
function physical_residual!(r,f::CounterflowDiffusionFlame,x,w;previous=nothing,dt=Inf)
    counterflow_residual!(r,f,x,w;previous,dt)
    # Two-point control has different auxiliary/boundary equations.
    isnothing(f.control_points) || return r
    n,N=f.gas.n_species,length(f.grid)
    iv,il,im=n+2,n+3,n+4
    p=w.properties; z=f.grid
    for j in 1:N
        boundary=j==1 || j==N
        r[1,j]*=boundary || !isempty(f.fixed_temperature) ? 1000 : 1000/Arrhenius._flame_timescale
        r[iv,j]*=boundary ? Arrhenius._counterflow_Vscale : Arrhenius._counterflow_Vscale/Arrhenius._flame_timescale
        r[il,j]*=j==1 ? -1 : Arrhenius._counterflow_Lscale
        j<N && (r[im,j]/=-(z[j+1]-z[j]))
        if boundary
            # Cantera eliminates the largest species independently at each
            # inlet; native eliminates one species throughout the grid.
            flux=j==1 ? f.fuel_mass_flux : f.oxidizer_mass_flux
            yin=j==1 ? f.fuel_Y : f.oxidizer_Y
            for k in 1:n
                r[k+1,j]=j==1 ? flux*yin[k]-x[im,j]*x[k+1,j]-p.flux[k,1] :
                    p.flux[k,N-1]+x[im,j]*x[k+1,j]+flux*yin[k]
            end
            excess=argmax(@view x[2:n+1,j])
            r[excess+1,j]=1-sum(@view x[2:n+1,j])
        else
            for k in 1:n
                r[k+1,j]/=Arrhenius._flame_timescale
            end
            # Restore the interior species PDE replaced by native sum(Y)=1.
            k=f.dependent_species
            jm,jp=x[im,j]>0 ? (j-1,j) : (j,j+1)
            cell=.5*(z[j+1]-z[j-1]); dz=z[jp]-z[jm]
            r[k+1,j]=(f.gas.MW[k]*p.source[k,j]-(p.flux[k,j]-p.flux[k,j-1])/cell-
                x[im,j]*(x[k+1,jp]-x[k+1,jm])/dz)/p.rho[j]
            previous===nothing || (r[k+1,j]-=(x[k+1,j]-previous[k+1,j])/dt)
        end
    end
    return r
end

const ROWS = [
    ("7MPa regrid=0",            7e6, 0,  200,  :fixed_growth,     1.5),
    ("7MPa regrid=1",            7e6, 1,  200,  :fixed_growth,     1.5),
    ("7MPa regrid=3",            7e6, 3,  200,  :fixed_growth,     1.5),
    ("15MPa budget=200",        15e6, 3,  200,  :fixed_growth,     1.5),
    ("15MPa budget=1250",       15e6, 3,  1250, :fixed_growth,     1.5),
    ("15MPa residual-ratio",    15e6, 3,  1250, :residual_ratio,   1.5),
    ("15MPa newton-iterations", 15e6, 3,  1250, :newton_iterations,1.5),
]

function build_flame(MECH,P)
    gas = CreateSolution(MECH)
    mdot_f = 0.3
    # In the reference example both inlet .phase properties return the same
    # gas object; their densities are equal, despite different inlet T and X.
    # Preserve the example's effective boundary condition, mdot_o = 3*mdot_f.
    mdot_o = mdot_f * 3.0
    grid = collect(range(0.0, 30e-3; length=20))
    return CounterflowDiffusionFlame(gas; fuel="H2:1", oxidizer="O2:1",
        T_fuel=800.0, T_oxidizer=711.0, P, mdot_fuel=mdot_f, mdot_oxidizer=mdot_o,
        grid)
end

function run_rows(MECH,TAG;SAVE_DTS=false)
    dts_io = nothing
    SAVE_DTS && (dts_io = open("$(TAG)_dts.csv", "w"); println(dts_io, "row,step,dt"))

    println("row,ok,pts,Tmax,steps,jac,newton_iters,ts_failures,regrids,wall_s")
    for (label, P, regrid, budget, strategy, growth) in ROWS[parse(Int,get(ENV,"ROW_START","1")):end]
        f = build_flame(MECH,P)
        snap = (Arrhenius.COUNTERFLOW_STATS.newton_iterations,
                Arrhenius.COUNTERFLOW_STATS.jacobians,
                Arrhenius.COUNTERFLOW_STATS.time_steps,
                Arrhenius.COUNTERFLOW_STATS.time_step_failures,
                Arrhenius.COUNTERFLOW_STATS.regrids)
        call0 = length(Arrhenius.COUNTERFLOW_STATS.calls)
        dts = Float64[]
        cb = SAVE_DTS ? (dt -> push!(dts, dt)) : nothing
        t0 = time()
        ok = true
        try
            solve!(f; refine_grid=true, ratio=2.0, slope=0.06, curve=0.08, prune=0.02,
                max_time_steps=budget, time_step_regrid=regrid,
                time_step_growth_strategy=strategy, time_step_growth_factor=growth,
                jacobian_mode=:finite_difference, auto=false,
                initial_time_step=1e-5, loglevel=0, time_step_callback=cb)
        catch err
            ok = false
            println(stderr, "$label failed: $err")
        end
        wall = time() - t0
        s = Arrhenius.COUNTERFLOW_STATS
        steps = s.time_steps - snap[3]
        jac = s.jacobians - snap[2]
        nit = s.newton_iterations - snap[1]
        tsf = s.time_step_failures - snap[4]
        rg = s.regrids - snap[5]
        tmax = maximum(temperature(f))
        println("$label,$ok,$(length(f.grid)),$(round(tmax,digits=1)),$steps,$jac,$nit,$tsf,$rg,$(round(wall,digits=1))")
        calls = s.calls[call0+1:end]
        for (i, c) in enumerate(calls)
            println("  attempt $i: grid=$(Int(c[1])) dt_in=$(c[2]) newton=$(Int(c[3])) jac=$(Int(c[4])) steps=$(Int(c[5])) ts_fail=$(Int(c[6])) wall=$(round(c[7],digits=1))s conv=$(Bool(c[8]))")
        end
        if SAVE_DTS
            for (i, dt) in enumerate(dts)
                println(dts_io, "$label,$i,$dt")
            end
            flush(dts_io)
        end
        if ok
            r=similar(f.state); w=CounterflowWorkspace(f)
            physical_residual!(r,f,f.state,w)
            sumy=maximum(abs.(sum(mass_fractions(f);dims=1).-1))
            continuity=maximum(abs.(r[end,:]))
            println("  invariants sumY=",sumy," continuity=",continuity," physical_residual=",maximum(abs.(r)))
            npzwrite("$(TAG)_$(replace(label,r"[^a-zA-Z0-9]"=>"_")).npz",Dict("grid"=>f.grid,"state"=>f.state,"T"=>temperature(f),"Y"=>mass_fractions(f),"physical_residual"=>r,"species_names_utf8"=>Vector{UInt8}(join(f.gas.species_names,"\n"))))
        end
        flush(stdout)
        f = nothing
        GC.gc()
    end
    SAVE_DTS && close(dts_io)

end

if abspath(PROGRAM_FILE)==@__FILE__
    run_rows(abspath(ARGS[1]),ARGS[2];SAVE_DTS="--dts" in ARGS)
end
