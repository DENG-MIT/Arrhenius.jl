# Prepare: python mechanism/export_surface.py methane_pox_on_pt.yaml Pt_surf --output pox.surface.npz
using Arrhenius
include("surface_solver.jl")
include("surface_flow_solver.jl")
include("surface_pfr_setup.jl")

function surf_pfr(path;output="surf_pfr.csv")
    flow = methane_surface_pfr(path)
    distance = collect(range(0.0,0.003;length=301))
    solution = solve_surface_flow(flow,0.003;integrator=native_surface_flow_bdf,
        coverage_integrator=native_surface_bdf,reltol=1e-9,abstol=1e-15,
        saveat=distance,tstops=distance[2:end],dt=1e-10,maxiters=1_000_000)
    state = reduce(hcat,solution.u)
    ng = flow.gas.n_species
    X = reduce(hcat,[Y2X(flow.gas,state[1:ng,k]) for k in axes(state,2)])
    pressure = [surface_flow_properties(flow,state[1:ng,k]).pressure for k in axes(state,2)]
    data = hcat(distance.*1e3,fill(flow.temperature-273.15,length(distance)),pressure./one_atm,
                transpose(X),transpose(state[ng+1:end,:]))
    header = vcat(["Distance (mm)","T (C)","P (atm)"],flow.gas.species_names,
                   flow.surface.species_names[1:flow.surface.n_surface])
    open(output,"w") do io
        println(io,join(header,','))
        for row in eachrow(data)
            println(io,join(row,','))
        end
    end
    return (;flow,solution,mole_fractions=X,pressure)
end
if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS)>=1 || error("supply the prepared methane/Pt surface archive")
    surf_pfr(ARGS[1];output=length(ARGS)>1 ? ARGS[2] : "surf_pfr.csv")
end
