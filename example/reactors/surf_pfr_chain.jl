# Prepare methane_pox_on_pt.yaml / Pt_surf with mechanism/export_surface.py first.
using Arrhenius
include("surface_solver.jl")
include("surface_pfr_setup.jl")

function surf_pfr_chain(path;output="surf_pfr_chain.csv",reactors=201)
    flow = methane_surface_pfr(path)
    result = surface_reactor_chain(flow,0.003;integrator=native_surface_bdf,reactors,
        reltol=1e-10,abstol=1e-19,maxiters=1_000_000,save_everystep=false)
    data = hcat(result.distance.*1e3,fill(flow.temperature-273.15,reactors),result.pressure./one_atm,
                transpose(result.mole_fractions),transpose(result.coverages))
    header = vcat(["Distance (mm)","T (C)","P (atm)"],flow.gas.species_names,
                   flow.surface.species_names[1:flow.surface.n_surface])
    open(output,"w") do io
        println(io,join(header,','))
        for row in eachrow(data)
            println(io,join(row,','))
        end
    end
    return result
end
if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS)>=1 || error("supply the prepared methane/Pt surface archive")
    surf_pfr_chain(ARGS[1];output=length(ARGS)>1 ? ARGS[2] : "surf_pfr_chain.csv")
end
