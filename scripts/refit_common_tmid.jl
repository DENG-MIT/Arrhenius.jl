using Arrhenius
length(ARGS) in (2,3) || error("usage: julia --project scripts/refit_common_tmid.jl <input.yaml> <output.yaml> [Tmid]")
input, output = abspath.(ARGS[1:2])
tmid = length(ARGS)==3 ? parse(Float64,ARGS[3]) : nothing
report,_ = refit_nasa7_common_tmid(input,output;Tmid=tmid)
println("species_changed=$(report.species_changed) skipped=$(length(report.skipped_species)) Tmid=$(report.Tmid) shared_interval=$(report.shared_interval)")
for key in (:max_abs_delta_cp_R,:max_abs_delta_h_RT,:max_abs_delta_s_R,:max_abs_delta_g_RT)
    value=maximum(getproperty(change,key) for change in report.changes;init=0.0)
    println("$(key)=$(value)")
end
println("output=$(report.output)")