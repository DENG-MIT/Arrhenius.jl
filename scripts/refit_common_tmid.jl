using Arrhenius
length(ARGS) in (2,3,4) || error("usage: julia --project scripts/refit_common_tmid.jl <input.yaml> <output.yaml> [Tmid] [--force]")
input, output = abspath.(ARGS[1:2])
rest=ARGS[3:end]
force="--force" in rest
numeric=filter(x->x!="--force",rest)
length(numeric)<=1 || error("provide at most one Tmid")
tmid=isempty(numeric) ? nothing : parse(Float64,first(numeric))
report,_ = refit_nasa7_common_tmid(input,output;Tmid=tmid,force)
println("species_changed=$(report.species_changed) refused=$(report.species_refused) skipped=$(length(report.skipped_species)) forced=$(force) Tmid=$(report.Tmid) shared_interval=$(report.shared_interval)")
for key in (:max_abs_delta_cp_R,:max_abs_delta_h_RT,:max_abs_delta_s_R,:max_abs_delta_g_RT)
    value=maximum(getproperty(change,key) for change in report.changes;init=0.0)
    println("$(key)=$(value)")
end
println("output=$(report.output)")
