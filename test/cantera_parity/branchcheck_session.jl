# Run a branchcheck's native cases in one persistent process. Cantera jobs and
# gates use the existing harness; BASELINE_REFINE supplies the main method.
using Arrhenius, YAML, NPZ
harness=joinpath(get(ENV,"NH3_HARNESS_DIR","") == "" ? error("set NH3_HARNESS_DIR to the validation harness directory") : ENV["NH3_HARNESS_DIR"],"h/nh3_native_server.jl")
Base.include_string(Main,split(read(harness,String),"\nwhile true\n")[1],harness)
function refine_expr(path)
    expr=Meta.parseall(read(path,String))
    method=only(filter(x->x isa Expr && x.head==:function &&
        startswith(string(x.args[1]),"_refine_axisymmetric!"),expr.args))
    # Older baseline methods predate this keyword; ignoring it preserves their
    # original interpolation while allowing the current solver to call them.
    parameters=only(filter(x->x isa Expr && x.head==:parameters,method.args[1].args))
    any(x->x isa Expr && x.head==:kw && x.args[1]==:interpolate_velocity,parameters.args) ||
        push!(parameters.args,Expr(:kw,:interpolate_velocity,false))
    return method
end
branch_out=ENV["BRANCHCHECK_OUT"]
manifest=split.(readlines(joinpath(branch_out,"manifest.csv"))[2:end],",")
main_expr=refine_expr(ENV["BASELINE_REFINE"])
branch_expr=refine_expr(joinpath(dirname(pathof(Arrhenius)),"PremixedCounterflowFlames.jl"))
for (label,expr) in (("main",main_expr),("branch",branch_expr))
    label in split(get(ENV,"BRANCHCHECK_MODES","main,branch"),",") || continue
    Core.eval(Arrhenius,expr)
    mkpath(joinpath(branch_out,label))
    open(joinpath(branch_out,label,"counts.tsv"),"w") do io
        println(io,"id\tok\tjacobians\tnewton\tsteps\tfailures\tpoints")
        for row in manifest
            id,spec,reference=strip.(row)
            out=joinpath(branch_out,label,id*".npz")
            Arrhenius.reset!(Arrhenius.COUNTERFLOW_STATS)
            ok=true
            try
                Base.invokelatest(run_case,spec,out,1,"1")
            catch err
                ok=false
                open(out*".error","w") do eio; showerror(eio,err); end
                println("BRANCHCHECK_ERROR $label $id ",sprint(showerror,err))
            end
            s=Arrhenius.COUNTERFLOW_STATS
            # The adapter performs one warm and one measured identical solve.
            # Raw counts include both; do not mix these with per-solve counts.
            # Other flame families do not increment CounterflowStats.
            points=ok ? length(npzread(out)["grid"]) : 0
            println(io,join((id,ok,s.jacobians,s.newton_iterations,s.time_steps,s.time_step_failures,points),'\t'))
            flush(io)
            println("BRANCHCHECK_DONE $label $id")
        end
    end
end
