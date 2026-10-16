# Compare both policies in one Julia process; reuse the production adapter.
using Arrhenius, YAML, NPZ, LinearAlgebra
BLAS.set_num_threads(1)
root=joinpath(ENV["HOME"],"work/hotwin")
out=joinpath(root,"round2")
adapter=joinpath(get(ENV,"NH3_HARNESS_DIR","") == "" ? error("set NH3_HARNESS_DIR to the validation harness directory") : ENV["NH3_HARNESS_DIR"],"h/nh3_native_server.jl")
Base.include_string(Main,split(read(adapter,String),"\nwhile true\n")[1],adapter)
repo=normpath(joinpath(dirname(pathof(Arrhenius)),".."))
source=haskey(ENV,"POLICY_A_SOURCE") ? read(ENV["POLICY_A_SOURCE"],String) :
    read(`git -C $repo show 593d693:src/PremixedCounterflowFlames.jl`,String)
function policy!(label)
    text=label=="A" ? source : replace(source,
        "interpolate_velocity=false)"=>"interpolate_velocity=true)",
        "refinement_interpolation=:mass_flux)"=>"refinement_interpolation=:velocity)")
    for ex in Meta.parseall(text).args
        ex isa Expr || continue
        # The solve! definition carries a docstring macro.
        isref=ex.head==:function && startswith(string(ex.args[1]),"_refine_axisymmetric!")
        issolve=ex.head==:macrocall && occursin("function solve!(f::AbstractPremixedCounterflow",string(ex))
        (isref || issolve) && Core.eval(Arrhenius,ex)
    end
end
rows=split.(readlines(joinpath(out,"manifest.csv"))[2:end],",")
handles=Dict(label=>open(joinpath(out,label,"counts.tsv"),"w") for label in ("A","B"))
for io in values(handles)
    println(io,"id\tok\tjacobians\tnewton\tsteps\tfailures\tpoints\trecoveries")
end
try
    # H09 first so a wrong variant switch is detected before the full set.
    cases=vcat([["H09",joinpath(root,"H09.json"),""]],rows)
    for row in cases, label in ("A","B")
        policy!(label)
        id,spec,reference=strip.(row)
        target=joinpath(out,label,id*".npz")
        Arrhenius.reset!(Arrhenius.COUNTERFLOW_STATS)
        ok=true
        try
            Base.invokelatest(run_case,spec,target,1,"1")
        catch err
            ok=false
            open(target*".error","w") do io; showerror(io,err,catch_backtrace()); end
        end
        s=Arrhenius.COUNTERFLOW_STATS
        points=ok ? length(npzread(target)["grid"]) : 0
        println(handles[label],join((id,ok,s.jacobians,s.newton_iterations,s.time_steps,s.time_step_failures,points,s.recoveries),'\t'))
        flush(handles[label])
        println("ROUND2_DONE $label $id ok=$ok jac=$(s.jacobians) points=$points recoveries=$(s.recoveries)")
        flush(stdout)
    end
finally
    foreach(close,values(handles))
end
println("ROUND2_NATIVE_COMPLETE")
