# Paired catalytic dispatch check using the existing catalytic fixture.
# BASELINE_REFINE: main source; CATALYTIC_CHECK_OUT: output directory.
# Tight refinement criteria exercise midpoint insertion in both sequences.
using Arrhenius, Test, NPZ
function catalytic_refine_expr(path)
    expr=Meta.parseall(read(path,String))
    only(filter(x->x isa Expr && x.head==:function && startswith(string(x.args[1]),"_refine_axisymmetric!"),expr.args))
end
out=ENV["CATALYTIC_CHECK_OUT"];mkpath(out)
source=read(joinpath(dirname(pathof(Arrhenius)),"..","test/catalytic_flames_native.jl"),String)
source=replace(source,"Arrhenius.solve!(f)"=>"begin local nstart=length(f.grid); Arrhenius.solve!(f;ratio=3.,slope=.1,curve=.2); push!(catalytic_snapshots,(nstart,copy(f.grid),copy(f.state),copy(f.coverages),Arrhenius.counterflow_residual!(similar(f.state),f))); end")
records=Dict()
for (label,path) in (("main",ENV["BASELINE_REFINE"]),
        ("branch",joinpath(dirname(pathof(Arrhenius)),"PremixedCounterflowFlames.jl")))
    Core.eval(Arrhenius,catalytic_refine_expr(path))
    global catalytic_snapshots=Any[]
    Arrhenius.reset!(Arrhenius.COUNTERFLOW_STATS)
    Base.include_string(Main,source,joinpath(dirname(pathof(Arrhenius)),"..","test/catalytic_flames_native.jl"))
    s=Arrhenius.COUNTERFLOW_STATS
    records[label]=(deepcopy(catalytic_snapshots),s.residual_evals-length(catalytic_snapshots))
    d=Dict{String,Any}("residual_evals"=>[s.residual_evals-length(catalytic_snapshots)],
        "diagnostic_residual_evals"=>[length(catalytic_snapshots)])
    for (i,(nstart,z,x,c,r)) in enumerate(catalytic_snapshots)
        d["$i/initial_points"]=[nstart]
        d["$i/grid"]=z;d["$i/state"]=x;d["$i/coverages"]=c;d["$i/residual"]=r
    end
    npzwrite(joinpath(out,label*".npz"),d)
end
@testset "catalytic dispatch branch parity" begin
    @test records["main"][2]==records["branch"][2]
    @test length(records["main"][1])==length(records["branch"][1])==2
    @test any(length(a[2])>a[1] for a in records["main"][1])
    for (a,b) in zip(records["main"][1],records["branch"][1])
        @test a[1]==b[1]
        @test a[2]==b[2] # positions, not just point counts
        @test a[3] ≈ b[3]
        @test a[4]==b[4]
        # The residual implementation is identical for the two refinement
        # policies. Evaluate it at both states against the existing Newton
        # tolerance (_counterflow_newton!, default tolerance=1e-8).
        @test maximum(abs,a[5]) < 1e-8
        @test maximum(abs,b[5]) < 1e-8
        println("CATALYTIC_COMPARE state_max_abs=",maximum(abs,a[3]-b[3]),
            " native_residual_main=",maximum(abs,a[5]),
            " native_residual_branch=",maximum(abs,b[5]))
    end
end
println("CATALYTIC_BRANCH_PARITY snapshots=",length(records["main"][1])," residual_evals=",records["main"][2]," numerical=true")
