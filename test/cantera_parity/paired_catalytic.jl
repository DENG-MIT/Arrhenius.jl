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
source=replace(source,"Arrhenius.solve!(f)"=>"begin local nstart=length(f.grid); Arrhenius.solve!(f;ratio=3.,slope=.1,curve=.2); push!(catalytic_snapshots,(nstart,copy(f.grid),copy(f.state),copy(f.coverages))); end")
records=Dict()
for (label,path) in (("main",ENV["BASELINE_REFINE"]),
        ("branch",joinpath(dirname(pathof(Arrhenius)),"PremixedCounterflowFlames.jl")))
    Core.eval(Arrhenius,catalytic_refine_expr(path))
    global catalytic_snapshots=Any[]
    Arrhenius.reset!(Arrhenius.COUNTERFLOW_STATS)
    Base.include_string(Main,source,joinpath(dirname(pathof(Arrhenius)),"..","test/catalytic_flames_native.jl"))
    s=Arrhenius.COUNTERFLOW_STATS
    records[label]=(deepcopy(catalytic_snapshots),s.residual_evals)
    d=Dict{String,Any}("residual_evals"=>[s.residual_evals])
    for (i,(nstart,z,x,c)) in enumerate(catalytic_snapshots)
        d["$i/initial_points"]=[nstart]
        d["$i/grid"]=z;d["$i/state"]=x;d["$i/coverages"]=c
    end
    npzwrite(joinpath(out,label*".npz"),d)
end
@testset "catalytic dispatch branch parity" begin
    @test records["main"][2]==records["branch"][2]
    @test length(records["main"][1])==length(records["branch"][1])==2
    @test any(length(a[2])>a[1] for a in records["main"][1])
    for (a,b) in zip(records["main"][1],records["branch"][1])
        @test a==b
    end
end
println("CATALYTIC_BRANCH_PARITY snapshots=",length(records["main"][1])," residual_evals=",records["main"][2]," exact=true")
