# Export Arrhenius residual and banded Jacobian at stored Cantera fixed-T states.
# Arguments are repeated spec.json state.npz out.npz hscale groups.
using Arrhenius, YAML, NPZ, LinearAlgebra
BLAS.set_num_threads(1)
const H_SCALE = Ref(1.0)
function Arrhenius._flame_perturbation(f::Arrhenius.AbstractPremixedFlame,u,k,j)
    H_SCALE[] * 1e-7 * max(abs(u[k,j]),k==1 ? .1 : 1e-5)
end
function run_case(specpath,statepath,outfile,hscale)
    H_SCALE[] = parse(Float64,hscale)
    spec=YAML.load_file(specpath); d=npzread(statepath); out=outfile; id=spec["id"]
    gas=CreateSolution(spec["mechanism"]*".yaml")
    z=vec(d["grid"]); T=vec(d["T"]); Yraw=Matrix{Float64}(d["Y"]); U=vec(d["velocity"])
    Y=size(Yraw,1)==length(z) ? Yraw : size(Yraw,2)==length(z) ? permutedims(Yraw) : error("Y shape does not match grid")
    X=Dict(String(k)=>parse(Float64,v) for (k,v) in split.(split(spec["X"],","),":"))
    transport=spec["transport"]=="multicomponent" ? :multicomponent : :mixture_averaged
    data=transport==:multicomponent ? MultiTransportData(spec["mechanism"]*".yaml.multicomponent.npz",gas;mechanism=spec["mechanism"]*".yaml") : nothing
    f=BurnerFlame(gas;mdot=spec["mdot"],T=spec["T"],P=spec["P_atm"]*one_atm,X,
        width=spec["width"],transport_model=transport,multicomponent_data=data,
        discretization=:finite_difference,grid=z)
    set_temperature_profile!(f,Float64.(spec["profile_z"]),Float64.(spec["profile_T"]);relative=false)
    length(f.grid)==length(z) || error("profile insertion changed stored grid")
    n=gas.n_species
    meanMW=1 ./ (Y*(1 ./ gas.MW))
    rho=spec["P_atm"]*one_atm .* meanMW ./ (Arrhenius.R .* T)
    f.state[1,:]=T./1000
    f.state[2:n+1,:]=permutedims(Y)
    f.state[end,:]=rho.*U
    r=similar(f.state); w=Arrhenius.FlameWorkspace(f)
    Arrhenius.flame_residual!(r,f,f.state,w)
    J=Arrhenius._flame_jacobian(f,f.state,w,r;analytic=true)
    B,N=size(f.state); m=B*N; bw=2B-1
    bands=zeros(2bw+1,m)
    for col in 1:m
        for row in max(1,col-bw):min(m,col+bw)
            bands[bw+1+row-col,col]=J[2bw+1+row-col,col]
        end
    end
    npzwrite(out,Dict("r_arr"=>vec(r),"bands_arr"=>bands,"bw"=>[bw],"rho"=>rho,
        "T"=>T,"velocity"=>U,"scale"=>[H_SCALE[]],"dependent_species"=>[f.dependent_species],"inlet_Y"=>f.inlet_Y))
    println("case=$id points=$N dofs=$m band=$bw hscale=$(H_SCALE[]) residual_inf=$(norm(r,Inf)) file=$out")
end

length(ARGS)%4==0 && !isempty(ARGS) || error("provide spec/state/output/hscale groups")
for i in 1:4:length(ARGS)
    run_case(ARGS[i],ARGS[i+1],ARGS[i+2],ARGS[i+3])
end