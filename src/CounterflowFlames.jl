# Governing equations and boundary discretization:
# https://cantera.org/dev/reference/onedim/governing-equations.html
# Flow1D.cpp, Boundary1D.cpp and CounterflowDiffusionFlame.set_initial_guess
# at Cantera commit 726522be4e2a13454d8415b7ef799d621f665cf3.
"Axisymmetric opposed fuel/oxidizer jets with optional optically thin gas radiation."
mutable struct CounterflowDiffusionFlame{G<:Solution}
    gas::G
    grid::Vector{Float64}
    pressure::Float64
    fuel_temperature::Float64
    oxidizer_temperature::Float64
    fuel_Y::Vector{Float64}
    oxidizer_Y::Vector{Float64}
    fuel_mass_flux::Float64
    oxidizer_mass_flux::Float64
    state::Matrix{Float64}
    anchor::Int
    dependent_species::Int
    converged::Bool
    transport_model::Symbol
    soret_enabled::Bool
    multicomponent_data::Union{Nothing,MultiTransportData}
    flux_gradient_basis::Symbol
    fixed_temperature::Vector{Float64}
    radiation_enabled::Bool
    boundary_emissivities::Tuple{Float64,Float64}
    # (left_z, left_T, right_z, right_T) when two-point control is enabled.
    control_points::Union{Nothing,NTuple{4,Float64}}
end

# Scales only condition the algebraic system; physical values are returned below.
const _counterflow_Vscale = 100.0
const _counterflow_Lscale = 10000.0
_flame_refine_threshold(f::CounterflowDiffusionFlame,k) =
    k == 1 ? sqrt(eps(Float64))/1000 :
    k == size(f.state,1)-2 ? sqrt(eps(Float64))/_counterflow_Vscale : sqrt(eps(Float64))
_flame_refine_active(f::CounterflowDiffusionFlame,k) =
    k != size(f.state,1)-1 && (isnothing(f.control_points) || k != f.gas.n_species+2)
function _flame_refine_protect!(keep,f::CounterflowDiffusionFlame)
    isnothing(f.control_points) && return nothing
    for z in (f.control_points[1],f.control_points[3])
        j=findfirst(==(z),f.grid)
        isnothing(j) && throw(ArgumentError("control coordinate must remain on the grid"))
        keep[j]=1
    end
    nothing
end
const _counterflow_stefan_boltzmann = 2*pi^5*(1.380649e-23)^4 /
    (15*(6.62607015e-34)^3*(299792458.)^2)
_counterflow_erf(x) = ccall((:erf,Base.Math.libm),Float64,(Float64,),x)

"""
    CounterflowDiffusionFlame(gas; fuel, oxidizer, mdot_fuel, mdot_oxidizer,
                             T_fuel=300, T_oxidizer=300, P=one_atm, width=.02)

Construct a native axisymmetric counterflow diffusion flame. Both inlet mass
fluxes are positive magnitudes [kg/(m² s)], directed into the domain. The axial
velocity is positive at the fuel inlet and negative at the oxidizer inlet.
Initialization uses a native stoichiometric HP equilibrium and an error-function
mixing layer. Supports mixture-averaged neutral ideal gas, with optional gray
optically thin H2O/CO2 radiation (`radiation=true`).
"""
function CounterflowDiffusionFlame(gas::Solution; fuel,oxidizer,mdot_fuel,mdot_oxidizer,
        T_fuel=300.,T_oxidizer=300.,P=one_atm,width=.02,grid=nothing,
        radiation=false,boundary_emissivities=(0.,0.))
    all(x->isfinite(x)&&x>0,(P,mdot_fuel,mdot_oxidizer,width)) ||
        throw(ArgumentError("positive finite pressure, width and inlet mass fluxes required"))
    all(x->isfinite(x)&&200<=x<=6000,(T_fuel,T_oxidizer)) ||
        throw(ArgumentError("inlet temperatures must lie in 200–6000 K"))
    gas.trans.poly_order == 5 || throw(ArgumentError("native transport sidecar required"))
    gas.thermo isa IdealGasThermo || throw(ArgumentError("counterflow requires ideal-gas thermodynamics"))
    z=isnothing(grid) ? width.*[0,.2,.4,.6,.8,1] : Float64.(grid)
    length(z)>=5 && all(isfinite,z) && all(>(0),diff(z)) ||
        throw(ArgumentError("at least five increasing grid points required"))
    xf=mole_fractions(gas,fuel); xo=mole_fractions(gas,oxidizer)
    yf=xf.*gas.MW./dot(xf,gas.MW); yo=xo.*gas.MW./dot(xo,gas.MW)
    demand=_oxygen_demand(gas)
    df=dot(demand,yf./gas.MW); dox=dot(demand,yo./gas.MW)
    df>0 && dox<0 || throw(ArgumentError("reactive fuel and oxygen-bearing oxidizer required"))
    zst=-dox/(df-dox)
    yst=zst.*yf.+(1-zst).*yo
    xst=mole_fractions(gas,yst;basis=:mass)
    eq=equilibrate(gas;T=.5*(T_fuel+T_oxidizer),P,X=xst,mode=:HP)
    rhof=P*dot(xf,gas.MW)/(R*T_fuel); rhoo=P*dot(xo,gas.MW)/(R*T_oxidizer)
    uf=mdot_fuel/rhof; uo=mdot_oxidizer/rhoo
    width=z[end]-z[1]
    strain=(uf+uo)/width
    curvature=-.5*(rhof+rhoo)*strain^2
    stagnation=sqrt(mdot_fuel*uf)*width/(sqrt(mdot_fuel*uf)+sqrt(mdot_oxidizer*uo))
    transport=TransportWorkspace(gas)
    mixture_transport!(transport,gas,P,eq.T,eq.X)
    ioxygen=findfirst(s->lowercase(s)=="o2",gas.species_names)
    isnothing(ioxygen) && throw(ArgumentError("O2 species required for mixing-layer initialization"))
    mixing_scale=sqrt(strain/(2*transport.diffusion[ioxygen]))
    n=gas.n_species
    state=zeros(n+4,length(z))
    for j in eachindex(z)
        x=z[j]-z[1]; fraction=x/width
        Z=.5*(1-_counterflow_erf(mixing_scale*(x-stagnation)))
        if Z>zst
            blend=(Z-zst)/(1-zst)
            Y=eq.Y.+blend.*(yf.-eq.Y)
            T=eq.T+blend*(T_fuel-eq.T)
        else
            blend=Z/zst
            Y=yo.+blend.*(eq.Y.-yo)
            T=T_oxidizer+blend*(eq.T-T_oxidizer)
        end
        state[1,j]=T/1000
        state[2:n+1,j]=Y
        spread=x<stagnation ? strain*x/stagnation : strain*(width-x)/(width-stagnation)
        state[n+2,j]=spread/_counterflow_Vscale
        state[n+3,j]=curvature/_counterflow_Lscale
        rho=P/(R*T*sum(Y./gas.MW))
        state[end,j]=rho*((1-fraction)*uf-fraction*uo)
    end
    state[1,1]=T_fuel/1000; state[1,end]=T_oxidizer/1000
    state[end,1]=mdot_fuel; state[end,end]=-mdot_oxidizer
    flame=CounterflowDiffusionFlame(gas,z,Float64(P),Float64(T_fuel),Float64(T_oxidizer),
        yf,yo,Float64(mdot_fuel),Float64(mdot_oxidizer),state,1,argmax(yst),false,
        :mixture_averaged,false,nothing,:mole,Float64[],false,(0.,0.),nothing)
    return set_radiation!(flame,radiation;boundary_emissivities)
end

"""
    set_radiation!(flame, enabled=true; boundary_emissivities=(0,0))

Enable or disable the gray optically thin H2O/CO2 radiation model. Emissivities
of the fuel and oxidizer boundaries must lie in [0,1]. The current solution is
retained for continuation; call `solve!(flame; refine_grid=false)` to repeat the
official nonradiating-to-radiating calculation on the same grid.
"""
function set_radiation!(f::CounterflowDiffusionFlame,enabled::Bool=true;
        boundary_emissivities=f.boundary_emissivities)
    length(boundary_emissivities)==2 && all(e->isfinite(e)&&0<=e<=1,boundary_emissivities) ||
        throw(ArgumentError("two boundary emissivities in [0,1] required"))
    f.radiation_enabled=enabled
    f.boundary_emissivities=(Float64(boundary_emissivities[1]),Float64(boundary_emissivities[2]))
    f.converged=false
    return f
end

# Nearest grid point to the first (fromleft=true) or last crossing of target.
function _counterflow_crossing(T,z,target,fromleft)
    N=length(z)
    for j in (fromleft ? (1:N-1) : (N-1:-1:1))
        a=T[j]-target; b=T[j+1]-target
        if a==0 || b==0 || a*b<0
            return fromleft ? (abs(a)<abs(b) ? j : j+1) : (abs(b)<abs(a) ? j+1 : j)
        end
    end
    return nothing
end

"""
    set_two_point_control!(flame; temperature, decrement=0)

Enable two-point temperature control: fixes the temperature [K] at the interior
grid points nearest the first (fuel side) and last (oxidizer side) crossing of
`temperature`, recording the actual selected grid temperatures minus `decrement`
[K]. A positive oxidizer mass-flux auxiliary row [kg/(m² s)] is inserted before
the spread-rate row on first enable; the inlet mass fluxes become unknowns.
`temperature` must lie in 200–6000 K and exceed both inlet temperatures.
"""
function set_two_point_control!(f::CounterflowDiffusionFlame;temperature,decrement=0.)
    isempty(f.fixed_temperature) || throw(ArgumentError("two-point control conflicts with fixed_temperature"))
    isfinite(temperature) && 200<=temperature<=6000 ||
        throw(ArgumentError("control temperature must lie in 200–6000 K"))
    isfinite(decrement) && decrement>=0 ||
        throw(ArgumentError("nonnegative finite decrement required"))
    temperature>f.fuel_temperature && temperature>f.oxidizer_temperature ||
        throw(ArgumentError("control temperature must exceed both inlet temperatures"))
    T=1000 .* vec(f.state[1,:]); z=f.grid; N=length(z)
    jl=_counterflow_crossing(T,z,temperature,true)
    jr=_counterflow_crossing(T,z,temperature,false)
    (!isnothing(jl) && !isnothing(jr) && 2<=jl && jl<jr && jr<=N-1) ||
        throw(ArgumentError("two distinct interior control crossings required"))
    left_T,right_T=T[jl]-decrement,T[jr]-decrement
    200<=left_T<=6000 && 200<=right_T<=6000 &&
        left_T>f.fuel_temperature && right_T>f.oxidizer_temperature ||
        throw(ArgumentError("decrement must leave control temperatures above their inlets and within 200–6000 K"))
    n=f.gas.n_species
    if isnothing(f.control_points)
        f.state=vcat(f.state[1:n+1,:],fill(f.oxidizer_mass_flux,1,N),f.state[n+2:end,:])
    end
    f.control_points=(z[jl],left_T,z[jr],right_T)
    f.converged=false
    return f
end

"""
    disable_two_point_control!(flame)

Remove the auxiliary mass-flux row, store the current axial mass-flux endpoints
[kg/(m² s)] as the inlet fluxes and restore inlet-flux boundary conditions.
Enable/disable without solving reproduces the initial state bit-exactly.
"""
function disable_two_point_control!(f::CounterflowDiffusionFlame)
    isnothing(f.control_points) && return f
    n=f.gas.n_species
    f.fuel_mass_flux=f.state[end,1]
    f.oxidizer_mass_flux=-f.state[end,end]
    f.state=f.state[[1:n+1;n+3:n+5],:]
    f.control_points=nothing
    f.converged=false
    return f
end

@inline function _counterflow_planck_polynomial(T,coefficients)
    t=1000/T
    value=0.
    @inbounds for n in 0:5
        value+=coefficients[n+1]*t^n
    end
    return value/one_atm
end

# Liu–Rogg gray optically thin approximation, with RADCAL/TNF Planck fits.
# Same coefficients and sign convention as Cantera Flow1D::computeRadiation.
function _counterflow_radiation!(loss,absorption,f,u,X)
    fill!(loss,0); fill!(absorption,0)
    f.radiation_enabled || return
    ico2=findfirst(==("CO2"),f.gas.species_names)
    ih2o=findfirst(==("H2O"),f.gas.species_names)
    cco2=(18.741,-121.310,273.500,-194.050,56.310,-5.8169)
    ch2o=(-.23093,-1.12390,9.41530,-2.99880,.51382,-1.86840e-5)
    sigma=_counterflow_stefan_boltzmann
    # Inlet temperatures are prescribed Dirichlet data, avoiding artificial
    # nonlocal Jacobian couplings from perturbing their constrained unknowns.
    boundary=sigma*(f.boundary_emissivities[1]*f.fuel_temperature^4+
        f.boundary_emissivities[2]*f.oxidizer_temperature^4)
    for j in axes(u,2)
        T=1000*u[1,j]
        kp=0.
        isnothing(ih2o) || (kp+=f.pressure*X[ih2o,j]*_counterflow_planck_polynomial(T,ch2o))
        isnothing(ico2) || (kp+=f.pressure*X[ico2,j]*_counterflow_planck_polynomial(T,cco2))
        absorption[j]=kp
        # Cantera leaves the last inlet entry zero; no energy equation uses it.
        j<size(u,2) && (loss[j]=2*kp*(2*sigma*T^4-boundary))
    end
end

"""
Return `(heat_loss, planck_absorption)` at all grid points, in W/m³ and 1/m.
Positive heat loss cools the gas. Values are zero when radiation is disabled;
the last inlet heat-loss entry is zero, matching Cantera's diagnostic convention.
"""
function radiation_source(f::CounterflowDiffusionFlame)
    n,N=f.gas.n_species,length(f.grid)
    X=zeros(n,N)
    for j in 1:N
        denominator=sum(max(f.state[k+1,j],0)/f.gas.MW[k] for k in 1:n)
        for k in 1:n
            X[k,j]=max(f.state[k+1,j],0)/f.gas.MW[k]/denominator
        end
    end
    loss=zeros(N); absorption=zeros(N)
    _counterflow_radiation!(loss,absorption,f,f.state,X)
    return (heat_loss=loss,planck_absorption=absorption)
end
radiative_heat_loss(f::CounterflowDiffusionFlame)=radiation_source(f).heat_loss

temperature(f::CounterflowDiffusionFlame)=1000 .* vec(f.state[1,:])
mass_fractions(f::CounterflowDiffusionFlame)=copy(f.state[2:f.gas.n_species+1,:])
spread_rate(f::CounterflowDiffusionFlame)=_counterflow_Vscale .* vec(f.state[end-2,:])
pressure_curvature(f::CounterflowDiffusionFlame)=_counterflow_Lscale .* vec(f.state[end-1,:])
"Whether the computed temperature rise is less than 10 K."
extinct(f::CounterflowDiffusionFlame)=maximum(temperature(f))-max(f.fuel_temperature,f.oxidizer_temperature)<10
function velocity(f::CounterflowDiffusionFlame)
    return [f.state[end,j]*R*(1000*f.state[1,j])*
        sum(f.state[k+1,j]/f.gas.MW[k] for k in 1:f.gas.n_species)/f.pressure
        for j in eachindex(f.grid)]
end

# Allocate the existing shared property storage without initializing a premixed flame.
struct _CounterflowPropertyAdapter{G} <: AbstractPremixedFlame
    gas::G
    grid::Vector{Float64}
    transport_model::Symbol
    multicomponent_data::Union{Nothing,MultiTransportData}
    soret_enabled::Bool
    flux_gradient_basis::Symbol
end
struct CounterflowWorkspace{F}
    properties::F
    viscosity::Vector{Float64}
    band::Matrix{Float64}
    perturbed::Matrix{Float64}
    residual_perturbed::Matrix{Float64}
    steps::Vector{Float64}
    radiative_heat_loss::Vector{Float64}
    planck_absorption::Vector{Float64}
end
function CounterflowWorkspace(f::CounterflowDiffusionFlame)
    adapter=_CounterflowPropertyAdapter(f.gas,f.grid,f.transport_model,
        f.multicomponent_data,f.soret_enabled,f.flux_gradient_basis)
    B,N=size(f.state)
    return CounterflowWorkspace(FlameWorkspace(adapter),zeros(N-1),zeros(6B-2,B*N),
        zeros(B,N),zeros(B,N),zeros(N),zeros(N),zeros(N))
end

"Evaluate continuity, radial momentum, pressure curvature, species and energy."
function counterflow_residual!(residual,f::CounterflowDiffusionFlame,u=f.state,
        w=CounterflowWorkspace(f);previous=nothing,dt=Inf,update_transport=true,
        nodes=eachindex(f.grid))
    p=w.properties
    _flame_properties!(p,f,u;update_transport,nodes)
    _counterflow_radiation!(w.radiative_heat_loss,w.planck_absorption,f,u,p.X)
    n,N=f.gas.n_species,length(f.grid)
    z,MW=f.grid,f.gas.MW
    ctrl=f.control_points
    iq=isnothing(ctrl) ? 0 : n+2
    iv,il,im=size(u,1)-2,size(u,1)-1,size(u,1)
    jl=jr=0
    if !isnothing(ctrl)
        jl=findfirst(==(ctrl[1]),z); jr=findfirst(==(ctrl[3]),z)
        (isnothing(jl) || isnothing(jr)) &&
            error("two-point control coordinates are not on the current grid")
    end
    if update_transport
        for j in 1:N-1
            Tmid=500*(u[1,j]+u[1,j+1])
            for k in 1:n
                p.ymid[k]=max(0.,.5*(u[k+1,j]+u[k+1,j+1]))
            end
            denominator=sum(p.ymid./MW)
            p.xmid .= p.ymid./MW./denominator
            w.viscosity[j],_=mixture_transport!(p.transport,f.gas,f.pressure,Tmid,p.xmid)
        end
    end
    @inbounds for j in 1:N
        if j<N
            residual[im,j]=u[im,j+1]-u[im,j]+(z[j+1]-z[j])*
                _counterflow_Vscale*(p.rho[j+1]*u[iv,j+1]+p.rho[j]*u[iv,j])
        else
            residual[im,j]=u[im,j]+(isnothing(ctrl) ? f.oxidizer_mass_flux : u[iq,j])
        end
        if isnothing(ctrl)
            residual[il,j]=j==1 ? u[im,j]-f.fuel_mass_flux : u[il,j]-u[il,j-1]
        elseif j==1
            residual[il,j]=u[il,2]-u[il,1]
        elseif j==jl
            residual[il,j]=u[1,j]-ctrl[2]/1000
        elseif j<jl
            residual[il,j]=u[il,j+1]-u[il,j]
        else
            residual[il,j]=u[il,j]-u[il,j-1]
        end
        if !isnothing(ctrl)
            if j==jr
                residual[iq,j]=u[1,j]-ctrl[4]/1000
            elseif j<jr
                residual[iq,j]=u[iq,j+1]-u[iq,j]
            else
                residual[iq,j]=u[iq,j]-u[iq,j-1]
            end
        end
        if j==1 || j==N
            residual[iv,j]=u[iv,j]
            residual[1,j]=u[1,j]-(j==1 ? f.fuel_temperature : f.oxidizer_temperature)/1000
            fuel_flux=isnothing(ctrl) ? f.fuel_mass_flux : u[im,1]
            oxidizer_flux=isnothing(ctrl) ? f.oxidizer_mass_flux : u[iq,N]
            for k in 1:n
                residual[k+1,j]=j==1 ? fuel_flux*f.fuel_Y[k]-u[im,j]*u[k+1,j]-p.flux[k,1] :
                    p.flux[k,N-1]+u[im,j]*u[k+1,j]+oxidizer_flux*f.oxidizer_Y[k]
            end
        else
            left=z[j]-z[j-1]; right=z[j+1]-z[j]; cell=.5*(left+right)
            jm,jp=u[im,j]>0 ? (j-1,j) : (j,j+1)
            dz=z[jp]-z[jm]
            V=_counterflow_Vscale*u[iv,j]
            dV=_counterflow_Vscale*(u[iv,jp]-u[iv,jm])/dz
            shear=_counterflow_Vscale*(w.viscosity[j]*(u[iv,j+1]-u[iv,j])/right-
                w.viscosity[j-1]*(u[iv,j]-u[iv,j-1])/left)/cell
            residual[iv,j]=_flame_timescale/(p.rho[j]*_counterflow_Vscale)*
                (shear-_counterflow_Lscale*u[il,j]-u[im,j]*dV-p.rho[j]*V^2)
            dT=1000*(u[1,jp]-u[1,jm])/dz
            cpmean=0.; chemical=0.; enthalpyflux=0.
            for k in 1:n
                cpmean+=u[k+1,j]*p.cp[k,j]/MW[k]
                chemical+=p.h[k,j]*p.source[k,j]
                enthalpyflux+=.5*(p.flux[k,j-1]+p.flux[k,j])*(p.h[k,jp]-p.h[k,jm])/(dz*MW[k])
                residual[k+1,j]=_flame_timescale/p.rho[j]*(MW[k]*p.source[k,j]-
                    (p.flux[k,j]-p.flux[k,j-1])/cell-u[im,j]*(u[k+1,jp]-u[k+1,jm])/dz)
            end
            conduction=1000*(p.conductivity[j]*(u[1,j+1]-u[1,j])/right-
                p.conductivity[j-1]*(u[1,j]-u[1,j-1])/left)/cell
            residual[1,j]=_flame_timescale/(1000*p.rho[j]*cpmean)*
                (conduction-chemical-enthalpyflux-u[im,j]*cpmean*dT-w.radiative_heat_loss[j])
            if previous !== nothing
                # Under two-point control the radial momentum and auxiliary
                # mass-flux equations are algebraic; T and species retain time terms.
                for k in 1:(isnothing(ctrl) ? n+2 : n+1)
                    residual[k,j]-=_flame_timescale/dt*(u[k,j]-previous[k,j])
                end
            end
        end
        if !isempty(f.fixed_temperature)
            residual[1,j]=u[1,j]-f.fixed_temperature[j]/1000
        end
        residual[f.dependent_species+1,j]=sum(@view(u[2:n+1,j]))-1
    end
    return residual
end

function _counterflow_jacobian!(f,u,w,r;previous=nothing,dt=Inf)
    B,N=size(u); kl=ku=2B-1
    fill!(w.band,0)
    w.perturbed .= u
    p=w.properties
    saved=(copy(p.X),copy(p.rho),copy(p.cp),copy(p.h),copy(p.source))
    for k in 1:B, color in 1:3
        for j in color:3:N
            w.steps[j]=1e-7*max(abs(u[k,j]),k==1 ? .1 : k<=f.gas.n_species+1 ? 1e-5 : .01)
            w.perturbed[k,j]=u[k,j]+w.steps[j]
        end
        counterflow_residual!(w.residual_perturbed,f,w.perturbed,w;
            previous,dt,update_transport=false,nodes=color:3:N)
        for j in color:3:N
            column=(j-1)*B+k
            for jj in max(1,j-1):min(N,j+1),kk in 1:B
                row=(jj-1)*B+kk
                w.band[kl+ku+1+row-column,column]=(w.residual_perturbed[kk,jj]-r[kk,jj])/w.steps[j]
            end
            w.perturbed[k,j]=u[k,j]
        end
        p.X .= saved[1]; p.rho .= saved[2]; p.cp .= saved[3]; p.h .= saved[4]; p.source .= saved[5]
    end
    return w.band
end

function _counterflow_residual_scales!(scales,J,u,bw)
    fill!(scales,0.)
    n=length(u)
    for j in 1:n, i in max(1,j-bw):min(n,j+bw)
        scales[i]+=abs(J[2bw+1+i-j,j])*abs(u[j])
    end
    for i in eachindex(scales)
        scales[i]=max(1.,scales[i])
    end
    scales
end

function _counterflow_constraints_converged(f,r,tol)
    B,N=size(r)
    maximum(abs,@view(r[B-1:B,:]))<=tol || return false
    maximum(abs,@view(r[f.dependent_species+1,:]))<=min(tol,1e-9) || return false
    maximum(abs,@view(r[:,1]))<=tol && maximum(abs,@view(r[:,N]))<=tol || return false
    isnothing(f.control_points) || maximum(abs,@view(r[f.gas.n_species+2,:]))<=tol
end

function _counterflow_newton!(f,w;previous=nothing,dt=Inf,maxiters=45,tolerance=1e-8,loglevel=0)
    u=f.state; r=similar(u); trial=similar(u); rt=similar(u)
    n=f.gas.n_species; bandwidth=2*size(u,1)-1
    pivots=LinearAlgebra.BlasInt[]; age=5; contraction=Inf
    correction_enabled=!isnothing(f.control_points)
    correction_enabled && (age=21)
    weights=correction_enabled ? _flame_correction_weights(u,previous!==nothing) : Float64[]
    trial_correction=correction_enabled ? Vector{Float64}(undef,length(u)) : Float64[]
    residual_scales=correction_enabled ? Vector{Float64}(undef,length(u)) : Float64[]
    residual_valid = false
    for iteration in 1:maxiters
        if !residual_valid
            counterflow_residual!(r,f,u,w;previous,dt)
        end
        residual_valid = false
        residualnorm=norm(r,Inf)
        residualnorm<tolerance && return true
        # Healthy correction-merit steps reuse the Jacobian; age-based refresh only.
        refresh=correction_enabled ? age>20 : (age>=5 || contraction>.7)
        step=try
            if refresh
                J=_counterflow_jacobian!(f,u,w,r;previous,dt)
                correction_enabled && _counterflow_residual_scales!(residual_scales,J,u,bandwidth)
                _,pivots=LinearAlgebra.LAPACK.gbtrf!(bandwidth,bandwidth,length(u),J)
                age=0
            end
            correction=-vec(copy(r))
            LinearAlgebra.LAPACK.gbtrs!('N',bandwidth,bandwidth,length(u),w.band,pivots,correction)
            reshape(correction,size(u))
        catch e
            e isa SingularException || rethrow()
            return false
        end
        all(isfinite,step) || return false
        alpha=1.
        within_bounds=true
        for j in axes(u,2),k in axes(u,1)
            # Cantera permits -1e-7 trace species during damped Newton steps.
            # Enforcing strict positivity stalls initially absent hydrocarbon radicals.
            low=k==1 ? .2 : k<=n+1 ? -1e-7 : -1e5
            high=k==1 ? 6. : k<=n+1 ? 1.00001 : 1e5
            if k==size(u,1)-1
                # Cantera's physical Lambda bounds, converted to stored Lambda/Lscale.
                low=-1e20/_counterflow_Lscale
                high=1e20/_counterflow_Lscale
            end
            if !isnothing(f.control_points)
                k==n+2 && (low=0.)
                k==size(u,1)-2 && (low=-1e-5/_counterflow_Vscale)
            end
            within_bounds &= low<=u[k,j]<=high
            if step[k,j]<0
                alpha=min(alpha,.99*(u[k,j]-low)/(-step[k,j]))
            elseif step[k,j]>0
                alpha=min(alpha,.99*(high-u[k,j])/step[k,j])
            end
        end
        step_merit=correction_enabled ? _flame_correction_norm(step,weights) : 0.
        # A fresh small correction and componentwise backward error can certify
        # a root when the absolute residual is limited by state resolution.
        if correction_enabled && refresh && within_bounds && step_merit<1 &&
                _counterflow_constraints_converged(f,r,tolerance) &&
                all(i->abs(r[i])<=tolerance*residual_scales[i],eachindex(r))
            return true
        end
        accepted=false
        trial_merit=Inf
        for backtrack in 1:24
            @. trial=u+alpha*step
            counterflow_residual!(rt,f,trial,w;previous,dt)
            contracts=false
            if all(isfinite,rt)
                if correction_enabled
                    copyto!(trial_correction,vec(rt))
                    LinearAlgebra.LAPACK.gbtrs!('N',bandwidth,bandwidth,length(u),w.band,pivots,trial_correction)
                    trial_merit=_flame_correction_norm(reshape(trial_correction,size(u)),weights)
                    contracts=trial_merit<1 || trial_merit<step_merit
                else
                    contracts=norm(rt)<norm(r)*(1-1e-4*alpha)
                end
            end
            if contracts
                u .= trial; accepted=true; break
            end
            alpha*=.5
        end
        loglevel>1 && println("Counterflow Newton ",iteration," residual=",residualnorm," damping=",alpha)
        if !accepted || alpha<=1e-10
            refresh && return false
            age=correction_enabled ? 21 : 5; continue
        end
        contraction=norm(rt)/norm(r)
        # Keep the old residual through the contraction test, then reuse the
        # accepted trial's residual and full properties at the identical state.
        copyto!(r, rt)
        residual_valid = true
        age=correction_enabled && trial_merit<1 ? 21 : age+1
    end
    return false
end

function _counterflow_steady!(f;max_time_steps=800,timestep=Ref(1e-6),loglevel=0)
    w=CounterflowWorkspace(f)
    _counterflow_newton!(f,w;loglevel) && return true
    dt=timestep[]
    for step in 1:max_time_steps
        previous=copy(f.state)
        if _counterflow_newton!(f,w;previous,dt,maxiters=25)
            dt=min(dt*1.5,.01); timestep[]=dt
            if step%10==0
                loglevel>0 && println("Counterflow transient ",step," dt=",dt)
                _counterflow_newton!(f,w;loglevel) && return true
            end
        else
            f.state .= previous; dt*=.25; timestep[]=dt
            dt>1e-12 || return false
        end
    end
    return false
end

function _counterflow_reset!(f,grid)
    initialized=CounterflowDiffusionFlame(f.gas;
        fuel=mole_fractions(f.gas,f.fuel_Y;basis=:mass),
        oxidizer=mole_fractions(f.gas,f.oxidizer_Y;basis=:mass),
        mdot_fuel=f.fuel_mass_flux,mdot_oxidizer=f.oxidizer_mass_flux,
        T_fuel=f.fuel_temperature,T_oxidizer=f.oxidizer_temperature,P=f.pressure,grid)
    f.grid=initialized.grid
    f.state=initialized.state
    f.anchor=1
    empty!(f.fixed_temperature)
    return f
end

"""
Solve the native axisymmetric counterflow flame. With `auto=true`,
recover failed or extinguished coarse-grid solves using the native prescribed
initial temperature profile, then progressively finer initial grids. Failure to
find a burning solution raises an error and leaves `converged=false`.
With `auto=false`, a converged extinguished state is returned for extinction
continuation. Two-point control always retains the current branch and disables
automatic reinitialization; call `set_two_point_control!` to choose new targets.
"""
function solve!(f::CounterflowDiffusionFlame;refine_grid=true,ratio=4.,slope=.2,curve=.3,prune=0.,
        max_points=1200,max_time_steps=800,loglevel=0,initial_time_step=1e-6,auto=true)
    isfinite(ratio) && ratio>1 && 0<slope<=1 && 0<curve<=1 || throw(ArgumentError("invalid refinement criteria"))
    isfinite(prune) && prune<=min(slope,curve) || throw(ArgumentError("invalid pruning criterion"))
    isfinite(initial_time_step) && initial_time_step>0 || throw(ArgumentError("positive finite initial time step required"))
    if !isnothing(f.control_points)
        isempty(f.fixed_temperature) || throw(ArgumentError("two-point control conflicts with fixed_temperature"))
        zL,_,zR,_=f.control_points
        f.grid[1]<zL<zR<f.grid[end] && zL in f.grid && zR in f.grid ||
            throw(ArgumentError("two distinct interior control coordinates must remain on the grid"))
        auto=false
    end
    f.converged=false
    timestep=Ref(Float64(initial_time_step))
    initial_N=length(f.grid)
    attempts=auto && refine_grid ? unique(vcat(initial_N,[N for N in (12,24,48) if N>initial_N])) : [initial_N]
    burning=false
    for (attempt,N) in enumerate(attempts)
        N<=max_points || error("counterflow initialization exceeds max_points=$max_points")
        if attempt>1
            _counterflow_reset!(f,collect(range(f.grid[1],f.grid[end];length=N)))
            timestep[]=initial_time_step
        end
        loglevel>0 && println("Counterflow solve on ",length(f.grid)," points")
        success=_counterflow_steady!(f;max_time_steps,timestep,loglevel)
        if auto && (!success || extinct(f))
            loglevel>0 && println("Recovering with the native initial temperature profile")
            _counterflow_reset!(f,copy(f.grid))
            f.fixed_temperature=temperature(f)
            success=_counterflow_steady!(f;max_time_steps,timestep,loglevel)
            empty!(f.fixed_temperature)
            if success
                success=_counterflow_steady!(f;max_time_steps,timestep,loglevel)
            end
        end
        burning=success && (!auto || !extinct(f))
        burning && break
    end
    burning || error(auto ? "native counterflow solver did not find a burning solution" :
        "native counterflow continuation solver did not converge")
    for pass in 1:40
        if !refine_grid || !_refine_flame!(f;ratio,slope,curve,prune,max_points)
            if !isnothing(f.control_points)
                f.fuel_mass_flux=f.state[end,1]
                f.oxidizer_mass_flux=-f.state[end,end]
            end
            f.converged=true
            return f
        end
        loglevel>0 && println("Counterflow solve on ",length(f.grid)," points")
        _counterflow_steady!(f;max_time_steps,timestep,loglevel) ||
            error("native counterflow solver failed on $(length(f.grid)) points")
        auto && extinct(f) && error("counterflow flame extinguished during refinement")
    end
    error("counterflow refinement did not finish")
end

export CounterflowDiffusionFlame, CounterflowWorkspace, counterflow_residual!, spread_rate, pressure_curvature, extinct
export set_radiation!, radiation_source, radiative_heat_loss

export set_two_point_control!, disable_two_point_control!
