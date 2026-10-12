# Local derivatives for the conservative premixed modified-Newton Jacobian.
# Transport coefficients remain frozen, as in the finite-difference Jacobian.
import ForwardDiff

struct _FlameDerivativeTag end
const _FlameDual = ForwardDiff.Dual{_FlameDerivativeTag,Float64,1}
_flame_dual(x,dx) = _FlameDual(x,ForwardDiff.Partials((dx,)))
_flame_partial(x::_FlameDual) = ForwardDiff.partials(x)[1]

struct _ConservativeFlameDerivatives
    kinetics::KineticsWorkspace{Float64}
    dual_kinetics::KineticsWorkspace{_FlameDual}
    concentration::Vector{Float64}
    entropy::Vector{Float64}
    dual_C::Vector{_FlameDual}
    dual_h::Vector{_FlameDual}
    dual_cp::Vector{_FlameDual}
    dual_s::Vector{_FlameDual}
    dual_source::Vector{_FlameDual}
    cp_temperature::Matrix{Float64}
    density_gradient::Matrix{Float64}
    enthalpy_gradient::Matrix{Float64}
    cp_gradient::Matrix{Float64}
    concentration_jacobian::Matrix{Float64}
    rate_gradient::Vector{Float64}
    source_direction::Vector{Float64}
    flux_left::Matrix{Float64}
    flux_right::Matrix{Float64}
    projected_left::Vector{Float64}
    projected_right::Vector{Float64}
    uncorrected_flux::Vector{Float64}
    ymid::Vector{Float64}
end

function _ConservativeFlameDerivatives(gas,N)
    n=gas.n_species;B=n+2
    return _ConservativeFlameDerivatives(KineticsWorkspace(gas.reaction),
        KineticsWorkspace(gas.reaction,_FlameDual),zeros(n),zeros(n),
        zeros(_FlameDual,n),zeros(_FlameDual,n),zeros(_FlameDual,n),
        zeros(_FlameDual,n),zeros(_FlameDual,n),zeros(n,N),zeros(B,N),
        zeros(B,N),zeros(B,N),zeros(n,n),zeros(n),zeros(n),zeros(B,B),
        zeros(B,B),zeros(n),zeros(n),zeros(n),zeros(n))
end

# Product differentiation does not divide by concentration, so first-order
# reactants retain their finite derivative at exactly zero concentration.
function _flame_product_gradient!(gradient,A,reaction_index,C,scale)
    rows,orders=rowvals(A),nonzeros(A)
    @inbounds for slot in nzrange(A,reaction_index)
        order=orders[slot];k=rows[slot]
        order==0 && continue
        derivative=order==1 ? 1.0 : order*_concentration_power(C[k],order-1)
        for other in nzrange(A,reaction_index)
            other==slot && continue
            derivative *= _concentration_power(C[rows[other]],orders[other])
        end
        gradient[k] += scale*derivative
    end
end

function _flame_base_rate(reaction,i,T,h)
    value=reaction.Arrhenius_coeffs[i,1]*exp(reaction.Arrhenius_coeffs[i,2]*log(T)-
        reaction.Arrhenius_coeffs[i,3]*(4184.0/R)/T)
    for (j,index) in enumerate(reaction.blowers_masel.reaction_indices)
        index==i || continue
        p=reaction.blowers_masel.coefficients
        barrier=_blowers_masel_barrier(p[j,3],p[j,4],dot(@view(reaction.vk[:,i]),h))
        value=p[j,1]*exp(p[j,2]*log(T)-barrier/(R*T))
    end
    return value
end

function _flame_falloff_derivative(reaction,j,i,T,collider,high)
    # These are the RHS's actual zero/negative-pressure branches.
    collider<=0 && return 0.0
    low=reaction.Arrhenius_0[j,1]*exp(reaction.Arrhenius_0[j,2]*log(T)-
        reaction.Arrhenius_0[j,3]*(4184.0/R)/T)
    reduced=low*_flame_dual(collider,1.0)/high
    reduced<=0 && return 0.0
    value=high*reduced/(1+reduced)
    index=reaction.index_falloff_Troe[j]
    if index>0
        p=reaction.Troe_
        log_center=log10((1-p[index,1])*exp(-T/p[index,4])+
            p[index,1]*exp(-T/p[index,2])+exp(-p[index,3]/T))
        c=-.4-.67*log_center;N=.75-1.27*log_center
        f=(log10(reduced)+c)/(N-.14*(log10(reduced)+c))
        value *= exp(log(10.0)*log_center/(1+f*f))
    end
    return _flame_partial(value)
end

function _flame_concentration_jacobian!(d,reaction,T,C,h,s,cache)
    # Pressure-interpolated rates retain the generic finite-difference path.
    isempty(reaction.plog.reaction_indices) || return false
    wdot!(d.source_direction,reaction,T,C,s,h,d.kinetics;
        temperature_cache=cache,get_rate_constants=true)
    J=d.concentration_jacobian;fill!(J,0)
    kf,kr=d.kinetics.kf,d.kinetics.kr
    vr,vv=rowvals(reaction.vk),nonzeros(reaction.vk)
    er,ev=rowvals(reaction.efficiencies_coeffs),nonzeros(reaction.efficiencies_coeffs)
    @inbounds for i in 1:reaction.n_reactions
        fill!(d.rate_gradient,0)
        _flame_product_gradient!(d.rate_gradient,reaction.reactant_orders,i,C,kf[i])
        if reaction.is_reversible[i]
            _flame_product_gradient!(d.rate_gradient,reaction.product_stoich_coeffs,i,C,-kr[i])
        end
        pf=1.0;pr=reaction.is_reversible[i] ? 1.0 : 0.0
        for slot in nzrange(reaction.reactant_orders,i)
            pf *= _concentration_power(C[rowvals(reaction.reactant_orders)[slot]],nonzeros(reaction.reactant_orders)[slot])
        end
        if reaction.is_reversible[i]
            for slot in nzrange(reaction.product_stoich_coeffs,i)
                pr *= _concentration_power(C[rowvals(reaction.product_stoich_coeffs)[slot]],nonzeros(reaction.product_stoich_coeffs)[slot])
            end
        end
        multiplier=pf-(reaction.is_reversible[i] ? pr/d.kinetics.equilibrium_constants[i] : 0.0)
        collider_derivative=0.0
        if i in reaction.index_three_body
            collider_derivative=_flame_base_rate(reaction,i,T,h)
        else
            index=findfirst(==(i),reaction.index_falloff)
            if index!==nothing
                collider=0.0
                for slot in nzrange(reaction.efficiencies_coeffs,i)
                    collider += ev[slot]*C[er[slot]]
                end
                collider_derivative=_flame_falloff_derivative(reaction,index,i,T,collider,
                    _flame_base_rate(reaction,i,T,h))
            end
        end
        if collider_derivative!=0
            for slot in nzrange(reaction.efficiencies_coeffs,i)
                d.rate_gradient[er[slot]] += multiplier*collider_derivative*ev[slot]
            end
        end
        for k in eachindex(C)
            value=d.rate_gradient[k]
            for slot in nzrange(reaction.vk,i)
                J[vr[slot],k] += vv[slot]*value
            end
        end
    end
    return all(isfinite,J)
end

@inline function _flame_band_add!(band,B,j,k,jj,kk,value)
    row=(j-1)*B+k;column=(jj-1)*B+kk
    @inbounds band[4B-1+row-column,column] += value
end

function _flame_node_derivatives!(d,f,u,w,r,previous,dt,previous_enthalpy)
    n,N=f.gas.n_species,length(f.grid);B=n+2;gas=f.gas;MW=gas.MW;c=w.conservative
    fill!(d.density_gradient,0);fill!(d.enthalpy_gradient,0);fill!(d.cp_gradient,0)
    @inbounds for j in 1:N
        T=1000u[1,j];Td=_flame_dual(T,1000.0)
        sumY=0.0;inverseMW=0.0
        for k in 1:n
            positive=max(u[k+1,j],0.0)
            sumY+=positive;inverseMW+=positive/MW[k]
            d.concentration[k]=f.pressure/(R*T)*w.X[k,j]
            d.dual_C[k]=_flame_dual(d.concentration[k],-1000d.concentration[k]/T)
        end
        cal_cp_R!(d.dual_cp,gas,Td,f.pressure,@view(w.X[:,j]))
        cpT=0.0
        for k in 1:n
            d.cp_temperature[k,j]=R*_flame_partial(d.dual_cp[k])
            cpT+=u[k+1,j]*d.cp_temperature[k,j]/MW[k]
            d.enthalpy_gradient[k+1,j]=w.h[k,j]/MW[k]
            d.cp_gradient[k+1,j]=w.cp[k,j]/MW[k]
            # A forward derivative at zero agrees with the original positive
            # finite-difference perturbation. Negative trial species stay clipped.
            if u[k+1,j]>=0
                d.density_gradient[k+1,j]=w.rho[j]*(1/sumY-1/(MW[k]*inverseMW))
            end
        end
        d.density_gradient[1,j]=-1000w.rho[j]/T
        d.enthalpy_gradient[1,j]=1000c.cp[j];d.cp_gradient[1,j]=cpT
        j==1 && continue
        cal_h_RT!(d.dual_h,gas,Td,f.pressure,@view(w.X[:,j]))
        cal_s0_R!(d.dual_s,gas,Td,f.pressure,@view(w.X[:,j]))
        for k in 1:n
            d.dual_h[k] *= R*Td;d.dual_s[k] *= R
            d.entropy[k]=ForwardDiff.value(d.dual_s[k])
        end
        wdot!(d.dual_source,gas.reaction,Td,d.dual_C,d.dual_s,d.dual_h,d.dual_kinetics)
        _flame_concentration_jacobian!(d,gas.reaction,T,d.concentration,
            @view(w.h[:,j]),d.entropy,w.rate_caches[j]) || return false
        mul!(d.source_direction,d.concentration_jacobian,@view(w.X[:,j]))
        cell=.5*(f.grid[min(j+1,N)]-f.grid[j-1]);factor=_flame_timescale/w.rho[j]
        energy_flux=-factor/(1000c.cp[j])*(c.enthalpy_flux[j]-c.enthalpy_flux[j-1])/cell
        energy=energy_flux
        if previous!==nothing
            energy-=_flame_timescale/(1000c.cp[j]*dt)*(c.enthalpy[j]-previous_enthalpy[j])
        end
        for l in 1:B
            density_ratio=d.density_gradient[l,j]/w.rho[j]
            for k in 1:n
                chemistry=l==1 ? _flame_partial(d.dual_source[k]) :
                    l==B || u[l,j]<0 ? 0.0 : f.pressure/(R*T*MW[l-1]*inverseMW)*
                        (d.concentration_jacobian[k,l-1]-d.source_direction[k])
                balance=r[k+1,j]
                if previous!==nothing
                    balance+=_flame_timescale/dt*(u[k+1,j]-previous[k+1,j])
                end
                value=factor*MW[k]*chemistry-balance*density_ratio
                previous!==nothing && l==k+1 && (value-=_flame_timescale/dt)
                _flame_band_add!(w.band,B,j,k+1,j,l,value)
            end
            value=-energy_flux*density_ratio-energy*d.cp_gradient[l,j]/c.cp[j]
            if previous!==nothing
                value-=_flame_timescale/(1000c.cp[j]*dt)*d.enthalpy_gradient[l,j]
            end
            _flame_band_add!(w.band,B,j,1,j,l,value)
        end
    end
    return true
end

@inline _flame_centering_derivative(Pe) = Pe<1e-4 ? -1/6+Pe^2/120 :
    -2/Pe^2+(Pe>700 ? 0.0 : .5/sinh(Pe/2)^2)

function _flame_face_derivatives!(d,f,u,w,j)
    n=f.gas.n_species;B=n+2;MW=f.gas.MW;c=w.conservative
    L,Rj=d.flux_left,d.flux_right;fill!(L,0);fill!(Rj,0)
    dz=f.grid[j+1]-f.grid[j];Tl=1000u[1,j];Tr=1000u[1,j+1];Tm=.5*(Tl+Tr)
    ysum=0.0;inverse_left=0.0;inverse_right=0.0
    @inbounds for k in 1:n
        d.ymid[k]=max(.5*(u[k+1,j]+u[k+1,j+1]),0.0);ysum+=d.ymid[k]
        inverse_left+=max(u[k+1,j],0.0)/MW[k]
        inverse_right+=max(u[k+1,j+1],0.0)/MW[k]
    end
    d.ymid ./= ysum
    if f.transport_model==:multicomponent
        @inbounds for k in 1:n
            pl=0.0;pr=0.0
            for l in 1:n
                pl+=w.multi_prefactor[k,l,j]*w.X[l,j]
                pr+=w.multi_prefactor[k,l,j]*w.X[l,j+1]
            end
            d.projected_left[k]=pl;d.projected_right[k]=pr
        end
        @inbounds for l in 1:n,k in 1:n
            L[k+1,l+1]=u[l+1,j]>=0 ? -(w.multi_prefactor[k,l,j]-d.projected_left[k])/(dz*MW[l]*inverse_left) : 0.0
            Rj[k+1,l+1]=u[l+1,j+1]>=0 ? (w.multi_prefactor[k,l,j]-d.projected_right[k])/(dz*MW[l]*inverse_right) : 0.0
        end
    else
        fluxsum=0.0
        @inbounds for k in 1:n
            gradient=f.flux_gradient_basis==:mass ? u[k+1,j+1]-u[k+1,j] : w.X[k,j+1]-w.X[k,j]
            d.uncorrected_flux[k]=-w.diffusion_prefactor[k,j]*gradient/dz
            fluxsum+=d.uncorrected_flux[k]
        end
        @inbounds for l in 1:n
            sl=0.0;sr=0.0
            for k in 1:n
                if f.flux_gradient_basis==:mass
                    dl=k==l ? w.diffusion_prefactor[k,j]/dz : 0.0;dr=-dl
                else
                    dl=u[l+1,j]>=0 ? w.diffusion_prefactor[k,j]*((k==l)-w.X[k,j])/(dz*MW[l]*inverse_left) : 0.0
                    dr=u[l+1,j+1]>=0 ? -w.diffusion_prefactor[k,j]*((k==l)-w.X[k,j+1])/(dz*MW[l]*inverse_right) : 0.0
                end
                L[k+1,l+1]=dl;Rj[k+1,l+1]=dr;sl+=dl;sr+=dr
            end
            midpoint_derivative=u[l+1,j]+u[l+1,j+1]>=0 ? .5/ysum : 0.0
            for k in 1:n
                correction=midpoint_derivative*((k==l)-d.ymid[k])*fluxsum
                L[k+1,l+1]-=d.ymid[k]*sl+correction
                Rj[k+1,l+1]-=d.ymid[k]*sr+correction
            end
        end
    end
    if f.soret_enabled
        @inbounds for k in 1:n
            L[k+1,1]=w.thermal_diffusion[k,j]*1000Tr/(dz*Tm^2)
            Rj[k+1,1]=-w.thermal_diffusion[k,j]*1000Tl/(dz*Tm^2)
        end
    end
    cpface=0.0;cpTl=0.0;cpTr=0.0;density_diffusion=Inf
    @inbounds for k in 1:n
        ysumk=u[k+1,j]+u[k+1,j+1]
        cpface+=.25*ysumk*(w.cp[k,j]+w.cp[k,j+1])/MW[k]
        cpTl+=.25*ysumk*d.cp_temperature[k,j]/MW[k]
        cpTr+=.25*ysumk*d.cp_temperature[k,j+1]/MW[k]
        density_diffusion=min(density_diffusion,c.density_diffusion[k,j])
    end
    thermal=w.conductivity[j]/cpface<density_diffusion
    thermal && (density_diffusion=w.conductivity[j]/cpface)
    mdot=.5*(u[B,j]+u[B,j+1]);Pe=abs(mdot)*dz/density_diffusion
    rw=.5*_flame_face_centering(Pe);mdot<0 && (rw=1-rw)
    sign_m=mdot<0 ? -1.0 : 1.0
    drdPe=.5*sign_m*_flame_centering_derivative(Pe)
    @inbounds for l in 1:B
        cpl=l==1 ? cpTl : l==B ? 0.0 : .25*(w.cp[l-1,j]+w.cp[l-1,j+1])/MW[l-1]
        cpr=l==1 ? cpTr : l==B ? 0.0 : cpl
        dpl=thermal ? Pe*cpl/cpface : 0.0;dpr=thermal ? Pe*cpr/cpface : 0.0
        if l==B
            dpl+=.5*sign_m*dz/density_diffusion;dpr+=.5*sign_m*dz/density_diffusion
        end
        drl=drdPe*dpl;drr=drdPe*dpr
        hdiff=c.enthalpy[j+1]-c.enthalpy[j]
        hl=mdot*((1-rw)*d.enthalpy_gradient[l,j]+drl*hdiff)
        hr=mdot*(rw*d.enthalpy_gradient[l,j+1]+drr*hdiff)
        if l==B
            hmean=(1-rw)*c.enthalpy[j]+rw*c.enthalpy[j+1]
            hl+=.5*hmean;hr+=.5*hmean
        elseif l==1
            hl+=1000w.conductivity[j]/dz;hr-=1000w.conductivity[j]/dz
        end
        for k in 1:n
            hmean=.5*(w.h[k,j]+w.h[k,j+1])/MW[k]
            hl+=hmean*L[k+1,l];hr+=hmean*Rj[k+1,l]
            if l==1
                hl+=500w.cp[k,j]*w.flux[k,j]/MW[k]
                hr+=500w.cp[k,j+1]*w.flux[k,j]/MW[k]
            end
            difference=u[k+1,j+1]-u[k+1,j]
            L[k+1,l]+=mdot*(drl*difference+(l==k+1 ? 1-rw : 0.0))
            Rj[k+1,l]+=mdot*(drr*difference+(l==k+1 ? rw : 0.0))
            if l==B
                mean_y=(1-rw)*u[k+1,j]+rw*u[k+1,j+1]
                L[k+1,l]+=.5*mean_y;Rj[k+1,l]+=.5*mean_y
            end
        end
        L[1,l]=hl;Rj[1,l]=hr
    end
end

function _conservative_flame_jacobian!(f,u,w,r;previous=nothing,dt=Inf,previous_enthalpy=nothing)
    isempty(f.gas.reaction.plog.reaction_indices) || return false
    n,N=f.gas.n_species,length(f.grid);B=n+2;d=w.conservative.derivatives;c=w.conservative
    previous!==nothing && previous_enthalpy===nothing &&
        (previous_enthalpy=_flame_previous_enthalpy!(c,f,previous))
    fill!(w.band,0)
    _flame_node_derivatives!(d,f,u,w,r,previous,dt,previous_enthalpy) || return false
    @inbounds for face in 1:N-1
        _flame_face_derivatives!(d,f,u,w,face)
        for side in 0:1
            node=face+side;sign=side==0 ? -1.0 : 1.0
            cell=node==1 ? 1.0 : .5*(f.grid[min(node+1,N)]-f.grid[node-1])
            factor=node==1 ? sign : sign*_flame_timescale/(w.rho[node]*cell)
            for other in 0:1
                matrix=other==0 ? d.flux_left : d.flux_right
                for l in 1:B,k in 1:n
                    _flame_band_add!(w.band,B,node,k+1,face+other,l,factor*matrix[k+1,l])
                end
                if node>1
                    for l in 1:B
                        _flame_band_add!(w.band,B,node,1,face+other,l,
                            factor/(1000c.cp[node])*matrix[1,l])
                    end
                end
            end
        end
    end
    # Natural outlet face F=mY and FH=mh.
    factor=-_flame_timescale/(w.rho[N]*.5*(f.grid[N]-f.grid[N-1]))
    @inbounds for k in 1:n
        _flame_band_add!(w.band,B,N,k+1,N,k+1,factor*u[B,N])
        _flame_band_add!(w.band,B,N,k+1,N,B,factor*u[k+1,N])
        _flame_band_add!(w.band,B,1,k+1,1,B,f.inlet_Y[k])
    end
    @inbounds for l in 1:B
        value=u[B,N]*d.enthalpy_gradient[l,N]+(l==B ? c.enthalpy[N] : 0.0)
        _flame_band_add!(w.band,B,N,1,N,l,factor/(1000c.cp[N])*value)
    end
    # Replace algebraic rows after assembling physical balances.
    @inbounds for j in 1:N
        for jj in max(1,j-1):min(N,j+1),l in 1:B
            row=(j-1)*B+f.dependent_species+1;column=(jj-1)*B+l
            w.band[4B-1+row-column,column]=0
            if j==1 || f isa BurnerFlame && !isempty(f.imposed_temperature)
                row=(j-1)*B+1;w.band[4B-1+row-column,column]=0
            end
        end
        for k in 1:n
            _flame_band_add!(w.band,B,j,f.dependent_species+1,j,k+1,1.0)
        end
        if j==1 || f isa BurnerFlame && !isempty(f.imposed_temperature)
            _flame_band_add!(w.band,B,j,1,j,1,1.0)
        end
        if f isa BurnerFlame
            _flame_band_add!(w.band,B,j,B,j,B,1.0)
            j>1 && _flame_band_add!(w.band,B,j,B,j-1,B,-1.0)
        elseif j==f.anchor
            _flame_band_add!(w.band,B,j,B,j,1,1.0)
        elseif j<f.anchor
            _flame_band_add!(w.band,B,j,B,j+1,B,1.0)
            _flame_band_add!(w.band,B,j,B,j,B,-1.0)
        else
            _flame_band_add!(w.band,B,j,B,j,B,1.0)
            _flame_band_add!(w.band,B,j,B,j-1,B,-1.0)
        end
    end
    return all(isfinite,w.band)
end
