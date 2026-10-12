# Planar ionized-flame balances, following Cantera
# Flow1D.cpp, Boundary1D.cpp and IonFlow.cpp
# at 726522be4e2a13454d8415b7ef799d621f665cf3. The translated equations retain
# Cantera's license:
# Copyright (c) 2001-2009, California Institute of Technology. All rights reserved.
# Copyright (c) 2009 Sandia Corporation. Under the terms of Contract
# AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain
# rights in this software.
# Copyright (c) 2011-2026, Cantera Developers. All rights reserved.
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# - Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# - Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# - Neither the name of the California Institute of Technology, Sandia
#   Corporation nor the names of other contributors may be used to endorse or
#   promote products derived from this software without specific prior written
#   permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


function flame_residual!(residual,f::IonizedFlame,u,w::IonizedFlameWorkspace;
        previous=nothing,dt=Inf,update_transport=true,nodes=eachindex(f.grid),
        previous_enthalpy=nothing)
    n,N=f.gas.n_species,length(f.grid)
    size(residual)==size(u)==(n+3,N) || throw(DimensionMismatch("ionized residual/state size mismatch"))
    previous===nothing || size(previous)==size(u) || throw(DimensionMismatch("previous state size mismatch"))
    dt>0 || throw(ArgumentError("positive time step required"))
    _ion_flame_properties!(w,f,u;update_transport,nodes)
    c=w.conservative
    if c !== nothing
        _ion_conservative_fluxes!(c,f,u,w)
        if previous !== nothing && previous_enthalpy === nothing
            previous_enthalpy=_flame_previous_enthalpy!(c,f,previous)
        end
    end
    z,MW,q=f.grid,f.gas.MW,f.ion_data.charges
    mdotin=f.kind==:free ? w.rho[1]*u[end,1] : f.mass_flux
    @inbounds for j in 1:N
        mdot=w.rho[j]*u[end,j]
        if f.kind==:burner
            residual[end,j]=j==1 ? mdot-f.mass_flux : mdot-w.rho[j-1]*u[end,j-1]
        elseif j==N
            residual[end,j]=mdot-w.rho[j-1]*u[end,j-1]
        elseif j==f.anchor
            residual[end,j]=u[1,j]-f.fixed_temperature/1000
        elseif j<f.anchor
            residual[end,j]=mdot-w.rho[j+1]*u[end,j+1]
        else
            residual[end,j]=w.rho[j-1]*u[end,j-1]-mdot
        end
        if j==1
            residual[1,j]=u[1,j]-f.inlet_temperature/1000
            for k in 1:n
                residual[k+1,j]=c===nothing ?
                    mdotin*f.inlet_Y[k]-mdot*u[k+1,j]-w.flux[k,1] : mdotin*f.inlet_Y[k]-c.species_flux[k,1]
            end
            residual[w.excess[1]+1,j]=1-sum(@view(u[2:n+1,j]))
            if f.field_enabled
                for k in 1:n
                    if q[k]!=0
                        residual[k+1,j]=u[k+1,j]-u[k+1,j+1]
                        k==w.excess[1] || (residual[k+1,j]+=mdotin*f.inlet_Y[k])
                    end
                end
            end
        elseif j==N && c===nothing
            residual[1,j]=u[1,j]-u[1,j-1]
            for k in 1:n
                residual[k+1,j]=u[k+1,j]-u[k+1,j-1]
            end
            residual[w.excess[2]+1,j]=1-sum(@view(u[2:n+1,j]))
        elseif c !== nothing
            _ion_conservative_cell!(residual,f,u,w,j;previous,dt,previous_enthalpy)
        else
            left=z[j]-z[j-1];right=z[j+1]-z[j];cell=(left+right)/2
            up=u[end,j]>0 ? j : j+1
            dz=z[up]-z[up-1]
            cpmean,chemical,enthalpyflux=0.0,0.0,0.0
            for k in 1:n
                cpmean+=u[k+1,j]*w.cp[k,j]/MW[k]
                chemical+=w.h[k,j]*w.source[k,j]
                enthalpyflux+=(w.flux[k,j-1]+w.flux[k,j])/2 *
                    (w.h[k,up]-w.h[k,up-1])/(dz*MW[k])
                residual[k+1,j]=_flame_timescale/w.rho[j] *
                    (MW[k]*w.source[k,j]-(w.flux[k,j]-w.flux[k,j-1])/cell -
                     mdot*(u[k+1,up]-u[k+1,up-1])/dz)
            end
            conduction=1000*(w.conductivity[j]*(u[1,j+1]-u[1,j])/right -
                w.conductivity[j-1]*(u[1,j]-u[1,j-1])/left)/cell
            dTdz=1000*(u[1,up]-u[1,up-1])/dz
            residual[1,j]=_flame_timescale/(1000*w.rho[j]*cpmean) *
                (conduction-chemical-enthalpyflux-mdot*cpmean*dTdz)
            if previous!==nothing
                for k in 1:n+1
                    residual[k,j]-=_flame_timescale/dt*(u[k,j]-previous[k,j])
                end
            end
        end
        if !f.field_enabled || j==1
            residual[end-1,j]=u[end-1,j]
        else
            charge=0.0
            for k in 1:n
                charge+=q[k]*u[k+1,j]/MW[k]
            end
            residual[end-1,j]=u[end-1,j]-u[end-1,j-1] -
                (z[j]-z[j-1])*_ION_FLAME_FARADAY*w.rho[j]*charge/(1000*_ION_FLAME_EPS0)
        end
    end
    return residual
end
