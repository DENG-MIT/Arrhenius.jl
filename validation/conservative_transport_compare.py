"""Independent CT4 face-property audit and matched-domain convergence estimates."""
import argparse,json
from pathlib import Path
import cantera as ct
import numpy as np
p=argparse.ArgumentParser(description=__doc__)
p.add_argument("mechanism",type=Path);p.add_argument("native",type=Path);p.add_argument("reference",type=Path)
p.add_argument("output",type=Path);a=p.parse_args()
g=ct.Solution(str(a.mechanism));mw=g.molecular_weights
audits={};convergence={}
for path in sorted(a.native.glob("*.npz")):
    if path.name.startswith("failed-"):continue
    case,rest=path.stem.split("-",1);mode,level=rest.rsplit("-",1)
    v=np.load(path);Y=v["Y"];T=v["T"];z=v["grid"];m=v["state"][-1];P=float(v["P"][0]);n,N=Y.shape
    multi=mode.startswith("multi");mass=mode.startswith("mass");soret=mode.endswith("soret")
    g.transport_model="multicomponent" if multi else "mixture-averaged"
    X=np.empty_like(Y);cp=np.empty_like(Y);h=np.empty_like(Y)
    for j in range(N):
        g.TPY=T[j],P,Y[:,j];X[:,j]=g.X;cp[:,j]=g.partial_molar_cp;h[:,j]=g.partial_molar_enthalpies
    F=np.empty_like(Y);FH=np.empty(N);Jall=np.empty((n,N-1))
    for j in range(N-1):
        dz=z[j+1]-z[j];ymid=.5*(Y[:,j]+Y[:,j+1]);tmid=.5*(T[j]+T[j+1])
        g.TPY=tmid,P,ymid;lam=g.thermal_conductivity
        cpface=np.sum(.25*(Y[:,j]+Y[:,j+1])*(cp[:,j]+cp[:,j+1])/mw)
        D=g.mix_diff_coeffs_mass if mass else g.mix_diff_coeffs
        Pe=abs(.5*(m[j]+m[j+1]))*dz/min(np.min(g.density*D),lam/cpface)
        c=1-Pe/6+Pe**3/360 if Pe<1e-4 else 1+2/Pe-1/np.tanh(Pe/2)
        if multi:J=g.density*mw/g.mean_molecular_weight**2*(g.multi_diff_coeffs@(mw*(X[:,j+1]-X[:,j])))/dz
        else:
            J=-g.density*D*(Y[:,j+1]-Y[:,j])/dz if mass else -g.density*mw/g.mean_molecular_weight*D*(X[:,j+1]-X[:,j])/dz
            J-=g.Y*sum(J)
        if soret:J-=g.thermal_diff_coeffs*(T[j+1]-T[j])/(tmid*dz)
        Jall[:,j]=J
        F[:,j]=.5*(m[j]+m[j+1])*((1-c/2)*Y[:,j]+c/2*Y[:,j+1])+J
        HL=np.sum(Y[:,j]*h[:,j]/mw);HR=np.sum(Y[:,j+1]*h[:,j+1]/mw)
        FH[j]=.5*(m[j]+m[j+1])*((1-c/2)*HL+c/2*HR)+np.sum(J*.5*(h[:,j]+h[:,j+1])/mw)-lam*(T[j+1]-T[j])/dz
    F[:,-1]=m[-1]*Y[:,-1];FH[-1]=m[-1]*np.sum(Y[:,-1]*h[:,-1]/mw)
    E=v["element_matrix"]
    d=dict(points=N,speed=float(v["velocity"][0]),Tmax=float(max(T)),
        element_outlet_error=float(np.max(np.abs(E@(Y[:,-1]-v["inlet_Y"])))),
        element_whole_domain_flux_error_over_mdot=float(np.max(np.abs(E@F-m[0]*(E@v["inlet_Y"])[:,None]))/m[0]),
        mass_diffusive_flux_error_over_mdot=float(np.max(np.abs(np.sum(Jall,axis=0)))/m[0]),
        native_CT_diffusive_flux_max_abs=float(np.max(np.abs(Jall-v["diffusive_flux"]))),
        native_CT_species_flux_max_abs=float(np.max(np.abs(F-v["species_flux"]))),
        native_CT_enthalpy_flux_max_abs=float(np.max(np.abs(FH-v["enthalpy_flux"]))))
    if case!="fixed":d["enthalpy_relative_flux_range"]=float(np.ptp(FH)/max(1.,np.max(np.abs(FH))))
    audits[path.stem]=d
for case in ["free","burner","fixed"]:
    for mode in ["mole","mass","mole-soret","mass-soret","multi","multi-soret"]:
        prefix=f"{case}-{mode}-"
        native=sorted(a.native.glob(prefix+"[0-9].npz"));reference=sorted(a.reference.glob(prefix+"[0-9].npz"))
        if len(native)<3 or len(reference)<3:continue
        def values(paths):
            out=[]
            for path in paths:
                v=np.load(path);out.append(dict(points=len(v["grid"]),domain=[float(v["grid"][0]),float(v["grid"][-1])],
                    speed=float(v["velocity"][0]),Tmax=float(max(v["T"]))))
            return out
        n,r=values(native),values(reference)
        assert n[-1]["domain"]==r[-1]["domain"],(case,mode,n[-1]["domain"],r[-1]["domain"])
        d=dict(native=n,reference=r)
        for key in ["speed","Tmax"]:
            # Estimates are labeled with their assumed order; sequences remain
            # visible so an unconverged estimate cannot masquerade as a limit.
            nl=(4*n[-1][key]-n[-2][key])/3;rl=2*r[-1][key]-r[-2][key]
            d[key+"_native_p2_estimate"]=nl;d[key+"_reference_p1_estimate"]=rl
            d[key+"_estimated_limit_difference"]=nl-rl
            d[key+"_native_extrapolation_change"]=nl-(4*n[-2][key]-n[-3][key])/3
            d[key+"_reference_extrapolation_change"]=rl-(2*r[-2][key]-r[-3][key])
        def profile(path):
            v=np.load(path);z=v["grid"].copy()
            # Free flames have an arbitrary spatial phase. Compare reaction
            # profiles after aligning the independently solved 700 K crossing.
            if case=="free":z-=np.interp(700.,v["T"],z)
            return z,v["T"],v["Y"]
        profiles=[profile(path) for path in [native[-2],native[-1],reference[-2],reference[-1]]]
        lo=max(v[0][0] for v in profiles);hi=min(v[0][-1] for v in profiles)
        q=np.unique(np.concatenate([v[0] for v in profiles]));q=q[(q>=lo)&(q<=hi)]
        interpolated=[(np.interp(q,z,T),np.stack([np.interp(q,z,y) for y in Y])) for z,T,Y in profiles]
        na,nb,ra,rb=interpolated
        d["profiles_aligned_at_K"]=700. if case=="free" else None
        d["raw_finest_temperature_profile_max_abs_K"]=float(np.max(np.abs(nb[0]-rb[0])))
        d["raw_finest_species_profile_max_abs"]=float(np.max(np.abs(nb[1]-rb[1])))
        d["extrapolated_temperature_profile_max_abs_K"]=float(np.max(np.abs((4*nb[0]-na[0])/3-(2*rb[0]-ra[0]))))
        d["extrapolated_species_profile_max_abs"]=float(np.max(np.abs((4*nb[1]-na[1])/3-(2*rb[1]-ra[1]))))
        baseline=profile(native[0]);coarse_z,coarse_T,coarse_Y=baseline
        mask=(q>=coarse_z[0])&(q<=coarse_z[-1]);qcoarse=q[mask]
        d["source_native_temperature_profile_error_against_reference_limit_K"]=float(np.max(np.abs(
            np.interp(qcoarse,coarse_z,coarse_T)-(2*rb[0]-ra[0])[mask])))
        d["source_native_species_profile_error_against_reference_limit"]=float(np.max(np.abs(
            np.stack([np.interp(qcoarse,coarse_z,y) for y in coarse_Y])-(2*rb[1]-ra[1])[:,mask])))
        convergence[f"{case}-{mode}"]=d
result=dict(source_commit="726522be4e2a13454d8415b7ef799d621f665cf3",cantera_version=ct.__version__,
    scope="Conservation/property correctness and equal-domain grid convergence; no benchmark qualification.",audits=audits,convergence=convergence)
a.output.write_text(json.dumps(result,indent=2))
for key in ["element_outlet_error","element_whole_domain_flux_error_over_mdot","native_CT_diffusive_flux_max_abs","native_CT_enthalpy_flux_max_abs","enthalpy_relative_flux_range"]:
    print(key,max(d.get(key,0) for d in audits.values()))
for name,d in convergence.items():
    print(name,"estimated limit differences",d["speed_estimated_limit_difference"],d["Tmax_estimated_limit_difference"],
        "last extrapolation changes",d["Tmax_native_extrapolation_change"],d["Tmax_reference_extrapolation_change"])
