# Compare fixed-T Cantera and Arrhenius Jacobians after mapping variables and rows.
# Usage: python fixed_t_matrix_cross.py cantera.npz arr_base.npz arr_h01.npz arr_h10.npz case_id
import sys
import numpy as np
from scipy.sparse import csc_matrix, coo_matrix, diags
from scipy.sparse.linalg import splu

ct=np.load(sys.argv[1]); ar=np.load(sys.argv[2]); arlo=np.load(sys.argv[3]); arhi=np.load(sys.argv[4]); name=sys.argv[5]; ctsteps=[np.load(p) for p in sys.argv[6:]]
shape=tuple(int(x) for x in ct["shape"])
Jct=csc_matrix((ct["J_data"],ct["J_indices"],ct["J_indptr"]),shape=shape)
Jctsteps=[csc_matrix((d["J_data"],d["J_indices"],d["J_indptr"]),shape=shape) for d in ctsteps]

def arr_matrix(d):
    bw=int(d["bw"][0]); bands=d["bands_arr"]; m=bands.shape[1]
    rows=[]; cols=[]; vals=[]
    for col in range(m):
        lo=max(0,col-bw); hi=min(m,col+bw+1)
        for row in range(lo,hi):
            value=bands[bw+row-col,col]
            if value != 0.0:
                rows.append(row); cols.append(col); vals.append(value)
    return coo_matrix((vals,(rows,cols)),shape=(m,m)).tocsc()

Jarr=arr_matrix(ar); Jlo=arr_matrix(arlo); Jhi=arr_matrix(arhi)
N=len(ct["grid"]); Y=np.asarray(ct["Y"]); T=np.asarray(ct["T"]); U=np.asarray(ct["velocity"])
MW=np.asarray(ct["molecular_weights"]); n=len(MW); B=n+2; m=N*B
assert shape==(m,m) and Jarr.shape==(m,m) and Y.shape==(N,n)
meanmw=1.0/np.sum(Y/MW[None,:],axis=1)
rho=np.asarray(ar["rho"]); assert len(rho)==N
# D = d(T/1000, Y, rho*U) / d(T, Y, U), point-major.
rows=[]; cols=[]; vals=[]
for j in range(N):
    b=j*B
    rows.append(b); cols.append(b); vals.append(1e-3)
    for k in range(n):
        rows.append(b+1+k); cols.append(b+1+k); vals.append(1.0)
    rows.append(b+B-1); cols.append(b); vals.append(-rho[j]*U[j]/T[j])
    for k in range(n):
        rows.append(b+B-1); cols.append(b+1+k); vals.append(-rho[j]*U[j]*meanmw[j]/MW[k])
    rows.append(b+B-1); cols.append(b+B-1); vals.append(rho[j])
D=coo_matrix((vals,(rows,cols)),shape=(m,m)).tocsc()
# Native residual rows relative to Cantera: fixed-T T rows x1e-3; interior
# species rows x1e-4; boundary species flux/gradient rows x1; closure rows x-1;
# burner mass-flow rows x1. This is r_arr = S * r_ct when equations match.
S=np.ones(m)
labels=["T",*list(ct["species_names"]),"velocity"]
for j in range(N):
    b=j*B; S[b]=1e-3
    if 0<j<N-1:
        S[b+1:b+1+n]=1e-4
    else:
        excess=int(np.argmax(Y[j]))
        native_excess=int(ar["dependent_species"][0])-1
        assert excess==native_excess, f"closure species mismatch at boundary {j}: CT={excess}, Arrhenius={native_excess}"
        S[b+1+excess]=-1.0
# The inlet row forms differ by the inlet continuity row: native mdot*(Yin-Y)-F
# equals Cantera mdot_in*Yin-rhoU*Y-F plus Yin*(rhoU-mdot_in).
Yin=np.asarray(ct["inlet_Y"]); native_excess=int(ar["dependent_species"][0])-1
assert len(Yin)==n and np.allclose(Yin,np.asarray(ar["inlet_Y"]),rtol=1e-12,atol=1e-15)
for j in (0,N-1):
    assert int(np.argmax(Y[j]))==native_excess, f"closure species mismatch at boundary {j}"
flow0=B-1
mix_rows=[]; mix_cols=[]; mix_vals=[]; inv_vals=[]
for k in range(n):
    if k != native_excess:
        mix_rows.append(1+k); mix_cols.append(flow0); mix_vals.append(Yin[k]); inv_vals.append(-Yin[k])
M=diags(S,format="csc")+coo_matrix((mix_vals,(mix_rows,mix_cols)),shape=(m,m)).tocsc()
Minv=diags(1.0/S,format="csc")+coo_matrix((inv_vals,(mix_rows,mix_cols)),shape=(m,m)).tocsc()
JA=Minv @ Jarr @ D
JLo=Minv @ Jlo @ D
JHi=Minv @ Jhi @ D
rct=np.asarray(ct["r_ct"]); rarr_ct=Minv @ np.asarray(ar["r_arr"])
# Residual/Jacobian 2x2, all steps expressed in Cantera physical unknowns.
def norm(dx): return float(np.linalg.norm(dx/ct["weights"])/np.sqrt(int(ct["dofs_full"][0])))
def solve_norm(J,rhs): return norm(splu(J.tocsc()).solve(-rhs))
print(f"case={name} points={N} dofs={m}")
print("row_map: T=1e-3; interior_species=1e-4; inlet/outlet_species_flux_or_gradient=1; closure=-1; flow=1; inlet_species_native=Cantera_species+Yin*Cantera_inlet_flow")
print("variable_map: T_native=T_K/1000; Y_native=Y; mdot_native=rho*U with dmdot/dT=-rho*U/T, dmdot/dYk=-rho*U*MWmix/MWk, dmdot/dU=rho")
rmap=M @ rct
rarr=np.asarray(ar["r_arr"])
resid_error=rarr-rmap
rmax=int(np.argmax(np.abs(resid_error)))
print(f"residual_map_rel_inf={np.max(np.abs(resid_error))/max(np.max(np.abs(rarr)),np.max(np.abs(rmap)),1e-300):.9g} max_at(point={rmax//B+1} component={labels[rmax%B]}) arr={rarr[rmax]:.9g} mapped_ct={rmap[rmax]:.9g}")
print(f"2x2 weighted_step_rms: Arr/Arr={solve_norm(JA,rarr_ct):.12g} CT/Arr={solve_norm(JA,rct):.12g} CT/CT={solve_norm(Jct,rct):.12g} Arr/CT={solve_norm(Jct,rarr_ct):.12g}")
# Causal intervention: replace only outlet rows versus previous-point
# species columns with the mapped Arrhenius finite-difference block.
from scipy.sparse import lil_matrix
jpatch=Jct.tolil(copy=True)
outlet_rows=np.arange(m-B,m)
previous_species=np.arange((N-2)*B+1,(N-1)*B-1)
jpatch[np.ix_(outlet_rows,previous_species)]=JA[np.ix_(outlet_rows,previous_species)]
Jpatch=jpatch.tocsc()
print(f"causal_outlet_previous_species_patch changed_entries={(Jpatch-Jct).nnz} CT_residual_step_after_patch={solve_norm(Jpatch,rct):.12g}")
# Separate intervention for the interior species-by-species block.
species_interior=np.concatenate([np.arange(j*B+1,j*B+1+n) for j in range(1,N-1)])
jspatch=Jct.tolil(copy=True)
jspatch[np.ix_(species_interior,species_interior)]=JA[np.ix_(species_interior,species_interior)]
Jspatch=jspatch.tocsc()
print(f"causal_interior_species_block_patch changed_entries={(Jspatch-Jct).nnz} CT_residual_step_after_patch={solve_norm(Jspatch,rct):.12g}")
# Split boundary intervention to isolate the inlet rows from outlet rows.
inlet_rows=np.arange(B)
outlet_rows_all=np.arange(m-B,m)
jipatch=Jct.tolil(copy=True)
jipatch[np.ix_(inlet_rows,species_interior)]=JA[np.ix_(inlet_rows,species_interior)]
print(f"causal_inlet_species_patch changed_entries={(jipatch.tocsc()-Jct).nnz} CT_residual_step_after_patch={solve_norm(jipatch.tocsc(),rct):.12g}")
jopatch=Jct.tolil(copy=True)
jopatch[np.ix_(outlet_rows_all,species_interior)]=JA[np.ix_(outlet_rows_all,species_interior)]
print(f"causal_outlet_all_species_patch changed_entries={(jopatch.tocsc()-Jct).nnz} CT_residual_step_after_patch={solve_norm(jopatch.tocsc(),rct):.12g}")

# Four disjoint row/column groups: interior T, interior species, interior flow,
# and all boundary rows/columns. A block is normalized by its largest Cantera entry.
interior=np.arange(B,m-B)
t_rows=interior[interior % B == 0]
flow_rows=interior[interior % B == B-1]
species_rows=interior[(interior % B != 0)&(interior % B != B-1)]
boundary=np.r_[np.arange(B),np.arange(m-B,m)]
groups=[("T_interior",t_rows),("species_interior",species_rows),("flow_interior",flow_rows),("boundary",boundary)]
noise=JLo-JA
noise2=JHi-JA
def coord(index):
    j,k=divmod(int(index),B)
    return f"point={j+1} component={labels[k]}"
for rn,ri in groups:
    for cn,ci in groups:
        # Intervention: replace exactly this mapped matrix block and solve
        # with the unchanged Cantera residual.
        jb=Jct.tolil(copy=True)
        jb[np.ix_(ri,ci)]=JA[np.ix_(ri,ci)]
        print(f"causal_block_patch {rn} x {cn} CT_residual_step_after_patch={solve_norm(jb.tocsc(),rct):.12g}")
        A=JA[ri,:][:,ci].toarray(); C=Jct[ri,:][:,ci].toarray()
        delta=np.abs(A-C); scale=max(float(np.max(np.abs(C))) if C.size else 0.0,1e-300)
        denom=max(float(np.max(np.abs(C))) if C.size else 0.0,float(np.max(np.abs(A))) if A.size else 0.0,1e-300)
        q=np.unravel_index(int(np.argmax(delta)),delta.shape) if delta.size else (0,0)
        global_row=int(ri[q[0]]) if len(ri) else 0; global_col=int(ci[q[1]]) if len(ci) else 0
        N0=np.abs(noise[ri,:][:,ci].toarray()); N1=np.abs(noise2[ri,:][:,ci].toarray())
        arr_noise=np.maximum(N0,N1); ct_noise=np.zeros_like(delta)
        for variant in Jctsteps:
            ct_noise=np.maximum(ct_noise,np.abs((variant-Jct)[ri,:][:,ci].toarray()))
        nmax=max(float(np.max(arr_noise)) if arr_noise.size else 0.0,float(np.max(ct_noise)) if ct_noise.size else 0.0)
        threshold=5*np.maximum(arr_noise,ct_noise)+1e-12*scale
        beyond=np.argwhere(delta>threshold)
        detail="none"
        if len(beyond):
            order=sorted(beyond.tolist(),key=lambda ij:delta[ij[0],ij[1]]/max(threshold[ij[0],ij[1]],1e-300),reverse=True)[:5]
            detail=";".join(f"{coord(ri[i])} <- {coord(ci[j])} diff={delta[i,j]:.4g} noise={max(arr_noise[i,j],ct_noise[i,j]):.4g} A={A[i,j]:.4g} C={C[i,j]:.4g}" for i,j in order)
        if rn == "boundary" and cn == "species_interior":
            orderabs=np.dstack(np.unravel_index(np.argsort(delta.ravel())[-10:],delta.shape))[0][::-1]
            print("boundary_species_largest=" + ";".join(f"{coord(ri[i])} <- {coord(ci[j])} A={A[i,j]:.8g} C={C[i,j]:.8g}" for i,j in orderabs))
        print(f"block {rn} x {cn}: max_rel_ct={float(np.max(delta))/scale:.9g} max_rel_symmetric={float(np.max(delta))/denom:.9g} max_abs={float(np.max(delta)):.9g} max_arr={A[q[0],q[1]]:.9g} max_ct={C[q[0],q[1]]:.9g} fd_step_sensitivity_rel={nmax/scale:.9g} arr_fd_step_sensitivity_rel={float(np.max(arr_noise))/scale if arr_noise.size else 0.0:.9g} ct_fd_step_sensitivity_rel={float(np.max(ct_noise))/scale if ct_noise.size else 0.0:.9g} max_at({coord(global_row)} <- {coord(global_col)}) beyond_5x_sensitivity={len(beyond)} [{detail}]")