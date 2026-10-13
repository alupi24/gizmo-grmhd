#!/usr/bin/env python3
"""
Face closure with GIZMO's HYDRO_TENSOR_FACE_CORRECTIONS (non-NGBITER variant), replicated from
hydro/gradients.c:1644-1676 (moments accumulated at gradients.c:2078-2105) and hydro/compute_finitevol_faces.h:42-57.

Per particle i, with dp = x_i - x_j and weights W_ij = W(r_ij, h_i) over r < h_i:
  S1 = sum W dp,  S2 = sum W dp dp,  S3.B = sum W dp dp (dp.B)
  B_i = -S2^{-1} S1                         -> corrected weight W (1 + B.dp) has zero first moment
  M_i = S2 - S3.B   (as coded)              [the corrected second moment sum W (1 + B.dp) dp dp would be S2 + S3.B]
  A_ij (from i) = V_i W_ij (1 + B_i.dp) M_i^{-1} dp + V_j W_ji (1 + B_j.(-dp)) M_j^{-1} dp
Reports the closure |sum_j A_ij| / sum_j |A_ij| for: standard MFM faces, tensor faces as coded, tensor faces with the
sign of the S3 term flipped, and the split of the standard closure into the own-side term (-V_i E_i^{-1} S1_i) and the rest.
"""
import os, importlib.util
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('p0', os.path.join(here, 'phase0_projection.py'))
p0 = importlib.util.module_from_spec(spec); spec.loader.exec_module(p0)

def closure_rel(geo, F):
    """F: face vectors as GIZMO's Face_Area_Vec seen from I (dp = x_I - x_J); outward A = -F."""
    C = -(p0.sadd(geo.N, geo.I, F) - p0.sadd(geo.N, geo.J, F))
    Al = np.linalg.norm(F, axis=1)
    S = p0.sadd(geo.N, geo.I, Al) + p0.sadd(geo.N, geo.J, Al)
    return np.sqrt(np.mean((np.linalg.norm(C, axis=1)/S)**2))

def tensor_faces(geo, sign_S3):
    dp = geo.dp
    N = geo.N
    # moments from both sides of each pair (own dp: +dp for I, -dp for J)
    S1 = p0.sadd(N, geo.I, geo.wkI[:, None]*dp) + p0.sadd(N, geo.J, geo.wkJ[:, None]*(-dp))
    o2 = dp[:, :, None]*dp[:, None, :]
    S2 = p0.sadd(N, geo.I, geo.wkI[:, None, None]*o2) + p0.sadd(N, geo.J, geo.wkJ[:, None, None]*o2)
    S2inv = np.linalg.inv(S2)
    B = -np.einsum('nab,nb->na', S2inv, S1)
    o3I = o2*(np.einsum('pa,pa->p', dp, B[geo.I]))[:, None, None]          # dp dp (dp.B_I)
    o3J = o2*(np.einsum('pa,pa->p', -dp, B[geo.J]))[:, None, None]         # (-dp)(-dp)(-dp.B_J)
    S3B = p0.sadd(N, geo.I, geo.wkI[:, None, None]*o3I) + p0.sadd(N, geo.J, geo.wkJ[:, None, None]*o3J)
    M = S2 + sign_S3*S3B
    Minv = np.linalg.inv(M)
    VI, VJ = geo.V[geo.I], geo.V[geo.J]
    fI = VI*geo.wkI*(1 + np.einsum('pa,pa->p', B[geo.I], dp))
    fJ = VJ*geo.wkJ*(1 - np.einsum('pa,pa->p', B[geo.J], dp))
    F = fI[:, None]*np.einsum('pab,pb->pa', Minv[geo.I], dp) + fJ[:, None]*np.einsum('pab,pb->pa', Minv[geo.J], dp)
    # first-moment check: sum_j W(1+B.dp) dp_own must vanish for every particle
    m1 = p0.sadd(N, geo.I, (geo.wkI*(1 + np.einsum('pa,pa->p', B[geo.I], dp)))[:, None]*dp) + \
         p0.sadd(N, geo.J, (geo.wkJ*(1 + np.einsum('pa,pa->p', B[geo.J], -dp)))[:, None]*(-dp))
    m1rel = np.max(np.linalg.norm(m1, axis=1)/np.linalg.norm(S1, axis=1).clip(1e-300))
    fdot = (F*dp).sum(1)
    return F, m1rel, int((fdot < 0).sum())

if __name__ == '__main__':
    L = ["# Face closure: standard MFM faces versus HYDRO_TENSOR_FACE_CORRECTIONS", "",
         "closure = rms_i |sum_j A_ij| / sum_j |A_ij|. `own-side` = the part of the standard closure coming from particle i's own "
         "kernel moments, -V_i E_i^{-1} S1_i (the term the tensor correction is designed to remove).", "",
         "| set | N | standard | own-side part | rest | tensor (as coded, S2 - S3.B) | tensor (S2 + S3.B) | first-moment check | faces with A.dp<0 (coded / flipped) |",
         "|---|---|---|---|---|---|---|---|---|"]
    for dname in ('lattice', 'glass', 'jitter'):
        for n in (16, 24):
            x = {'lattice': lambda: p0.lattice(n), 'jitter': lambda: p0.jittered(n, 0.3), 'glass': lambda: p0.glass(n)}[dname]()
            geo = p0.Geometry(x, 32, verbose=False)
            F_std = -geo.A
            # own-side split of the standard faces: F_own = wkI V_I NV_I dp (I side) ; the J side of the same pair is J's own term
            F_ownI = (geo.wkI*geo.V[geo.I])[:, None]*np.einsum('pab,pb->pa', geo.NV[geo.I], geo.dp)
            F_ownJ = (geo.wkJ*geo.V[geo.J])[:, None]*np.einsum('pab,pb->pa', geo.NV[geo.J], geo.dp)
            Cown = -(p0.sadd(geo.N, geo.I, F_ownI) - p0.sadd(geo.N, geo.J, F_ownJ))
            Crest = -(p0.sadd(geo.N, geo.I, F_ownJ) - p0.sadd(geo.N, geo.J, F_ownI))
            Al = np.linalg.norm(F_std, axis=1); S = p0.sadd(geo.N, geo.I, Al) + p0.sadd(geo.N, geo.J, Al)
            own = np.sqrt(np.mean((np.linalg.norm(Cown, axis=1)/S)**2)); rest = np.sqrt(np.mean((np.linalg.norm(Crest, axis=1)/S)**2))
            F_c, m1c, nbad_c = tensor_faces(geo, -1.0)
            F_f, m1f, nbad_f = tensor_faces(geo, +1.0)
            L.append(f"| {dname} | {n}^3 | {closure_rel(geo, F_std):.3e} | {own:.3e} | {rest:.3e} | {closure_rel(geo, F_c):.3e} | "
                     f"{closure_rel(geo, F_f):.3e} | {max(m1c, m1f):.1e} | {nbad_c} / {nbad_f} |")
            print(L[-1], flush=True)
    outp = os.path.join(here, "RESULTS_tensor_faces_ngb32.md")
    open(outp, 'w').write("\n".join(L) + "\n")
    print(f"# report -> {outp}")
