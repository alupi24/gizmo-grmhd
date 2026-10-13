#!/usr/bin/env python3
"""
Follow-up to phase0_projection.py: how much of GIZMO's face-sum divergence is face non-closure?

  D_i = sum_j 1/2 (B'_i + B'_j) . A_ij  =  B_i . C_i  +  R_i,      C_i = sum_j A_ij (closure vector, zero for closed faces)

B_i . C_i is the divergence a uniform field equal to B_i would show (zeroth order, does not converge with resolution);
R_i is the part produced by the field variation and its reconstruction. Also reports, for the exact projection, how much
the face normal fluxes are changed (relative to rms|B|), split the same way by running it on D and on R alone.

Usage: python3 phase0_closure_split.py [--ns 16 24 32] [--dists lattice jitter glass] [--fields uniform tg3d abc grad_long bump]
"""
import argparse, os, importlib.util
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('p0', os.path.join(here, 'phase0_projection.py'))
p0 = importlib.util.module_from_spec(spec); spec.loader.exec_module(p0)

def rmsn(geo, D, Brms):
    return np.sqrt(np.mean((geo.dx*D/geo.V)**2))/Brms

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--ns', type=int, nargs='+', default=[16, 24, 32])
    ap.add_argument('--ngb', type=int, default=32)
    ap.add_argument('--dists', nargs='+', default=['lattice', 'jitter', 'glass'])
    ap.add_argument('--fields', nargs='+', default=['uniform', 'tg3d', 'abc', 'grad_long', 'bump'])
    a = ap.parse_args()
    L = ["# Face non-closure versus reconstruction in GIZMO's face-sum divergence", "",
         "`div(D)`: rms dx|D_i|/V_i/rms|B| of the full face-sum divergence (gradient std, limiters on). "
         "`div(B.C)`: same for the closure term B_i . sum_j A_ij. `div(R)`: same for the remainder R = D - B.C. "
         "`dface(exact)`: rms change of the face normal field B_face.n produced by the exact projection of D, relative to rms|B|; "
         "`dface(exact on R)`: the same if only R had to be removed (i.e. with closed faces).", ""]
    L.append("| dist | N | field | div(D) | div(B.C) | div(R) | dface(exact) | dface(exact on R) |")
    L.append("|---|---|---|---|---|---|---|---|")
    for dname in a.dists:
        for n in a.ns:
            x = {'lattice': lambda: p0.lattice(n), 'jitter': lambda: p0.jittered(n, 0.3), 'glass': lambda: p0.glass(n)}[dname]()
            geo = p0.Geometry(x, a.ngb, verbose=False)
            for fname in a.fields:
                fld = p0.make_field(fname, geo)
                Brms = np.sqrt(np.mean((fld['B']**2).sum(1)))
                grads, _ = p0.limited_gradients(geo, fld['B'])
                QR, QL = p0.face_states(geo, fld['B'], grads['std'])
                D, _ = p0.face_divergence(geo, QR, QL)
                BC = (fld['B']*geo.closure).sum(1)
                R = D - BC
                def dface(rhs):
                    if np.linalg.norm(rhs) < 1e-12*np.linalg.norm(D) + 1e-300: return 0.0
                    c, _, _ = p0.pcg(geo, rhs)
                    dBn = -0.5*(c[geo.I] - c[geo.J])*geo.Dlen2*geo.Alen      # change of B_face.n (notes eq. 4)
                    return np.sqrt(np.mean(dBn**2))/Brms
                L.append(f"| {dname} | {n}^3 | {fname} | {rmsn(geo, D, Brms):.3e} | {rmsn(geo, BC, Brms):.3e} | {rmsn(geo, R, Brms):.3e} | {dface(D):.3e} | {dface(R):.3e} |")
                print(L[-1], flush=True)
    outp = os.path.join(here, f"RESULTS_closure_split_ngb{a.ngb}.md")
    open(outp, 'w').write("\n".join(L) + "\n")
    print(f"# report -> {outp}")
