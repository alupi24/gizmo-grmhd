#!/usr/bin/env python3
"""
Follow-up to phase0_projection.py: does the projection remove the right thing?

For B = B0 + grad(phi) the solenoidal part is B0. A consistent projection should leave face values B_face.n -> B0.n
as resolution increases. Two choices of the face increment (both one scalar c_i per particle, applied after the limiter):
  Tu et al.: delta B'_i(ij) = -c_i |D_ij|^2 A_ij           -> w_ij = 1/2 |D_ij|^2 |A_ij|^2     (notes eqs. 3-6)
  TPFA     : delta B'_i(ij) = -(2 c_i / r_ij) n_ij         -> w_ij = |A_ij| / r_ij             (two-point flux Laplacian)
Both are solved exactly (PCG 1e-12). `recovery` = rms|B_face.n - B0.n| / rms|grad(phi)(x_mid).n| (1 = untouched, 0 = removed).
On irregular particle sets the full face sum also contains the face-closure term B_i . sum_j A_ij; we therefore also project
only the remainder R = D - B.C (what would be left with closed faces).

Usage: python3 phase0_weights.py [--ns 16 24 32] [--dists lattice glass]
"""
import argparse, os, copy, importlib.util
import numpy as np
import scipy.sparse as sp

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('p0', os.path.join(here, 'phase0_projection.py'))
p0 = importlib.util.module_from_spec(spec); spec.loader.exec_module(p0)

def reweight(geo, mode):
    g = copy.copy(geo)
    if mode == 'tpfa':
        g.Dlen2 = 2.0/(geo.r*geo.Alen)              # increment vector = c * Dlen2 * A = (2c/r) n
    g.w = 0.5*g.Dlen2*g.Alen**2
    g.Wsum = p0.sadd(g.N, g.I, g.w) + p0.sadd(g.N, g.J, g.w)
    g.Wm = sp.coo_matrix((np.concatenate([g.w, g.w]), (np.concatenate([g.I, g.J]), np.concatenate([g.J, g.I]))), shape=(g.N, g.N)).tocsr()
    return g

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--ns', type=int, nargs='+', default=[16, 24, 32])
    ap.add_argument('--ngb', type=int, default=32)
    ap.add_argument('--dists', nargs='+', default=['lattice', 'glass'])
    ap.add_argument('--fields', nargs='+', default=['grad_short', 'grad_long', 'tg3d'])
    a = ap.parse_args()
    L = ["# Weighting of the projection: Tu et al. versus two-point (TPFA)", "",
         "Exact solves. `recovery` (gradient fields): rms|B_face.n - B0.n| / rms|grad(phi).n|, 0 = gradient part fully removed. "
         "`face err` (tg3d): rms|B_face.n - B_exact.n| / rms|B|. `rhs = D`: full face-sum divergence; `rhs = R`: closure term removed (closed-face proxy). "
         "Field widths are fixed in units of the mean spacing (grad_short: 1.5 dx, grad_long: 4 dx), so the gradient part is equally resolved at every N.", ""]
    L.append("| dist | N | field | no projection | Tu, rhs=D | TPFA, rhs=D | Tu, rhs=R | TPFA, rhs=R |")
    L.append("|---|---|---|---|---|---|---|---|")
    for dname in a.dists:
        for n in a.ns:
            x = {'lattice': lambda: p0.lattice(n), 'jitter': lambda: p0.jittered(n, 0.3), 'glass': lambda: p0.glass(n)}[dname]()
            geo = p0.Geometry(x, a.ngb, verbose=False)
            gw = {'tu': reweight(geo, 'tu'), 'tpfa': reweight(geo, 'tpfa')}
            for fname in a.fields:
                fld = p0.make_field(fname, geo)
                grads, _ = p0.limited_gradients(geo, fld['B'])
                QR, QL = p0.face_states(geo, fld['B'], grads['std'])
                D, _ = p0.face_divergence(geo, QR, QL)
                R = D - (fld['B']*geo.closure).sum(1)
                vals = [p0.face_metric(geo, fld, QR, QL)]
                for rhs in (D, R):
                    for m in ('tu', 'tpfa'):
                        g = gw[m]
                        c, _, _ = p0.pcg(g, rhs)
                        QRm, QLm = p0.apply_correction(g, QR, QL, c)
                        vals.append(p0.face_metric(geo, fld, QRm, QLm))
                L.append(f"| {dname} | {n}^3 | {fname} | " + " | ".join(f"{v:.3f}" for v in vals) + " |")
                print(L[-1], flush=True)
    outp = os.path.join(here, f"RESULTS_weights_ngb{a.ngb}.md")
    open(outp, 'w').write("\n".join(L) + "\n")
    print(f"# report -> {outp}")
