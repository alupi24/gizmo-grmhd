#!/usr/bin/env python3
"""
Exact face closure by an antisymmetric correction of the effective faces (follow-up to phase0_projection.py).

Given GIZMO's faces A_ij (outward from i, A_ji = -A_ij) with closure error C_i = sum_j A_ij, find the correction
    A~_ij = A_ij - k_ij (lambda_i - lambda_j),     lambda_i a vector per particle, k_ij = k_ji > 0,
that is the minimum of sum_ij |dA_ij|^2 / k_ij subject to sum_j A~_ij = 0 for every i. The constraint gives
    sum_j k_ij (lambda_i - lambda_j) = C_i        (one weighted graph Laplacian, three right-hand sides),
which is solvable because sum_i C_i = 0 exactly. A~ stays antisymmetric (conservation) and closes to round-off.
Weights tested: k = |A_ij| (smallest relative change of each face) and k = 1.

Diagnostics: closure before/after; rms |dA|/|A|; faces whose orientation flips (A~.dp > 0 in GIZMO's convention);
first-order (linear) consistency  E_i = sum_j 1/2 (x_j - x_i) (x) A_ij  versus  V_i I  (the property that makes the
face sum exact for the divergence of a linear field); face-sum divergence of tg3d with the corrected faces; Jacobi
sweeps needed to reduce the closure by 10x and 100x from a cold start.
"""
import os, copy, importlib.util
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('p0', os.path.join(here, 'phase0_projection.py'))
p0 = importlib.util.module_from_spec(spec); spec.loader.exec_module(p0)

class Lap:
    """weighted graph Laplacian with the interface p0.pcg expects (Rc, Wsum, N)."""
    def __init__(self, geo, k):
        self.N, self.I, self.J, self.w = geo.N, geo.I, geo.J, k
        self.Wsum = p0.sadd(geo.N, geo.I, k) + p0.sadd(geo.N, geo.J, k)
    def Rc(self, c):
        return self.Wsum*c - (p0.sadd(self.N, self.I, self.w*c[self.J]) + p0.sadd(self.N, self.J, self.w*c[self.I]))

def closure_vec(geo, A):
    return p0.sadd(geo.N, geo.I, A) - p0.sadd(geo.N, geo.J, A)

def closure_rel(geo, A):
    C = closure_vec(geo, A)
    Al = np.linalg.norm(A, axis=1); S = p0.sadd(geo.N, geo.I, Al) + p0.sadd(geo.N, geo.J, Al)
    return np.sqrt(np.mean((np.linalg.norm(C, axis=1)/S)**2))

def linear_consistency(geo, A):
    r = -geo.dp                                      # x_J - x_I
    o = 0.5*r[:, :, None]*A[:, None, :]              # pair contribution to I; to J: 1/2 (x_I-x_J) (x) A_JI = same
    E = p0.sadd(geo.N, geo.I, o) + p0.sadd(geo.N, geo.J, o)
    dev = E - geo.V[:, None, None]*np.eye(3)[None]
    return np.sqrt(np.mean((np.linalg.norm(dev, axis=(1, 2))/geo.V)**2))

def corrected_faces(geo, A, k):
    lap = Lap(geo, k)
    C = closure_vec(geo, A)
    lam = np.stack([p0.pcg(lap, C[:, a], tol=1e-13)[0] for a in range(3)], 1)
    return A - k[:, None]*(lam[geo.I] - lam[geo.J]), lap, C

def jacobi_sweeps(geo, A, lap, C, targets=(0.1, 0.01)):
    lam = np.zeros((geo.N, 3)); c0 = closure_rel(geo, A); out = {}
    for s in range(1, 400):
        R = C - np.stack([lap.Rc(lam[:, a]) for a in range(3)], 1)
        lam = lam + R/lap.Wsum[:, None]
        cr = closure_rel(geo, A - lap.w[:, None]*(lam[geo.I] - lam[geo.J]))
        for t in targets:
            if t not in out and cr < t*c0: out[t] = s
        if len(out) == len(targets): break
    return out

def tg3d_div(geo, A):
    fld = p0.make_field('tg3d', geo)
    Brms = np.sqrt(np.mean((fld['B']**2).sum(1)))
    g = copy.copy(geo); g.A = A
    grads, _ = p0.limited_gradients(geo, fld['B'])
    QR, QL = p0.face_states(geo, fld['B'], grads['std'])
    D, _ = p0.face_divergence(g, QR, QL)
    return p0.div_metric(geo, D, Brms)

if __name__ == '__main__':
    L = ["# Exact face closure by an antisymmetric correction A~_ij = A_ij - k_ij (lambda_i - lambda_j)", "",
         "`closure` = rms |sum_j A_ij| / sum_j |A_ij|; `dA` = rms |A~ - A| / rms |A|; `flipped` = faces whose orientation reverses; "
         "`lin. cons.` = rms || sum_j 1/2 (x_j - x_i) (x) A_ij - V_i I || / V_i; `div tg3d` = face-sum divergence of the "
         "divergence-free tg3d field, rms dx|D|/V/rms|B| (no B projection); `Jacobi 10x/100x` = sweeps to reduce the closure "
         "10x / 100x from a cold start (k = |A|).", "",
         "| set | N | faces | closure | dA | flipped | lin. cons. | div tg3d | Jacobi 10x / 100x |",
         "|---|---|---|---|---|---|---|---|---|"]
    for dname in ('glass', 'jitter'):
        for n in (16, 24, 32):
            x = {'jitter': lambda: p0.jittered(n, 0.3), 'glass': lambda: p0.glass(n)}[dname]()
            geo = p0.Geometry(x, 32, verbose=False)
            A0 = geo.A
            rows = [('standard', A0, 0.0, 0, None)]
            for kname, k in (('closed, k=|A|', geo.Alen.copy()), ('closed, k=1', np.ones(geo.P))):
                A1, lap, C = corrected_faces(geo, A0, k)
                dA = np.sqrt(np.mean(((A1 - A0)**2).sum(1)))/np.sqrt(np.mean((A0**2).sum(1)))
                flipped = int(np.sum((A1*A0).sum(1) <= 0))
                js = jacobi_sweeps(geo, A0, lap, C) if kname == 'closed, k=|A|' else None
                rows.append((kname, A1, dA, flipped, js))
            for kname, A, dA, flipped, js in rows:
                jtxt = f"{js.get(0.1)} / {js.get(0.01)}" if js else ""
                L.append(f"| {dname} | {n}^3 | {kname} | {closure_rel(geo, A):.2e} | {dA:.2e} | {flipped} | {linear_consistency(geo, A):.3e} | {tg3d_div(geo, A):.3e} | {jtxt} |")
                print(L[-1], flush=True)
    outp = os.path.join(here, "RESULTS_exact_closure_ngb32.md")
    open(outp, 'w').write("\n".join(L) + "\n")
    print(f"# report -> {outp}")
