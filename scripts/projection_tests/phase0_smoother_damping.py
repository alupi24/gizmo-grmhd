#!/usr/bin/env python3
"""
Follow-up to phase0_projection.py: damping of the smoothers at equal cost.

Weighted Jacobi:      c <- c + omega   r / W                     (cost 1 PPE per sweep)
Damped stencil solve: c <- c + omega_S (c_local - c)              (cost 1 + intra/pairs PPE per pass)
where c_local is the restricted-Schwarz value of phase0_projection.stencil (keep only the centre of each local solve).
Reports div(variant)/div(none) as a function of cost, with GIZMO faces/gradients/limiters as in phase0_projection.py.

Usage: python3 phase0_smoother_damping.py [--n 24] [--dists glass jitter] [--fields tg3d grad_long bump]
"""
import argparse, os, sys, importlib.util
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('p0', os.path.join(here, 'phase0_projection.py'))
p0 = importlib.util.module_from_spec(spec); spec.loader.exec_module(p0)
spec_w = importlib.util.spec_from_file_location('pw', os.path.join(here, 'phase0_weights.py'))
pw = importlib.util.module_from_spec(spec_w); spec_w.loader.exec_module(pw)

def jacobi_w(geo, D, K, omega):
    c = np.zeros(geo.N); out = []
    for k in range(K):
        c = c + omega*(D - geo.Rc(c))/geo.Wsum
        out.append(c.copy())
    return out

def stencil_w(geo, D, passes, omega_s):
    c = np.zeros(geo.N); out = []
    for k in range(passes):
        cl = p0.stencil(geo, D, 1, c0=c)
        c = c + omega_s*(cl - c)
        out.append(c.copy())
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=24)
    ap.add_argument('--ngb', type=int, default=32)
    ap.add_argument('--dists', nargs='+', default=['glass', 'jitter'])
    ap.add_argument('--fields', nargs='+', default=['tg3d', 'grad_long', 'bump'])
    ap.add_argument('--weights', nargs='+', default=['tu', 'tpfa'])
    a = ap.parse_args()
    lines = [f"# Smoother damping at equal cost: N = {a.n}^3, N_ngb = {a.ngb}, GIZMO limiters on, gradient std",
             "", "Entries: div(variant) / div(no projection) (smaller is better); cost in pair-pass equivalents (PPE).", ""]
    for dname in a.dists:
        x = {'lattice': lambda: p0.lattice(a.n), 'jitter': lambda: p0.jittered(a.n, 0.3), 'glass': lambda: p0.glass(a.n)}[dname]()
        geo0 = p0.Geometry(x, a.ngb)
        for wname in a.weights:
          geo = pw.reweight(geo0, wname); geo._stencil = None
          st = geo.stencil_setup()
          for fname in a.fields:
              fld = p0.make_field(fname, geo)
              Brms = np.sqrt(np.mean((fld['B']**2).sum(1)))
              grads, _ = p0.limited_gradients(geo, fld['B'])
              QR, QL = p0.face_states(geo, fld['B'], grads['std'])
              D0, _ = p0.face_divergence(geo, QR, QL)
              d0 = p0.div_metric(geo, D0, Brms)
              def rel(c):
                  Dm, _ = p0.face_divergence(geo, *p0.apply_correction(geo, QR, QL, c))
                  return p0.div_metric(geo, Dm, Brms)/d0
              lines.append(f"## {dname} / weights {wname} / {fname}   (no projection: div = {d0:.3e}; stencil pass = {st['cost']:.1f} PPE)\n")
              lines.append("| smoother | " + " | ".join(f"{k} PPE" for k in (1, 2, 3, 5, 10, 15, 30)) + " |")
              lines.append("|---|" + "---|"*7)
              for om in (0.5, 2/3, 0.8, 1.0):
                  cs = jacobi_w(geo, D0, 30, om)
                  lines.append(f"| Jacobi omega={om:.2f} | " + " | ".join(f"{rel(cs[k-1]):.3f}" for k in (1, 2, 3, 5, 10, 15, 30)) + " |")
              for oms in (1.0, 0.8, 0.6, 0.5):
                  cs = stencil_w(geo, D0, 2, oms)
                  cost1, cost2 = st['cost'], 2*st['cost']
                  lines.append(f"| stencil omega_S={oms:.1f} | 1 pass ({cost1:.1f} PPE): {rel(cs[0]):.3f} | 2 passes ({cost2:.1f} PPE): {rel(cs[1]):.3f} | | | | | |")
              lines.append("")
              print(f"  done {dname}/{wname}/{fname}", flush=True)
    outp = os.path.join(here, f"RESULTS_smoother_damping_n{a.n}_ngb{a.ngb}.md")
    open(outp, 'w').write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"# report -> {outp}")
