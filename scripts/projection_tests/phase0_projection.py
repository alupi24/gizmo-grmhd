#!/usr/bin/env python3
"""
Phase-0 static comparison for the localized divergence projection
(divfree/hydro/docs/LOCAL_PROJECTION_notes.pdf).

Everything is evaluated with GIZMO's own MFM discretisation, replicated from the source:
  - smoothing length: 4pi/3 h^3 sum_j W(r_ij,h) = N_ngb, self included        (hydro/density.c:314,546)
  - V_i = 1/sum_j W(r_ij,h_i)                                                    (hydro/density.c:1005)
  - E_i = sum_{j!=i, r<h_i} W(r_ij,h_i) dp dp^T, NV_T = E^{-1} (+ conditioning)  (hydro/density.c:344-353,567-581)
  - pairs: r_ij < max(h_i,h_j)                                                   (hydro/hydro_evaluate.h:190, gradients.c:1899)
  - gradient: G_i[a,:] = NV_T_i . sum_{r<h_i} W(r,h_i) (x_j-x_i) (B_j-B_i)_a     (hydro/gradients.c:2270-2275, construct_gradient)
  - B slope limiter per row, alim=0.25, h_lim=max(h,MaxDistance), stol=0       (hydro/gradients.c:1448,1530-1544, local_slopelimiter)
  - faces: Face_Area_Vec = wk_i wt_i NV_T_i.dp + wk_j wt_j NV_T_j.dp, wt rule,  (hydro/compute_finitevol_faces.h:29-40,93-101)
           RSPH fallback if Face.dp < 0
  - face states at the midpoint, pairwise limiter of reconstruct_face_states, mode 1 (MAGNETIC: 0.5 / 0.375)  (hydro/reimann.h:141-245)
  - divergence: D_i = sum_j 1/2 (B'_i + B'_j) . A_ij with A_ij = -Face_Area_Vec (outward from i)  (hydro_core_meshless.h:488, hydro_evaluate.h:496)
Flat space (sqrt(gamma)=1): the densitized and primitive fields coincide.

Variants (notes, Sect. 4): gradient std (GIZMO) or tl (traceless + scalar limiter), each followed by
  none | Jacobi K sweeps (omega=2/3) | star (Phase 2b) | stencil-local solve, 1 or 2 passes (Phase 2) | exact (PCG to 1e-12, Phase 3).
Costs are given in pair-pass equivalents (PPE): one pass over all pairs with face evaluation, i.e. one constrained-gradient-type iteration.

Usage: python3 phase0_projection.py [--n 32] [--ngb 32] [--dists lattice jitter glass] [--fields ...] [--no-limiter]
"""
import argparse, sys, time, csv, os
import numpy as np
from scipy.spatial import cKDTree
import scipy.sparse as sp

KNORM = 8.0/np.pi                  # cubic spline, 3D (kernel.h, KERNEL_FUNCTION 3)
COND_DANGER = 1.0e7               # CONDITION_NUMBER_DANGER for MAGNETIC && !COOLING (allvars.h:1355)
OMEGA = 2.0/3.0
KAPPA_CAP = 0.5
rng = np.random.default_rng(20261013)

# ---------------------------------------------------------------- kernel
def W(r, h):
    q = r/h
    s = np.where(q < 0.5, 1.0 + 6.0*(q - 1.0)*q*q, np.where(q < 1.0, 2.0*(1.0-q)**3, 0.0))
    return KNORM/h**3 * s

def dWdr(r, h):
    q = r/h
    s = np.where(q < 0.5, q*(18.0*q - 12.0), np.where(q < 1.0, -6.0*(1.0-q)**2, 0.0))
    return KNORM/h**4 * s

# ---------------------------------------------------------------- scatter helpers
def sadd(N, idx, vals):
    vals = np.asarray(vals, dtype=float)
    if vals.ndim == 1:
        return np.bincount(idx, weights=vals, minlength=N)
    flat = vals.reshape(len(idx), -1)
    out = np.empty((N, flat.shape[1]))
    for c in range(flat.shape[1]):
        out[:, c] = np.bincount(idx, weights=flat[:, c], minlength=N)
    return out.reshape((N,) + vals.shape[1:])

def smax(N, idx, vals, init):
    out = np.full((N,) + vals.shape[1:], init, dtype=float)
    np.maximum.at(out, idx, vals)
    return out

def smin(N, idx, vals, init):
    out = np.full((N,) + vals.shape[1:], init, dtype=float)
    np.minimum.at(out, idx, vals)
    return out

# ---------------------------------------------------------------- particle sets
def lattice(n):
    g = (np.arange(n) + 0.5)/n
    return np.stack(np.meshgrid(g, g, g, indexing='ij'), -1).reshape(-1, 3)

def jittered(n, amp=0.3):
    return (lattice(n) + amp/n*rng.uniform(-1, 1, (n**3, 3))) % 1.0

def glass(n, iters=150, k=32, step=0.05):
    x = jittered(n, amp=0.5)
    dx_mean = 1.0/n
    for it in range(iters):
        tree = cKDTree(x, boxsize=1.0)
        d, j = tree.query(x, k=k+1)
        d, j = d[:, 1:], j[:, 1:]
        dvec = x[:, None, :] - x[j]
        dvec -= np.round(dvec)
        dn = np.linalg.norm(dvec, axis=2)
        f = (dx_mean/np.maximum(dn, 1e-3*dx_mean))**2
        push = (dvec*(f/dn)[:, :, None]).sum(1)
        pn = np.linalg.norm(push, axis=1)
        fac = step*dx_mean/np.maximum(pn, 1e-30)*np.minimum(1.0, pn/pn.mean())
        x = (x + push*fac[:, None]) % 1.0
    return x

# ---------------------------------------------------------------- geometry (GIZMO MFM)
class Geometry:
    def __init__(self, x, nngb, periodic=True, verbose=True):
        t0 = time.time()
        self.x, self.N, self.periodic = x, len(x), periodic
        box = 1.0 if periodic else None
        tree = cKDTree(x, boxsize=box)
        # --- smoothing lengths (density.c): 4pi/3 h^3 sum_j W(r_ij,h) = N_ngb, self included
        kq = min(self.N, 4*nngb if periodic else 16*nngb)
        d, _ = tree.query(x, k=kq)
        lo, hi = d[:, max(1, nngb//4)], d[:, -1]*0.999999
        f_hi = 4*np.pi/3*hi**3*W(d, hi[:, None]).sum(1)
        if np.any(f_hi < nngb):
            raise RuntimeError("neighbour search radius too small for the requested N_ngb")
        for _ in range(100):
            mid = 0.5*(lo + hi)
            f = 4*np.pi/3*mid**3*W(d, mid[:, None]).sum(1)
            lo = np.where(f < nngb, mid, lo); hi = np.where(f >= nngb, mid, hi)
        self.h = hi
        self.Neff = 4*np.pi/3*self.h**3*W(d, self.h[:, None]).sum(1)
        self.V = 1.0/W(d, self.h[:, None]).sum(1)
        self.dx = self.V**(1.0/3.0)
        # --- pairs r < max(h_i, h_j)
        lists = tree.query_ball_point(x, r=self.h*(1 - 1e-12), return_sorted=False)
        ii = np.repeat(np.arange(self.N), [len(l) for l in lists]); jj = np.concatenate([np.asarray(l, dtype=np.int64) for l in lists])
        keep = ii != jj
        a, b = np.minimum(ii[keep], jj[keep]), np.maximum(ii[keep], jj[keep])
        key = np.unique(a*self.N + b)
        self.I, self.J = key//self.N, key % self.N
        dp = x[self.I] - x[self.J]
        if periodic: dp -= np.round(dp)
        self.dp = dp                                   # kernel.dp = x_i - x_j (i = I)
        self.r = np.linalg.norm(dp, axis=1)
        self.P = len(self.I)
        hI, hJ = self.h[self.I], self.h[self.J]
        self.wkI = np.where(self.r < hI, W(self.r, hI), 0.0)
        self.wkJ = np.where(self.r < hJ, W(self.r, hJ), 0.0)
        self.dwkI = np.where(self.r < hI, dWdr(self.r, hI), 0.0)
        self.dwkJ = np.where(self.r < hJ, dWdr(self.r, hJ), 0.0)
        # --- NV_T (density.c): E = sum W dp dp^T over r<h_i, then conditioning on the dimensionless matrix E*h, NV_T = E^{-1}
        outer = dp[:, :, None]*dp[:, None, :]
        E = sadd(self.N, self.I, self.wkI[:, None, None]*outer) + sadd(self.N, self.J, self.wkJ[:, None, None]*outer)
        Ed = E*self.h[:, None, None]
        thr = 10.0*COND_DANGER
        self.n_condfix = 0
        for _ in range(5):
            cond = np.linalg.cond(Ed)
            bad = cond >= thr
            if not bad.any(): break
            self.n_condfix += int(bad.sum())
            tr = np.trace(Ed[bad], axis1=1, axis2=2)
            Ed[bad] += (1.05*(tr/3.0)/thr)[:, None, None]*np.eye(3)[None]
        self.cond = np.linalg.cond(Ed)
        self.NV = np.linalg.inv(Ed)*self.h[:, None, None]
        # --- faces (compute_finitevol_faces.h)
        VI, VJ = self.V[self.I], self.V[self.J]
        PsI, PsJ = VI**(1/3.), VJ**(1/3.)
        wtI, wtJ = VI.copy(), VJ.copy()
        big = (np.abs(VI - VJ)/np.minimum(VI, VJ)/3.0 > 5.5) & ((self.r > hI) | (self.r > hJ))
        wsym = (VI*PsJ + VJ*PsI)/(PsI + PsJ)
        wtI[big] = wsym[big]; wtJ[big] = wsym[big]
        self.n_bigwt = int(big.sum())
        NVdpI = np.einsum('pab,pb->pa', self.NV[self.I], dp)
        NVdpJ = np.einsum('pab,pb->pa', self.NV[self.J], dp)
        F = (self.wkI*wtI)[:, None]*NVdpI + (self.wkJ*wtJ)[:, None]*NVdpJ
        fdot = (F*dp).sum(1)
        cnum2 = COND_DANGER**2 - self.cond[self.I]**2
        fb = (self.cond[self.J]**2 > 1e12 + cnum2) | (fdot < 0)
        if fb.any():
            Fn = -(wtI*VI*self.dwkI + wtJ*VJ*self.dwkJ)/self.r
            F[fb] = (Fn[:, None]*dp)[fb]
        self.n_fallback = int(fb.sum())
        # S1: the same face evaluated from J's side (roles swapped, dp -> -dp) must be exactly -F
        F_swap = (self.wkJ*wtJ)[:, None]*np.einsum('pab,pb->pa', self.NV[self.J], -dp) + (self.wkI*wtI)[:, None]*np.einsum('pab,pb->pa', self.NV[self.I], -dp)
        if fb.any():
            Fn2 = -(wtJ*VJ*self.dwkJ + wtI*VI*self.dwkI)/self.r
            F_swap[fb] = (Fn2[:, None]*(-dp))[fb]
        self.S1 = np.max(np.abs(F_swap + F))/np.max(np.abs(F))
        self.A = -F                                    # outward from I
        self.Alen = np.linalg.norm(self.A, axis=1)
        self.nrm = self.A/self.Alen[:, None]
        self.closure = sadd(self.N, self.I, self.A) - sadd(self.N, self.J, self.A)   # sum_j A_ij (outward)
        # --- projection weights (notes eq. 5): w = 1/2 |D|^2 |A|^2, |D| = r/2
        self.Dlen2 = 0.25*self.r**2
        self.w = 0.5*self.Dlen2*self.Alen**2
        self.Wsum = sadd(self.N, self.I, self.w) + sadd(self.N, self.J, self.w)
        self.incr_max = np.maximum(smax(self.N, self.I, self.Dlen2*self.Alen, 0.0), smax(self.N, self.J, self.Dlen2*self.Alen, 0.0))
        self.MaxDistance = np.maximum(smax(self.N, self.I, self.r, 0.0), smax(self.N, self.J, self.r, 0.0))
        xm = x[self.I] - 0.5*dp
        self.xmid = xm % 1.0 if periodic else xm
        self.Wm = sp.coo_matrix((np.concatenate([self.w, self.w]), (np.concatenate([self.I, self.J]), np.concatenate([self.J, self.I]))), shape=(self.N, self.N)).tocsr()
        self.nngb_pairs = np.diff(self.Wm.indptr)
        self.t_geo = time.time() - t0
        self._stencil = None
        if verbose:
            print(f"  geometry: N={self.N} pairs={self.P} (mean {2*self.P/self.N:.1f} partners/particle) "
                  f"h/dx={np.mean(self.h/self.dx):.3f} cond median={np.median(self.cond):.2f} max={self.cond.max():.2f} "
                  f"condfix={self.n_condfix} RSPH-fallback faces={self.n_fallback} sym-weight faces={self.n_bigwt} [{self.t_geo:.1f}s]", flush=True)

    # ---- Laplacian of the projection system (notes eq. 6)
    def Rc(self, c):
        return self.Wsum*c - (sadd(self.N, self.I, self.w*c[self.J]) + sadd(self.N, self.J, self.w*c[self.I]))

    # ---- stencil-local solve, geometry-only precomputation (z = M_i^{-1} e_0, y = sub_i z)
    def stencil_setup(self):
        if self._stencil is not None: return self._stencil
        t0 = time.time()
        rows, cols, zf, yf = [], [], [], []
        intra = 0
        solve_err = 0.0
        Wm = self.Wm
        for i in range(self.N):
            S = np.concatenate(([i], Wm.indices[Wm.indptr[i]:Wm.indptr[i+1]]))
            sub = Wm[S][:, S].toarray()
            M = np.diag(self.Wsum[S]) - sub
            e0 = np.zeros(len(S)); e0[0] = 1.0
            z = np.linalg.solve(M, e0)
            if i % 997 == 0:   # S10: local solve accuracy on a sample
                solve_err = max(solve_err, np.max(np.abs(M @ z - e0)))
            rows.append(np.full(len(S), i)); cols.append(S); zf.append(z); yf.append(sub @ z)
            intra += np.count_nonzero(sub)//2
        self._stencil = dict(rows=np.concatenate(rows), cols=np.concatenate(cols), z=np.concatenate(zf), y=np.concatenate(yf),
                             cost=1.0 + intra/self.P, S10=solve_err, t=time.time()-t0)
        return self._stencil

# ---------------------------------------------------------------- fields
def make_field(name, geo):
    x, xm, dxm = geo.x, geo.xmid, np.mean(geo.dx)
    k = 2*np.pi
    B0 = np.array([0.3, -0.5, 1.0])
    out = dict(name=name, kind='divfree', B0=None)
    if name == 'uniform':
        f = lambda p: np.tile(B0, (len(p), 1))
        out['B0'] = B0
    elif name == 'tg3d':
        f = lambda p: np.stack([np.sin(k*p[:, 0])*np.cos(k*p[:, 1])*np.cos(k*p[:, 2]),
                                np.cos(k*p[:, 0])*np.sin(k*p[:, 1])*np.cos(k*p[:, 2]),
                                -2*np.cos(k*p[:, 0])*np.cos(k*p[:, 1])*np.sin(k*p[:, 2])], 1)
    elif name == 'abc':
        f = lambda p: np.stack([np.sin(k*p[:, 2]) + np.cos(k*p[:, 1]), np.sin(k*p[:, 0]) + np.cos(k*p[:, 2]),
                                np.sin(k*p[:, 1]) + np.cos(k*p[:, 0])], 1)
    elif name in ('grad_short', 'grad_long'):
        sig = (1.5 if name == 'grad_short' else 4.0)*dxm
        amp = 0.5*np.linalg.norm(B0)*sig*np.exp(0.5)
        def gphi(p):
            d = p - 0.5
            if geo.periodic: d -= np.round(d)
            phi = amp*np.exp(-(d**2).sum(1)/(2*sig**2))
            return -(d/sig**2)*phi[:, None]
        f = lambda p: B0[None, :] + gphi(p)
        out.update(kind='gradient', B0=B0, gphi=gphi, sigma=sig)
    elif name == 'bump':
        R = 5.0*dxm
        def f(p):
            d = p - 0.5
            if geo.periodic: d -= np.round(d)
            s = np.clip(1.0 - (d**2).sum(1)/R**2, 0.0, None)**2
            return np.stack([s, 0*s, 0*s], 1)
        out.update(kind='bump', R=R)
    else:
        raise ValueError(name)
    out['f'] = f
    out['B'] = f(x)
    out['Bmid'] = f(xm)
    return out

# ---------------------------------------------------------------- gradients and limiters
def mls_gradient(geo, B):
    dB = B[geo.J] - B[geo.I]
    m = -geo.dp                                         # x_J - x_I for the I side; for the J side (x_I-x_J)(B_I-B_J) = (-dp)(dB) too
    S = sadd(geo.N, geo.I, geo.wkI[:, None, None]*m[:, :, None]*dB[:, None, :]) + \
        sadd(geo.N, geo.J, geo.wkJ[:, None, None]*m[:, :, None]*dB[:, None, :])   # S[i, b(space), a(component)]
    G = np.einsum('ibc,ica->iab', geo.NV, S)           # G[i, a, b] = d_b B_a
    return G, dB

def minmax_dB(geo, dB):
    mx = np.maximum(smax(geo.N, geo.I, dB, 0.0), smax(geo.N, geo.J, -dB, 0.0))
    mn = np.minimum(smin(geo.N, geo.I, dB, 0.0), smin(geo.N, geo.J, -dB, 0.0))
    return mx, mn

def row_limiter_factors(geo, G, mx, mn):
    """local_slopelimiter(grad, valmax, valmin, alim2, h_lim, stol=0): factor per row (component)."""
    alim = np.where(geo.cond > 100, np.minimum(0.5, 0.25 + 0.25*(geo.cond - 100)/100), 0.25)   # alim2 with divB_prev = 0
    hlim = np.maximum(geo.h, geo.MaxDistance)
    dabs = np.linalg.norm(G, axis=2)
    absmin = np.minimum(np.abs(mx), np.abs(mn))
    with np.errstate(divide='ignore', invalid='ignore'):
        cfac = np.where(dabs > 0, absmin/(alim[:, None]*hlim[:, None]*dabs), 1.0)
    return np.minimum(cfac, 1.0)

def limited_gradients(geo, B, limiter=True):
    G, dB = mls_gradient(geo, B)
    tr = np.trace(G, axis1=1, axis2=2)
    Gtl = G - (tr/3.0)[:, None, None]*np.eye(3)[None]
    if limiter:
        mx, mn = minmax_dB(geo, dB)
        fs = row_limiter_factors(geo, G, mx, mn)                  # GIZMO: per row
        Gstd = G*fs[:, :, None]
        ft = row_limiter_factors(geo, Gtl, mx, mn).min(1)         # traceless: one scalar for the whole tensor
        Gtl = Gtl*ft[:, None, None]
    else:
        Gstd = G
    return {'std': Gstd, 'tl': Gtl}, G

def face_limiter(Qi, Qj, QR, QL):
    """reconstruct_face_states, mode 1, MAGNETIC (fac_minmax=0.5, fac_meddev=0.375); R = i side, L = j side."""
    Qmed = 0.5*(Qi + Qj)
    Qmax = np.maximum(Qi, Qj); Qmin = np.minimum(Qi, Qj)
    fac = 0.5*(Qmax - Qmin)
    Qmax_eff = Qmax + fac; Qmin_eff = Qmin - fac
    with np.errstate(divide='ignore', invalid='ignore'):
        c1 = (Qmax < 0) & (Qmax_eff > 0)
        Qmax_eff = np.where(c1, Qmax*Qmax/(Qmax - (Qmax_eff - Qmax)), Qmax_eff)
        c2 = (Qmin > 0) & (Qmin_eff < 0)
        Qmin_eff = np.where(c2, Qmin*Qmin/(Qmin + (Qmin - Qmin_eff)), Qmin_eff)
    fac2 = 0.375*(Qmax - Qmin)
    Qmed_max = np.minimum(Qmed + fac2, Qmax_eff)
    Qmed_min = np.maximum(Qmed - fac2, Qmin_eff)
    lt = Qi < Qj
    R1 = np.where(QR < Qmin_eff, Qmin_eff, QR); R1 = np.where(R1 > Qmed_max, Qmed_max, R1)
    L1 = np.where(QL > Qmax_eff, Qmax_eff, QL); L1 = np.where(L1 < Qmed_min, Qmed_min, L1)
    R2 = np.where(QR > Qmax_eff, Qmax_eff, QR); R2 = np.where(R2 < Qmed_min, Qmed_min, R2)
    L2 = np.where(QL < Qmin_eff, Qmin_eff, QL); L2 = np.where(L2 > Qmed_max, Qmed_max, L2)
    QR = np.where(lt, R1, R2); QL = np.where(lt, L1, L2)
    eq = Qi == Qj
    return np.where(eq, Qi, QR), np.where(eq, Qi, QL)

def face_states(geo, B, G, limiter=True):
    D_i = -0.5*geo.dp                                   # x_mid - x_I
    QR = B[geo.I] + np.einsum('pab,pb->pa', G[geo.I], D_i)
    QL = B[geo.J] + np.einsum('pab,pb->pa', G[geo.J], -D_i)
    if limiter:
        QR, QL = face_limiter(B[geo.I], B[geo.J], QR, QL)
    else:
        eq = B[geo.I] == B[geo.J]
        QR = np.where(eq, B[geo.I], QR); QL = np.where(eq, B[geo.I], QL)
    return QR, QL

def face_divergence(geo, QR, QL):
    flux = 0.5*((QR + QL)*geo.A).sum(1)
    return sadd(geo.N, geo.I, flux) - sadd(geo.N, geo.J, flux), flux

def apply_correction(geo, QR, QL, c):
    QRm = QR - (c[geo.I]*geo.Dlen2)[:, None]*geo.A       # notes eq. (12)
    QLm = QL + (c[geo.J]*geo.Dlen2)[:, None]*geo.A
    return QRm, QLm

# ---------------------------------------------------------------- solvers
def local_Bscale(geo, B):
    """max |B| over the particle and its partners (the notes' cap uses |B_i| alone, which vanishes where B does)."""
    Bmag = np.linalg.norm(B, axis=1)
    return np.maximum(Bmag, np.maximum(smax(geo.N, geo.I, Bmag[geo.J], 0.0), smax(geo.N, geo.J, Bmag[geo.I], 0.0)))

def increment_ratio(geo, c, Bscale):
    """largest face increment |c_i| |D|^2 |A| relative to the local field scale: what a cap would have to allow."""
    inc = np.abs(c)*geo.incr_max
    ok = Bscale > 0
    rat = inc[ok]/Bscale[ok]
    n_zeroB = int(np.sum((~ok) & (inc > 0)))
    return (rat.max() if rat.size else 0.0), (np.percentile(rat, 99) if rat.size else 0.0), n_zeroB

def jacobi(geo, D, K, c0=None, cap=None, history=None):
    c = np.zeros(geo.N) if c0 is None else c0.copy()
    nD = np.linalg.norm(D)
    for k in range(K):
        r = D - geo.Rc(c)
        c = c + OMEGA*r/geo.Wsum
        if cap is not None: c = np.clip(c, -cap, cap)
        if history is not None: history.append(np.linalg.norm(D - geo.Rc(c))/nD)
    return c

def star(geo, D, c0=None, cap=None):
    cold = np.zeros(geo.N) if c0 is None else c0
    T = geo.Wm @ cold
    cJ = (D + T)/geo.Wsum
    w2W_I = geo.w**2/geo.Wsum[geo.J]; w2W_J = geo.w**2/geo.Wsum[geo.I]
    s_w2W = sadd(geo.N, geo.I, w2W_I) + sadd(geo.N, geo.J, w2W_J)
    s_wcJ = sadd(geo.N, geo.I, geo.w*cJ[geo.J]) + sadd(geo.N, geo.J, geo.w*cJ[geo.I])
    c = (D + s_wcJ - s_w2W*cold)/(geo.Wsum - s_w2W)
    if cap is not None: c = np.clip(c, -cap, cap)
    return c

def stencil(geo, D, passes=1, c0=None, cap=None):
    st = geo.stencil_setup()
    c = np.zeros(geo.N) if c0 is None else c0.copy()
    for _ in range(passes):
        T = geo.Wm @ c
        c = np.bincount(st['rows'], weights=st['z']*(D + T)[st['cols']], minlength=geo.N) - \
            np.bincount(st['rows'], weights=st['y']*c[st['cols']], minlength=geo.N)
        if cap is not None: c = np.clip(c, -cap, cap)
    return c

def pcg(geo, D, tol=1e-12, maxit=20000):
    b = D - D.mean()                                    # remove the round-off null-space component (sum D = 0 up to round-off)
    x = np.zeros(geo.N); r = b.copy(); Minv = 1.0/geo.Wsum
    z = Minv*r; p = z.copy(); rz = r @ z; nb = np.linalg.norm(b)
    marks = {1e-1: None, 1e-2: None, 1e-3: None, 1e-6: None, 1e-12: None}
    for it in range(1, maxit + 1):
        Ap = geo.Rc(p); alpha = rz/(p @ Ap)
        x += alpha*p; r -= alpha*Ap
        rel = np.linalg.norm(r)/nb
        for t in marks:
            if marks[t] is None and rel < t: marks[t] = it
        if rel < tol: break
        z = Minv*r; rz_new = r @ z; p = z + (rz_new/rz)*p; rz = rz_new
    return x, marks, it

# ---------------------------------------------------------------- metrics
def div_metric(geo, D, Brms):
    return np.sqrt(np.mean((geo.dx*D/geo.V)**2))/Brms

def face_metric(geo, fld, QR, QL):
    Bn = 0.5*((QR + QL)*geo.nrm).sum(1)
    if fld['kind'] == 'gradient':
        gpart = (fld['gphi'](geo.xmid)*geo.nrm).sum(1)
        return np.sqrt(np.mean((Bn - (fld['B0'][None, :]*geo.nrm).sum(1))**2))/np.sqrt(np.mean(gpart**2))
    exact = (fld['Bmid']*geo.nrm).sum(1)
    return np.sqrt(np.mean((Bn - exact)**2))/np.sqrt(np.mean((fld['B']**2).sum(1)))

# ---------------------------------------------------------------- one distribution
GEOSTATS = []

def run_distribution(name, geo, fields, limiter, checks):
    rows = []
    st = geo.stencil_setup()
    sumA = sadd(geo.N, geo.I, geo.Alen) + sadd(geo.N, geo.J, geo.Alen)
    clo = np.linalg.norm(geo.closure, axis=1)
    GEOSTATS.append(dict(dist=name, N=geo.N, pairs=geo.P, partners=2*geo.P/geo.N, h_dx=np.mean(geo.h/geo.dx),
                         cond_med=np.median(geo.cond), cond_max=geo.cond.max(), condfix=geo.n_condfix, fallback=geo.n_fallback,
                         closure_rms=np.sqrt(np.mean((clo/geo.dx**2)**2)), closure_rel=np.sqrt(np.mean(clo**2))/np.mean(sumA),
                         stencil_cost=st['cost']))
    print(f"  stencil setup: cost {st['cost']:.1f} PPE per pass, local-solve accuracy {st['S10']:.1e} [{st['t']:.1f}s]", flush=True)
    checks.append((name, 'S0 h solve: max |N_eff - N_ngb|/N_ngb', np.max(np.abs(geo.Neff - geo.Neff.mean()))/geo.Neff.mean(), 1e-6))
    checks.append((name, 'S1 face antisymmetry: max|A_ij + A_ji|/max|A|', geo.S1, 1e-15))
    checks.append((name, 'S10 stencil local solve: max|M z - e0|', st['S10'], 1e-10))
    for fname in fields:
        fld = make_field(fname, geo)
        Brms = np.sqrt(np.mean((fld['B']**2).sum(1)))
        grads, Graw = limited_gradients(geo, fld['B'], limiter)
        Bscale = local_Bscale(geo, fld['B'])
        if fname == 'uniform':
            QR, QL = face_states(geo, fld['B'], grads['std'], limiter)
            D0, flux = face_divergence(geo, QR, QL)
            Fabs = sadd(geo.N, geo.I, np.abs(flux)) + sadd(geo.N, geo.J, np.abs(flux))
            ref = fld['B0'] @ geo.closure.T
            checks.append((f"{name}/{fname}", 'S6 uniform field: max|D_i - B0.sum_j A_ij| / max_i sum_j|flux_ij|', np.max(np.abs(D0 - ref))/np.max(Fabs), 1e-12))
        for gname, G in grads.items():
            if gname == 'tl':
                trn = np.abs(np.trace(G, axis1=1, axis2=2))/np.maximum(np.linalg.norm(G, axis=(1, 2)), 1e-300)
                checks.append((f"{name}/{fname}", 'S7 traceless gradient after scalar limiter: max|tr G|/|G|', trn.max(), 1e-12))
            QR, QL = face_states(geo, fld['B'], G, limiter)
            D0, flux = face_divergence(geo, QR, QL)
            # non-degenerate scale for the checks: total absolute flux through each particle's faces
            Fabs = sadd(geo.N, geo.I, np.abs(flux)) + sadd(geo.N, geo.J, np.abs(flux))
            Fmax, Ftot = np.max(Fabs), np.sum(Fabs)
            if gname == 'std':
                checks.append((f"{name}/{fname}", 'S2 |sum_i D_i| / sum_i sum_j |flux_ij| (no projection)', abs(D0.sum())/Ftot, 1e-14))
            base_div = div_metric(geo, D0, Brms)
            if np.linalg.norm(D0) < 1e-10*np.linalg.norm(Fabs):
                rows.append(dict(dist=name, field=fname, grad=gname, variant='none (D0 at round-off: faces close exactly, no projection needed)',
                                 cost_PPE=0.0, div=base_div, reduction=1.0, face=face_metric(geo, fld, QR, QL),
                                 incr_max=0.0, incr_p99=0.0, incr_where_B0=0, S3=0.0))
                print(f"    {fname:10s} {gname:3s}: base div {base_div:.3e}  (round-off: projection skipped)", flush=True)
                continue
            cex, marks, its = pcg(geo, D0)
            variants = [('none', np.zeros(geo.N), 0.0)]
            for K in (1, 2, 3, 5, 10):
                variants.append((f'jacobi K={K}', jacobi(geo, D0, K), float(K)))
            variants.append(('star', star(geo, D0), 2.0))
            variants.append(('stencil 1 pass', stencil(geo, D0, 1), st['cost']))
            variants.append(('stencil 2 passes', stencil(geo, D0, 2), 2*st['cost']))
            variants.append(('exact (PCG 1e-12)', cex, 1.0 + its))
            for vname, c, cost in variants:
                QRm, QLm = apply_correction(geo, QR, QL, c)
                Dm, _ = face_divergence(geo, QRm, QLm)
                r_lap = D0 - geo.Rc(c)
                s3 = np.max(np.abs(r_lap - Dm))/Fmax
                imax, i99, nzb = increment_ratio(geo, c, Bscale)
                dv = div_metric(geo, Dm, Brms)
                rows.append(dict(dist=name, field=fname, grad=gname, variant=vname, cost_PPE=round(cost, 1),
                                 div=dv, reduction=base_div/max(dv, 1e-300), face=face_metric(geo, fld, QRm, QLm),
                                 incr_max=imax, incr_p99=i99, incr_where_B0=nzb, S3=s3))
                checks.append((f"{name}/{fname}/{gname}/{vname}", 'S3 residual identity: max|r_Laplacian - D_direct| / max_i sum_j|flux_ij|', s3, 1e-13))
                if vname == 'exact (PCG 1e-12)':
                    checks.append((f"{name}/{fname}/{gname}", 'S4 exact mode: ||D_direct|| / ||D0|| (face sum recomputed)', np.linalg.norm(Dm)/np.linalg.norm(D0), 1e-10))
                    checks.append((f"{name}/{fname}/{gname}", 'S2 |sum_i D_i| / sum_i sum_j |flux_ij| (exact mode)', abs(Dm.sum())/Ftot, 1e-14))
                    rows[-1]['cg_marks'] = marks
            # S9 fixed point: starting from the exact solution, one step of each smoother must return it
            nc = np.linalg.norm(cex - cex.mean())
            for sname, fn in (('jacobi', lambda c0: jacobi(geo, D0, 1, c0=c0)), ('star', lambda c0: star(geo, D0, c0=c0)),
                              ('stencil', lambda c0: stencil(geo, D0, 1, c0=c0))):
                cn = fn(cex)
                checks.append((f"{name}/{fname}/{gname}", f'S9 exact solution is a fixed point of {sname}', np.linalg.norm(cn - cex)/nc, 1e-9))
            print(f"    {fname:10s} {gname:3s}: base div {base_div:.3e}  exact {rows[-1]['div']:.1e}  CG its to 1e-3/1e-12: {marks[1e-3]}/{marks[1e-12]}", flush=True)
    return rows

def self_check_jacobi_convergence(nngb, checks):
    """S5: weighted Jacobi converges to the exact solution (small glass, many sweeps)."""
    x = glass(12, iters=60)
    geo = Geometry(x, nngb, verbose=False)
    fld = make_field('grad_short', geo)
    grads, _ = limited_gradients(geo, fld['B'])
    QR, QL = face_states(geo, fld['B'], grads['std'])
    D0, _ = face_divergence(geo, QR, QL)
    cex, _, _ = pcg(geo, D0)
    hist = []
    c = jacobi(geo, D0, 6000, history=hist)
    h = np.array(hist)
    above = h[:-1] > 1e-12
    nonmono = np.mean(np.diff(h)[above] > 0.0) if above.any() else 0.0
    err = np.linalg.norm((c - c.mean()) - (cex - cex.mean()))/np.linalg.norm(cex - cex.mean())
    checks.append(('glass12/grad_short', 'S5 Jacobi residual after 6000 sweeps', h[-1], 1e-8))
    checks.append(('glass12/grad_short', 'S5 Jacobi residual: non-monotone sweeps while residual > 1e-12', nonmono, 0.0))
    checks.append(('glass12/grad_short', 'S5 Jacobi solution vs exact (mean removed)', err, 1e-6))

def self_check_linear(nngb, checks):
    """S8: GIZMO's MLS gradient reproduces a linear field exactly (non-periodic, interior particles)."""
    x = jittered(16, 0.3)
    geo = Geometry(x, nngb, periodic=False, verbose=False)
    G0 = rng.normal(size=(3, 3)); G0 -= np.trace(G0)/3*np.eye(3)
    B = np.array([0.2, -0.4, 0.7])[None, :] + x @ G0.T
    G, _ = mls_gradient(geo, B)
    inner = np.all((x > 2*geo.h[:, None]) & (x < 1 - 2*geo.h[:, None]), 1)
    err = np.max(np.abs(G[inner] - G0[None]))/np.max(np.abs(G0))
    checks.append(('jitter16-nonperiodic/linear', 'S8 MLS gradient of a linear field: max|G - G0|/|G0| (interior)', err, 1e-10))

def write_markdown(path, args, rows, checks, limiter):
    L = []
    L.append(f"# Phase-0 projection comparison: N = {args.n}^3, N_ngb = {args.ngb}, limiters {'ON' if limiter else 'OFF'}\n")
    L.append("Static test of the localized divergence projection (divfree/hydro/docs/LOCAL_PROJECTION_notes.pdf) with GIZMO's MFM faces, "
             "gradients and limiters replicated in Python (flat space). Generated by `phase0_projection.py`.\n")
    L.append("## Definitions\n")
    L.append("- `div` = rms_i( dx_i |D_i| / V_i ) / rms|B|, with D_i the face-sum divergence the hydro loop accumulates (dx_i = V_i^(1/3)).")
    L.append("- `reduction` = div(no projection) / div(variant), same gradient type.")
    L.append("- `face`: divergence-free fields (uniform, tg3d, abc): rms over faces of |B_face.n - B_exact(x_mid).n| / rms|B| (face reconstruction error); "
             "gradient fields (grad_short, grad_long: B = B0 + grad(phi)): rms|B_face.n - B0.n| / rms|grad(phi)(x_mid).n| (fraction of the gradient part left, 1 = untouched, 0 = removed); "
             "bump (GIZMO monopole-like): rms|B_face.n - B_bump(x_mid).n| / rms|B| (how much the faces were changed).")
    L.append("- `incr max / p99`: largest and 99th-percentile face increment |c_i| |D|^2 |A| relative to the local field scale max(|B|) over the stencil (what a cap would have to allow); "
             "`n(B=0)`: particles with a non-zero increment where the local field is exactly zero.")
    L.append("- `cost` in pair-pass equivalents (PPE): one pass over all pairs with face evaluation, i.e. one constrained-gradient-type iteration. "
             "Jacobi K sweeps = K (the first sweep is fused with the pass that computes D); star = 2; stencil = 1 + (intra-stencil pairs)/(pairs) per pass; exact = 1 + PCG iterations.")
    L.append("- gradient `std`: GIZMO MLS gradient with the per-row limiter; `tl`: traceless projection G - tr(G)/3 I with a single scalar limiter factor (keeps tr G = 0).\n")
    L.append("## Geometry\n")
    L.append("| dist | N | pairs | partners/particle | h/dx | cond median/max | RSPH faces | closure rms |sum_j A_ij|/dx^2 | closure / sum|A| | stencil cost (PPE/pass) |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for g in GEOSTATS:
        L.append(f"| {g['dist']} | {g['N']} | {g['pairs']} | {g['partners']:.1f} | {g['h_dx']:.2f} | {g['cond_med']:.2f} / {g['cond_max']:.2f} | {g['fallback']} | {g['closure_rms']:.2e} | {g['closure_rel']:.2e} | {g['stencil_cost']:.1f} |")
    L.append("")
    nfail = sum(1 for (_, _, v, t) in checks if not v <= t)
    L.append(f"## Self-checks: {len(checks)} evaluated, {nfail} FAIL\n")
    L.append("Each check is listed once per kind with the worst value over all cases (full list in the companion `phase0_selfchecks_*.txt`).\n")
    L.append("| check | worst value | threshold | status | worst case |")
    L.append("|---|---|---|---|---|")
    worst = {}
    for where, what, val, thr in checks:
        if what not in worst or val > worst[what][0]: worst[what] = (val, thr, where)
    for what in sorted(worst):
        v, t, w = worst[what]
        L.append(f"| {what} | {v:.2e} | {t:.0e} | {'PASS' if v <= t else 'FAIL'} | {w} |")
    L.append("")
    L.append("## Results\n")
    keyf = lambda r: (r['dist'], r['field'], r['grad'])
    seen = []
    for r in rows:
        if keyf(r) not in seen: seen.append(keyf(r))
    for (d, f, g) in seen:
        sel = [r for r in rows if keyf(r) == (d, f, g)]
        L.append(f"### {d} / {f} / gradient {g}\n")
        L.append("| variant | cost (PPE) | div | reduction | face | incr max | incr p99 | n(B=0) |")
        L.append("|---|---|---|---|---|---|---|---|")
        for r in sel:
            L.append(f"| {r['variant']} | {r['cost_PPE']} | {r['div']:.3e} | {r['reduction']:.3g} | {r['face']:.3e} | {r['incr_max']:.2e} | {r['incr_p99']:.2e} | {r['incr_where_B0']} |")
        m = [r for r in sel if 'cg_marks' in r]
        if m:
            mk = m[0]['cg_marks']
            L.append("\nPCG iterations to relative residual " + ", ".join(f"{k:.0e}: {mk[k]}" for k in sorted(mk, reverse=True)) + "\n")
        else:
            L.append("")
    open(path, 'w').write("\n".join(L) + "\n")

# ---------------------------------------------------------------- main
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=32)
    ap.add_argument('--ngb', type=int, default=32)
    ap.add_argument('--dists', nargs='+', default=['lattice', 'jitter', 'glass'])
    ap.add_argument('--fields', nargs='+', default=['uniform', 'tg3d', 'abc', 'grad_short', 'grad_long', 'bump'])
    ap.add_argument('--no-limiter', action='store_true')
    ap.add_argument('--skip-selfchecks', action='store_true')
    ap.add_argument('--out', default=os.path.dirname(os.path.abspath(__file__)))
    a = ap.parse_args()
    limiter = not a.no_limiter
    tag = f"n{a.n}_ngb{a.ngb}" + ("_nolim" if not limiter else "")
    print(f"# Phase-0 projection comparison: N = {a.n}^3, N_ngb = {a.ngb}, limiters {'ON' if limiter else 'OFF'}", flush=True)
    checks, rows = [], []
    if not a.skip_selfchecks:
        t0 = time.time(); self_check_linear(a.ngb, checks); self_check_jacobi_convergence(a.ngb, checks)
        print(f"  standalone self-checks done [{time.time()-t0:.1f}s]", flush=True)
    for dname in a.dists:
        print(f"== {dname}", flush=True)
        x = {'lattice': lambda: lattice(a.n), 'jitter': lambda: jittered(a.n, 0.3), 'glass': lambda: glass(a.n)}[dname]()
        geo = Geometry(x, a.ngb)
        rows += run_distribution(dname, geo, a.fields, limiter, checks)

    # ---------------- outputs
    csvp = os.path.join(a.out, f"phase0_results_{tag}.csv")
    keys = ['dist', 'field', 'grad', 'variant', 'cost_PPE', 'div', 'reduction', 'face', 'incr_max', 'incr_p99', 'incr_where_B0', 'S3']
    with open(csvp, 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=keys, extrasaction='ignore'); wr.writeheader(); wr.writerows(rows)
    chkp = os.path.join(a.out, f"phase0_selfchecks_{tag}.txt")
    nfail = 0
    with open(chkp, 'w') as fh:
        for where, what, val, thr in checks:
            ok = val <= thr
            nfail += (not ok)
            fh.write(f"{'PASS' if ok else 'FAIL'}  {what:68s} {val:10.2e} (<= {thr:.0e})  [{where}]\n")
    mdp = os.path.join(a.out, f"RESULTS_phase0_{tag}.md")
    write_markdown(mdp, a, rows, checks, limiter)
    print(f"# self-checks: {len(checks)} evaluated, {nfail} FAIL -> {chkp}")
    print(f"# report  -> {mdp}")
    print(f"# results -> {csvp}")
