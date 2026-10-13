# Route B: hyperbolic divergence cleaning with a true cleaning speed c_h in GR (revision 2)

Plain-text companion of ROUTE_B_cleaning_notes.pdf (typeset from ROUTE_B_cleaning_notes.tex,
which is the reference version). Derivation notes to be checked before any code change.
Conventions follow Fedrigo & Lupi (2025, A&A 700, A112; arXiv:2506.15775), hereafter FL25:
signature (-,+,+,+), n^mu = (1/alpha, -beta^i/alpha), n_mu = (-alpha,0,0,0),
gamma^{mu nu} = g^{mu nu} + n^mu n^nu (so gamma^{0 mu} = 0),
curly-B^i = sqrt(gamma) B^i, vtilde^i = alpha v^i - beta^i, D = sqrt(gamma) rho W,
K_ij the extrinsic curvature with K_ij = -alpha Gamma^0_ij (as in gr/metric.c),
K = gamma^{ij} K_ij, so that nabla_mu n^mu = -K, and a^mu = n^nu nabla_nu n^mu = D^mu ln(alpha)
(purely spatial).

Cleaning scalars: psi is Dedner's scalar (d_t psi = -c_h^2 div B in flat space);
psihat = psi / c_h is the Tricco-Price variable, which is FL25's psihat and the quantity
stored by the code (Phi). Revision 1 of these notes put the hat on the field multiplying the
cleaning tensor and called psi/c_h "Xi"; that was a mislabelling, corrected here. For c_h = 1
the two coincide.

## 1. The problem with eq. (27) of FL25

    nabla_mu ( F*^{mu nu} + g^{mu nu} psihat ) = sigma c_h psihat n^nu .

g^{mu nu} is covariantly constant and Lorentz covariant, so the (B_perp, psihat) subsystem
propagates at lambda = -beta^perp +/- alpha sqrt(gamma^{perp perp}): the coordinate speed
of light, for any c_h. c_h survives only through psihat = psi/c_h and the damping term.
The code's live DIVB_VARYING_CH branch implements these light-speed states, while the time
step and the damping use c_h = v_ms (1+f)/(1+f v_ms^2) < 1.

## 2. Modified covariant equation

    h^{mu nu} = gamma^{mu nu} - c_h^{-2} n^mu n^nu = g^{mu nu} + q n^mu n^nu,   q = 1 - c_h^{-2}

    nabla_mu ( F*^{mu nu} + h^{mu nu} psi ) = sigma psihat n^nu = (sigma/c_h) psi n^nu        (1)

Dedner's psi sits next to the tensor. For c_h = 1: h = g, psi = psihat, and (1) is FL25
eq. (27) including the damping. The right-hand side is normalised so that the damping rate
is alpha sigma c_h for every c_h (FL25's sigma c_h psihat n^nu = sigma psi n^nu would give
alpha sigma c_h^3 in the h-system).

    h^{00} = -c_h^{-2}/alpha^2,  h^{0j} = c_h^{-2} beta^j/alpha^2,  h^{ij} = gamma^{ij} - c_h^{-2} beta^i beta^j/alpha^2
    nabla_mu h^{mu nu} = n^nu (n^mu d_mu q) - q K n^nu + q a^nu                                  (2)

Flat-space check (c_h const): d_t B^j = -d_j psi, d_t psi = -c_h^2 d_i B^i: Dedner's
system, speeds +/- c_h. This fixes the field next to h as Dedner's psi.

## 3. 3+1 reduction

X = nabla_mu F*^{mu 0} = -d_i curlyB^i/(alpha sqrt(gamma)),
alpha sqrt(gamma) Y^j = d_t curlyB^j + d_i (vtilde^i curlyB^j - vtilde^j curlyB^i).

### 3.1 nu = 0: scalar equation (covariant form)

    (d_t - beta^i d_i) psi = - alpha c_h^2 d_i curlyB^i/sqrt(gamma) - alpha sigma c_h psi
                             + psi [ 2 (d_t - beta^i d_i) ln c_h - alpha (c_h^2 - 1) K ]          (3)

In terms of psihat = psi/c_h:

    (d_t - beta^i d_i) psihat = - alpha c_h d_i curlyB^i/sqrt(gamma) - alpha sigma c_h psihat
                                - alpha (c_h^2 - 1) K psihat + psihat (d_t - beta^i d_i) ln c_h  (4)

Fluid advection: FL25 advect psihat with the fluid (their eqs. 29-30, mass-weighted); the
-sqrt(gamma) beta^j v^i d_i psihat term of their induction equation is the fingerprint of
that choice (the pure covariant equation advects with -beta only). To be confirmed. Kept here.

Energy-consistent closure (recommended): the last term of (4) injects cleaning energy
where c_h decreases in time. A local c_h (max magnetosonic speed over the neighbours,
re-evaluated every step) is time dependent by construction and, being a max-filter,
discontinuous in time, so the term can be neither neglected nor evaluated. Tricco, Price &
Bate (2016) require instead that the hyperbolic exchange conserve int [B^2/2 + psihat^2/2]
for arbitrary c_h(x,t), which amounts to dropping that term. Proposed for implementation:

    d_t (D psihat) + d_i ( vtilde^i D psihat )
      = - alpha c_h D d_i curlyB^i/sqrt(gamma) - alpha sigma c_h D psihat - alpha (c_h^2 - 1) K D psihat   (5)

No derivative of c_h appears. For c_h = 1, (5) is FL25 eq. (30).

### 3.2 nu = j: induction equation

Eliminating d_t psi with the fluid-advected (3): the beta^i beta^j terms cancel, the Z-terms
cancel against psi n^j Z, the damping cancels against the right-hand side. Writing psi = c_h psihat:

    d_t curlyB^j + d_i ( vtilde^i curlyB^j - vtilde^j curlyB^i )
      = - sqrt(gamma) ( c_h^{-2} beta^j v^i + alpha gamma^{ij} ) d_i (c_h psihat) + beta^j d_i curlyB^i
        - sqrt(gamma) (1 - c_h^{-2}) c_h psihat gamma^{jk} d_k alpha                                  (6)

For c_h = 1, (6) is FL25 eq. (29). d_i (c_h psihat) needs no explicit grad c_h: it is built
from face values of psi = c_h psihat, as the phi gradient is today. Re-doing the elimination
with the closure (5) would generate a term ~ beta^j psihat (d_t - beta.d) ln c_h; it is dropped
with the other c_h-derivative terms.

## 4. Riemann states (Rusanov, |lambda| = alpha c_h sqrt(gamma^{perp perp}))

With psi_{L,R} = c_h^{L,R} psihat_{L,R} and the face c_h the larger of the two particle values:

    curlyBbar^perp = 1/2 (curlyB_L + curlyB_R) - 1/2 sqrt(gamma gamma^{perp perp}) (psi_R - psi_L) / c_h      (7a)
    psibar         = 1/2 (psi_L + psi_R) - 1/2 c_h (gamma gamma^{perp perp})^{-1/2} (curlyB_R - curlyB_L)    (7b)

Reduce to FL25 eq. (32) for c_h = 1. Powers: 1/c_h in (7a), c_h in (7b), not c_h^2 as in the
dead `#ifndef DIVB_VARYING_CH` branch. In the solver (primitive B): psi_{L,R} = c_h^{L,R} phi_{L,R};
`corr_p /= c_h`; `phi_normal_db *= c_h`; `phi_normal_mean` from psi_{L,R}.

## 5. Energy and momentum

FL25 eq. (34) unchanged in form (chain rule on the full B source Xi^k); tau is not conserved
by the cleaning, now as before. Cleaning energy E = int [B^2/2 + psihat^2/2]: with (4) the
exchange injects int psihat^2 d_t ln c_h; with the closure (5) the exchange is conservative for
any c_h(x,t) and the damping dissipative. Discretely the exchange is only approximately
antisymmetric when c_h^L != c_h^R (same approximation as Newtonian GIZMO). With psihat stored,
a jump of c_h between steps leaves psihat^2/2 unchanged.

## 6. Choice of c_h

With (5) no smoothness is required: FL25's local choice (max over the pair of
v_ms (1+f)/(1+f v_ms^2), re-evaluated every step) is admissible and keeps the time-step
hierarchy. Suggested additions: (i) a floor c_h -> 1 near the horizon, since the K term of (5)
is a growth term for c_h < 1 where K > 0 (ingoing Kerr-Schild) and the lapse term of (6) is
largest there; (ii) one function gr_cleaning_speed() used by solver, scalar source, damping
and Courant condition.

## 7. Mapping to the code (implemented 12 Oct 2026, for review)

| Item | File | Change |
|---|---|---|
| single definition of c_h | `gr_cleaning_speed()` in gr/gr.h (implemented) | replaces VEL_CLEANING boost in hydro/hydro_evaluate.h:273-277, `c_h2` in hydro/hydro_toplevel.c:803-807, damping in predict.c:550-557 |
| face states (7) | gr/gr_mhd_riemann_solver.h: HLL :118-131, GForce :433-446, HLLC :756-763 | psi = c_h phi per side; `corr_p /= c_h`, `phi_normal_db *= c_h`; delete the `ch2` logic; c_h^{L,R} via Riemann_vec |
| psihat source (5) | hydro/hydro_toplevel.c:892-900 | `tmp_ded*tmp_ded` -> `c_h` (linear) |
| K term (5) | hydro/hydro_toplevel.c cleaning block | `DtPhi -= alpha (c_h^2-1) K Phi`; K from get_metric_derivatives (gr/metric.c:365) |
| beta^j v^i term (6) | hydro/hydro_toplevel.c:871-876 | multiply by c_h^{-2} |
| lapse-gradient term (6) | hydro/hydro_toplevel.c cleaning block | `DtB[j] -= gdet (1-c_h^-2) c_h Phi gamma^jk dalpha_k`, through corr_tau/corr_S |
| damping | predict.c:556 | rate alpha sigma c_h, sigma = K/Delta x (restores the 1/h) |
| time step | hydro/hydro_evaluate.h / timestep.c | vsig = max(vsig_fluid, 2 c_h) |
| retire | allvars.h:946, gr/gr.h:24-26 | DIVB_VARYING_CH, VEL_CLEANING_MAX; keep f |

Regression: with c_h = 1 everywhere all new terms vanish; results must match the current code.

## Points to confirm

1. The beta^j v^i term of FL25 eq. (29) comes from eliminating d_t psihat with the
   fluid-advected equation (30).
2. Damping normalisation: sigma psihat n^nu keeps the rate alpha sigma c_h.
3. K sign convention (K_ij = -alpha Gamma^0_ij -> nabla.n = -K); K > 0 near the horizon in
   ingoing Kerr-Schild, so the K term grows psihat for c_h < 1: hence the floor c_h -> 1.
4. The closure drops every c_h-derivative term in both equations in exchange for energy
   consistency for arbitrary c_h(x,t); for constant c_h it coincides with the covariant system.
