# Phase 0: static test of the localized divergence projection — summary (13 Oct 2026)

Everything here is static (one evaluation of the face-sum divergence, no time evolution), flat space, with GIZMO's MFM
discretisation replicated in Python: smoothing lengths, NV_T, effective faces, MLS gradients, the row-wise B slope limiter
and the pairwise face limiter (references to the C source in the header of `phase0_projection.py`). N_ngb = 32.
Particle sets: perfect lattice, lattice jittered by ±0.3 dx, relaxed glass. Box [0,1)^3, periodic.

## Files

| file | content |
|---|---|
| `phase0_projection.py` | main comparison: all variants of the notes on 6 fields x 3 particle sets, with self-checks |
| `RESULTS_phase0_n32_ngb32.md` | full tables at 32^3 + geometry statistics + self-check summary |
| `phase0_selfchecks_n32_ngb32.txt` | all 547 individual self-checks |
| `phase0_results_n32_ngb32.csv` | raw numbers |
| `phase0_closure_split.py`, `RESULTS_closure_split_ngb32.md` | face-sum divergence split into face non-closure and reconstruction, 16^3-32^3 |
| `phase0_weights.py`, `RESULTS_weights_ngb32.md` | Tu et al. weighting versus two-point (TPFA) weighting, exact solves, 16^3-32^3 |
| `phase0_smoother_damping.py`, `RESULTS_smoother_damping_n24_ngb32.md` | Jacobi and stencil solve with different damping, at equal cost, both weightings |
| `phase0_tensor_faces.py`, `RESULTS_tensor_faces_ngb32.md` | face closure with `HYDRO_TENSOR_FACE_CORRECTIONS` (non-NGBITER) versus standard faces |

Addendum (tensor faces): the correction reduces the closure error about 3x (glass 3.7% -> 1.2%, jitter 19% -> 6%), not to zero;
on the jittered lattice it produces faces with A.dp < 0 (1068 of ~4e5 pairs at 24^3), which GIZMO would revert to RSPH faces.
The sign of the S3 term (code: S2 - S3.B; corrected second moment: S2 + S3.B) changes the closure by < 2%.
(The "first-moment check" column reads 1.0 on the lattice only because S1 = 0 there; it is 1e-14 elsewhere.)

Addendum (exact closure, `phase0_exact_closure.py`, `RESULTS_exact_closure_ngb32.md`): correcting the faces as
A~_ij = A_ij - |A_ij| (lambda_i - lambda_j), with sum_j |A_ij| (lambda_i - lambda_j) = sum_j A_ij (vector Laplacian solve),
closes them to 1e-15, keeps them antisymmetric, flips no face (k = 1 instead of |A| flips 1e4-1e5 faces), changes the faces
by 5% rms (glass) / 27% (jitter), leaves the first-order consistency essentially as it was (8.4% -> 8.7% on the glass),
and brings the face-sum divergence of tg3d down to the reconstruction level (glass 32^3: 0.112 -> 0.0098; jitter: 0.61 -> 0.028),
now converging with resolution. Jacobi from a cold start: 2 sweeps reduce the closure 10x, 5 sweeps 100x.

## Self-checks (547 evaluated, 0 FAIL)

| check | worst value | threshold |
|---|---|---|
| S0 smoothing length solve (N_eff = N_ngb) | 1.3e-15 | 1e-6 |
| S1 face antisymmetry A_ji = -A_ij (both sides computed) | 0 (bitwise) | 1e-15 |
| S2 sum_i D_i / total absolute flux (before and after projection) | 1.4e-17 | 1e-14 |
| S3 Laplacian residual = face-sum divergence recomputed with the modified faces | 1.5e-15 | 1e-13 |
| S4 exact mode: recomputed face-sum divergence / initial | 1.2e-12 | 1e-10 |
| S5 Jacobi: monotone above round-off, converges to the exact solution | 2e-12 | 1e-6 |
| S6 uniform field: D_i = B0 . sum_j A_ij | 3e-16 | 1e-12 |
| S7 traceless gradient stays traceless after the scalar limiter | 1e-14 | 1e-12 |
| S8 MLS gradient exact for linear fields (interior, non-periodic) | 3e-15 | 1e-10 |
| S9 exact solution is a fixed point of Jacobi, star and stencil solve | 2.4e-12 | 1e-9 |
| S10 stencil local solves accurate | 1e-15 | 1e-10 |

Normalisation note: checks are normalised by the total absolute face flux, not by D itself, because on the perfect lattice the
faces close exactly and D is at round-off for uniform and ABC fields (those cases are reported and the projection is skipped).

## Findings

**1. On irregular particle sets the face-sum divergence is dominated by face non-closure, not by the reconstruction.**
GIZMO's effective faces do not satisfy sum_j A_ij = 0: |sum_j A_ij| / sum_j |A_ij| = 3.7% on the glass, 19% on the jittered
lattice (round-off on the perfect lattice). The resulting term B_i . sum_j A_ij is zeroth order and does not converge:

| set | div(D) | closure part | reconstruction part, 16^3 -> 24^3 -> 32^3 (tg3d) |
|---|---|---|---|
| lattice | 3.5e-3 (32^3) | round-off | 2.1e-2 -> 7.4e-3 -> 3.5e-3 |
| glass | 0.11 | 0.11 (all N, all fields) | 2.2e-2 -> 1.3e-2 -> 9.0e-3 |
| jitter | 0.61 | 0.61 (all N, all fields) | 5.0e-2 -> 3.0e-2 -> 2.2e-2 |

(`div` = rms dx|D|/V / rms|B|.) Removing the closure term by projection distorts the face normal fluxes of a perfectly
uniform field by 2% of |B| (glass) and 11% (jitter), at every resolution; projecting only the reconstruction part changes
the faces by 0.16% (glass, tg3d, 32^3).

**2. The Tu et al. weighting does not converge to the correct projection; the two-point weighting does.**
For B = B0 + grad(phi), fraction of the gradient part left at the faces after the exact solve (0 = removed, 1 = untouched):

| set (32^3) | field | Tu, w = |D|^2|A|^2/2 | TPFA, w = |A|/r |
|---|---|---|---|
| lattice | grad, width 1.5 dx | 0.495 | 0.088 |
| lattice | grad, width 4 dx | 0.497 | 0.023 |
| glass (closure term removed) | width 1.5 dx | 0.681 | 0.129 |
| glass (closure term removed) | width 4 dx | 0.680 | 0.062 |

Both zero the face-sum divergence to round-off, and neither degrades the divergence-free tg3d field (face error unchanged).
The Tu increment ∝ r^2|A| (c_i - c_j) is not a consistent discrete gradient across neighbour shells at different r, so it
removes the divergence by distorting the field and leaves about half of the actual gradient part. With TPFA the face
increment is (2 c_i / r_ij) along n_ij: the same implementation (one scalar per particle), different weights.

**3. At equal cost, weighted Jacobi beats the stencil-local solve; omega ≈ 1 is best.**
Divergence relative to no projection, glass 24^3, TPFA weights (Tu weights give nearly identical numbers):

| field | Jacobi omega=1, 3 passes | Jacobi omega=1, 10 passes | stencil, 1 pass (14.7 passes) | stencil, 2 passes (29 passes) |
|---|---|---|---|---|
| tg3d (closure-dominated, short wavelength) | 0.031 | 0.003 | 0.121 | 0.021 |
| grad, width 4 dx | 0.431 | 0.326 | 0.428 | 0.351 |
| bump (GIZMO monopole-like) | 0.628 | 0.334 | 0.551 | 0.375 |

The stencil solve's cost is dominated by the neighbour-neighbour faces (13-16 pair passes per pass, measured); damping it
(omega_S < 1) only makes it worse. The star approximation (2 passes) performs like Jacobi with 2-3 passes.
Long-wavelength errors (the 4 dx gradient, the bump) decay slowly with any smoother, as expected for a Laplacian.

**4. Exact mode is affordable for validation, not obviously for production.**
Jacobi-preconditioned CG at 32^3: 11-42 iterations to 1e-3 and 108-149 to 1e-12 (count grows like N^(1/3)).

**5. The traceless gradient gives no measurable benefit.** Differences with the standard gradient are within a few percent,
and on the lattice it creates a small divergence for the ABC field where the standard gradient has none.

**6. A cap based on the local |B| is unusable.** For the bump (B = 0 outside), the correction potential extends outside the
source as for any Poisson problem, so face increments are non-zero where B = 0; relative to the local |B| they reach 1e3-1e5.
For smooth fields the increments stay below ~1.5 |B| (jitter) and ~0.25 |B| (glass).

## Consequences for the implementation plan

1. Fix the face closure first. GIZMO already has options for this (`KERNEL_CRK_FACES`, `HYDRO_TENSOR_FACE_CORRECTIONS`);
   they should be tested in this harness before anything else, because without closed faces any divergence-zeroing scheme
   (projection, constrained gradient, cleaning, Powell) acts mainly on a zeroth-order geometric artefact.
2. Use the TPFA weighting w_ij = |A_ij|/r_ij, increment (2 c_i/r_ij) n_ij, instead of the Tu et al. form.
3. Use weighted Jacobi (omega ≈ 0.8-1) as the smoother, a few passes per step; drop the stencil-local solve and the
   traceless gradient.
4. Long-wavelength residuals need either the exact mode (CG) or the hyperbolic cleaning: Phase 1 should keep Route B on
   and measure what remains dynamically.
5. Replace the cap of eq. (13) of the notes by one based on a larger-scale field estimate, or drop it.

Caveats: static test only, flat space, three idealised particle sets; the limiters use divB_prev = 0 (alim2 = a_limiter).
