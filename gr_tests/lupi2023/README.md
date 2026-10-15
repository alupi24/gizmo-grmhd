# GR hydro tests of Lupi (2023, MNRAS 519, 1115) on the Kokkos port

Figures in the layout of the paper's own plotting scripts (`old_cluster/testIC/*`), produced by
`make_figures.py <run_base> [figure ...]` from the runs listed in its `RUNS` table. All runs use the
Kokkos branch `grmhd-kokkos` with `GENERAL_RELATIVITY=0` (HLL), the Kerr-Schild metric unless stated,
and the drift option `GR_DRIFT_HOLD_INTERNAL_ENERGY` (second-order geodesics in Kerr-Schild coordinates).
The exact special-relativistic Riemann solution comes from the `srrp` package used by the paper's script.

| figure | test (paper section) | setup | result |
|---|---|---|---|
| `fig_shocktube_marti1` | Marti & Muller problem 1, t=0.3 | 1D, 4000 cells on [0,2), MFM and MFV | L2(rho) 2.0e-2 (MFM), 1.3e-2 (MFV) |
| `fig_shocktube_marti2` | Marti & Muller problem 2, t=0.2 (inset: the thin shell) | as above | L2(rho) 1.4e-1 (MFM), 8.2e-2 (MFV); shell peak 10.5 vs 10.4 |
| `fig_shocktube_marti2vt` | problem 2 with v_t = 0.9, t=0.6 | as above | L2(rho) 0.25 (both): the shell between contact and shock is smeared at this resolution |
| `fig_shocktube_marti2vt16k` | the same at the paper's 16000 cells | 1.3e5 steps | L2(rho) 0.13, L2(v_x) 0.15, L2(v_y) 0.02 (both schemes) |
| `fig_orbits` | circular geodesics, a=0 r=10 and a=1 r=2 (Kerr-Schild) | 32 test particles, dt 0.0061, 4 orbits | max radial error 2.9e-6 (a=0), 1.7e-4 (a=1) |
| `fig_epicyclic` | epicyclic frequency kappa(r), a=-1, 0, 1 | circular orbits perturbed by 1e-5 in v, t=2500, dt 0.05 (a=1: 0.01) | median error 2e-4 to 4e-4 vs the Kerr value |
| `fig_apsidal` | apsidal precession of the r=90 eccentric orbit, a=0.1, -0.1, 0 | 32 particles, first 500 snapshots, dt 0.5 | 82.3 (dt 0.5), 82.36 (dt 0.1) deg/orbit vs 82.37 from the numerically integrated geodesic; paper marks 82.4 |
| `fig_freefall` | radial free fall from rest, Schwarzschild | 32 particles from r0 = 5..40, Boyer-Lindquist coordinates, excision at 1.05 horizon radii | |v_r| against the exact Schwarzschild solution |
| `fig_bondi` | Michel accretion, a=0, t=30 (paper's Fig. 12) | 2.7e5 cells (nx=72 in the paper's generator), excision at 1.05 horizon radii as the IC, MFM; top row Boyer-Lindquist with the original IC, bottom row Kerr-Schild with the IC converted to Kerr-Schild time | median |log rho - Michel| 0.005 dex in 4.2<r<30; accretion rate within 2.5 per cent of Michel for r>2.75 |
| `fig_polytrope`, `fig_polytrope_profiles` | (not in the paper) Newtonian n=1.5 polytrope, M=R=G=1, Lane-Emden profile, 18 t_dyn | 3.4e4 particles, positivity-limiter range factor 1.0, 0.5, 0.25 and step off | 1.0: outer layers drain to 0.37 of Lane-Emden, core contracts (rho_c 1.6); 0.5: same drift, slower (0.7 at 0.8 R); 0.25 and step off identical and clean (L1 0.007 inside 0.8 R) |
| `fig_tov_rhoc`, `fig_tov_profiles_MFM` | TOV star, K=1, Gamma=2, rho_c=0.129 (paper Figs 13-14) | 3.1e4 particles, tabulated metric (`GR_METRIC=-1`, areal radius), limiter range factor 0.25, 7.2 t_dyn | central density 0.95-1.02, interior within 2-3 per cent, no inward wave; the IC's surface overdensity at 0.85 R (1.2) is kept rather than smoothed |

Notes
- The shock tubes are run at 4000 cells (the transverse-velocity problem also at the paper's 16000): the GR time-step criteria
  (10 per cent of tau and of |S|, velocity change) make the strong-blast problems cost 1e5 steps; the cold side of problem 1
  (P=1e-6) needs `MinSizeTimestep` 1e-14 for the first steps.
- Test-particle runs use `GR_TEST_PARTICLE` (no hydro). For a=1 the orbits inside r=6 need a smaller step: at dt 0.05 the
  r=2.2 orbit decays at dr/dt = -4e-5 (the integration error, which scales as dt^2), so the a=1 panel combines a dt 0.01 run for
  r>6 with a dt 0.002 run of 32 orbits between r=1.6 and 6 (the two innermost, r<1.7, are captured at the marginally stable orbit).
- Bondi: the paper's IC is built in Schwarzschild (Boyer-Lindquist) time; for a Kerr-Schild run the conserved density and the
  coordinate velocity must be converted (D_KS = D_S (1 + 2M v_S/(r-2M)), v_KS = v_S/(1 + 2M v_S/(r-2M)); 17 per cent in D at r=5).
  With the code's default excision of 2.1 horizon radii (r=4.2) the flow inside r~10 sits 5-9 per cent below the Michel accretion
  rate; the old C code reproduces that to 1 per cent, so it is the boundary, not the port. With the excision at 1.05 (the IC's
  inner edge) the profile follows Michel to r~3. The old and the new drift give identical Bondi profiles.
- Positivity-preserving slope limiter (all hydro, not only GR): over the full neighbour range d_max its midpoint rule halves the
  gradient of every smooth stratified profile, and TOV stars and Newtonian polytropes develop an inward density wave (outer layers
  drain, the core contracts). The C code's `d_max *= 0.5` only slows it; the default is now `SLOPELIMITER_POSITIVITY_RANGE_FAC=0.25`,
  at which the step leaves smooth profiles alone. The relativistic magnetic rotor (`BOX_PERIODIC` is required for its IC) survives
  with 0.25, with the step fully off and with the GR face floor (`GR_FACE_FLOOR_FAC`, default 0.5) off; the floor is kept as a safety
  net and never binds in these tests. The rotor also runs without the GR velocity-change time-step criterion, which binds only
  below a fraction of ~0.03 (`GR_TIMESTEP_DV_FRAC`, default now 0.03).
- Limiter settings of the figures: the test-particle figures are independent of the limiter; the shock tubes were run with range
  factor 1.0, the Bondi runs with 0.5, the TOV figure with 0.25.
