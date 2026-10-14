/* gr/gr.h -- general-relativistic (magneto)hydrodynamics: compile-time constants, data structures and host prototypes.
 *
 * GR hydro: Lupi 2023, MNRAS, 519, 1115. GR MHD and the divergence cleaning: Fedrigo & Lupi 2025, A&A, 700, A112
 * (the "Route B" cleaning speed is documented in gr/ROUTE_B_cleaning_notes.pdf).
 *
 * This header is included from declarations/precompiler_logic.h, i.e. before the Kokkos macros and the typedefs are
 * available, so it holds only plain-C declarations: no KOKKOS_INLINE_FUNCTION bodies here. The device-callable
 * physics (metric, tensor algebra, EOS, conserved variables, con2prim) lives in gr/gr_functions.h and
 * gr/gr_con2prim_functions.h, included from declarations/cell_data.h; the host-only routines are in gr/gr_utils.cc.
 *
 * Port to the Kokkos code of the GR module of the C/C++ GIZMO branches (A. Lupi, G. Fedrigo). */
#ifndef GR_GR_H
#define GR_GR_H

/* kick/drift of the conserved variables: implicit solver tolerance and iteration cap (otherwise a predictor-corrector) */
#define IMPLICIT_SOLVER (1e-8)
#define IMPLICIT_MAX_NITER 100
#define GR_EPS_TOLERANCE (1e-8) /* finite-difference step of the metric derivatives */

#define PARTICLE (25487837)           /* ID of one particle for which the GR solver prints diagnostics (debug) */
#define GR_GAMMA_LIMITER (1000)       /* maximum Lorentz factor enforced inside GAMMA_LIMITER_RADIUS */
#define GAMMA_LIMITER_RADIUS (0.0)    /* radius (from the origin, or from the nearest sink for the binary metric) inside which the Lorentz factor is limited */
#define ZETA_CLEANING (1.0)
#ifndef GR_TIMESTEP_DV_FRAC
#define GR_TIMESTEP_DV_FRAC (0.01)    /* time step: the change of the coordinate velocity per step, |a| dt, is limited to this fraction of the local signal speed (core/timestep.cc) */
#endif

#ifdef DIVBCLEANING_GR
/* Divergence cleaning with a genuine cleaning speed c_h ("Route B", see gr/ROUTE_B_cleaning_notes.pdf):
   c_h = v_ms (1+f)/(1+f v_ms^2) from the local magnetosonic speed (Fedrigo & Lupi 2025), capped at the speed of light;
   each particle uses its local value (CellP.CleaningSpeed), faces use the maximum over the pair.
   KAPPA_CLEANING is the dimensionless damping constant K: damping rate = alpha*K*c_h/Delta_x. */
#define KAPPA_CLEANING (0.1)  /* K of Fedrigo & Lupi (2025); the pre-Route-B scheme (GR_CLEANING_LEGACY) uses the same K with MaxSignalVel/2 in place of c_h */
#define SIGMAH (1.0)          /* multiplier of the hyperbolic source -alpha*c_h*D*divB */
#define VEL_CLEANING (1.0)    /* f: relativistic boost of the magnetosonic speed */
#ifdef GR_CLEANING_LEGACY
/* GR_CLEANING_LEGACY: regression mode reproducing the pre-Route-B scheme with the same binary: light-speed cleaning states in the solver
   (c_h = 1 everywhere, so all Route B terms vanish), the damping rate KAPPA_CLEANING_LEGACY*alpha*(MaxSignalVel/2)/Delta_x (as in
   Fedrigo & Lupi 2025), and the old relativistic boost of the fluid signal velocity for the time step. HLL and GForce are reproduced exactly. */
#define KAPPA_CLEANING_LEGACY (0.1)
#endif
#endif

/* --- Metric coefficients at one point (3+1 split) --- */
struct gr_metric
{
    double gamma[6];           /*!< covariant spatial metric: xx, xy, xz, yy, yz, zz */
    double ugamma[6];          /*!< contravariant spatial metric, same ordering */
    double beta[3];            /*!< contravariant shift */
    double alpha;              /*!< lapse */
    double gdet;               /*!< sqrt of the spatial metric determinant */
    int limiter;               /*!< 1 where the Lorentz factor is limited to GR_GAMMA_LIMITER (inside GAMMA_LIMITER_RADIUS): set with the metric, so that no caller needs the position */
};

/* --- Conserved and primitive variables of one cell (the conserved ones are densitized, i.e. include sqrt(gamma)) --- */
struct gr_conservatives
{
    double D;
    double S[3];
    double tau;
    double p;                  /*!< total (gas + magnetic) pressure, carried along for the Riemann solvers */
#ifdef MAGNETIC
    double B[3];
#ifdef DIVBCLEANING_GR
    double phi;
#endif
#endif
};

struct gr_primitives
{
    double rho0;               /*!< rest-mass density */
    double v[3];               /*!< velocity measured by the Eulerian observer (contravariant) */
    double u;                  /*!< specific internal energy */
    double p;                  /*!< gas pressure */
#ifdef MAGNETIC
    double b[4];               /*!< rest-frame field (alpha*b^0, b_i) */
    double bsq;                /*!< b^mu b_mu */
    double B[3];               /*!< Eulerian field (primitive, not densitized) */
#endif
};

/* scalars shared by the con2prim Newton-Raphson iteration */
struct kernel_NR
{
    double Bsq, QdotBsq, Qtsq, Qdotn, D, half_Bsq;
};

#if (GR_METRIC>0)
#define EHradius(M,a) ((M)*(1+sqrt(1-(a)[0]*(a)[0]-(a)[1]*(a)[1]-(a)[2]*(a)[2]))) /* event-horizon radius of a Kerr hole of mass M and spin vector a (dimensionless) */
#endif

/* --- host routines (gr/gr_utils.cc) --- */
#ifdef PARTICLE_EXCISION
/* apply_excision() (gravity/analytic_gravity.h, called once per step) uses check_boundary() to remove the cells inside GR_EXCISION_RADIUS horizon radii */
void check_boundary(int target);
void smooth_conservatives_inside_horizon(int i, double r, double r_EH, double M);
#endif
void gr_con2prim_report_failure(const char *where, int code, const struct gr_conservatives *cons, const struct gr_primitives *prims, const struct gr_metric *gr); /* prints the state of a failed inversion (host) */

/* --- host routines of the conserved-variable integration (gr/gr_kickdrift.cc) --- */
void compute_initial_accelerations(void);
void do_kick_general_relativity(int i, double dtime, double mass_old, double mass_new, int mode);
void do_drift_general_relativity(int p, double dt);
void do_predict_general_relativity(int p, double dt);

#endif /* GR_GR_H */
