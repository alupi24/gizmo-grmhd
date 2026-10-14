#ifndef GR_GR_H
#define GR_GR_H
//Enable the following in case you want to solve the kick/drift operations using an implicit solver (otherwise uses a simple predictor-corrector scheme)
#define IMPLICIT_SOLVER (1e-8)
#define IMPLICIT_MAX_NITER 100
//#define GR_TEST_PARTICLE
#define GR_EPS_TOLERANCE (1e-8) //Arbitrary dx used to compute metric derivatives

#if (GR_METRIC==2)
//#define GR_INCLUDE_GAS_CONTRIBUTION
#define GR_INCLUDE_PN_CORRECTIONS (4) //1: 1PN only, 2: 2PN only, 3: 1PN and 2PN, 4: 2.5PN only, 5: 1PN and 2.5PN, 6: 2PN and 2.5PN, 7: all on
//#ifdef BOX_PERIODIC
#define PRINT_OUTFLOW_Z (1480)
//#endif
//#define GR_INCLUDE_STELLAR_HARDENING
#endif

#define PARTICLE (25487837)
#define GR_GAMMA_LIMITER (1000)
#define GAMMA_LIMITER_RADIUS (0.0)
#define RIGHT_GAMMA_BOUNDARY (3)
#define LEFT_GAMMA_BOUNDARY (-1)
#define ZETA_CLEANING (1.0)
#ifndef GR_TIMESTEP_DV_FRAC
#define GR_TIMESTEP_DV_FRAC (0.01) /* time step: the change of the coordinate velocity per step, |a| dt, is limited to this fraction of the local signal speed (core/timestep.cc) */
#endif
#ifdef DIVBCLEANING_GR
/* Divergence cleaning with a genuine cleaning speed c_h ("Route B", see gr/ROUTE_B_cleaning_notes.pdf):
   c_h = v_ms (1+f)/(1+f v_ms^2) from the local magnetosonic speed (Fedrigo & Lupi 2025), capped at the speed of light;
   each particle uses its local value (SphP.CleaningSpeed), faces use the maximum over the pair.
   KAPPA_CLEANING is the dimensionless damping constant K: damping rate = alpha*K*c_h/Delta_x. */
#define KAPPA_CLEANING (0.1)  /* K of Fedrigo & Lupi (2025); the pre-Route-B scheme (GR_CLEANING_LEGACY) uses the same K with MaxSignalVel/2 in place of c_h */
#define SIGMAH (1.0)          /* multiplier of the hyperbolic source -alpha*c_h*D*divB */
#define VEL_CLEANING (1.0)    /* f: relativistic boost of the magnetosonic speed */
#ifndef GR_CLEANING_LEGACY
static inline double gr_cleaning_speed(double vms) {double ch = vms*(1.+VEL_CLEANING)/(1.+VEL_CLEANING*vms*vms); if(ch > 1.) {ch = 1.;} if(ch < 0.) {ch = 0.;} return ch;}
#else
/* GR_CLEANING_LEGACY: regression mode reproducing the pre-Route-B scheme with the same binary: light-speed cleaning states in the solver
   (c_h = 1 everywhere, so all Route B terms vanish), the old damping rate KAPPA_CLEANING_LEGACY*alpha*(MaxSignalVel/2)/Delta_x (as in Fedrigo & Lupi 2025),
   and the old relativistic boost of the fluid signal velocity for the time step. HLL and GForce are reproduced exactly; the old HLLC
   carried inconsistent c_h^2 factors and is not reproduced. */
#define KAPPA_CLEANING_LEGACY (0.1)
static inline double gr_cleaning_speed(double vms) {return 1.;}
#endif
#endif
#ifdef GR_BETACOOLING
#define GR_BETA (1.0)
#define HR_COOLING (0.1)
#endif
/* --- Metric coefficients --- */
struct gr_metric
{
       double gamma[6];           /*!< covariant metric coefficients */
       double ugamma[6];          /*!< contravariant metric coefficients */
       double beta[3];            /*!< contravariant metric shift */
       double alpha;              /*!< metric lapse term */
       double gdet;               /*!< square root of the metric determinant */
};


/* --- Primitive/conservative variables (for convenience) ---*/
struct gr_conservatives{
	double D;
	double S[3];
	double tau;
	double p;
#ifdef MAGNETIC
	double B[3];
#ifdef DIVBCLEANING_GR
        double phi;
#endif
#endif
};

struct gr_primitives{
	double rho0;
	double v[3];
	double u;
	double p;
#ifdef MAGNETIC
        double b[4];
        double bsq;
        double B[3];
#endif
};

/* NR kernel variables */
struct kernel_NR {
  double Bsq,QdotBsq,Qtsq,Qdotn,D,half_Bsq ;
} ;
/* --- */

/* --- main routines */
void compute_initial_accelerations(void);
void do_kick_general_relativity(int i,double dtime,double mass_old,double mass_new,int mode);
void do_drift_general_relativity(int p,double dt);
void do_predict_general_relativity(int p,double dt);
/* --- */

/* --- EOS functions */
double gr_get_primitive_density(double density,double velocity[],struct gr_metric *coeffs, int gr_limiter);
double gr_eos_info(double W, double vsq, double *dpdw, double *dpdvsq, struct kernel_NR *kernel);
double gr_pressure_rho0_u(double rho0, double u);
double gr_pressure_rho0_w(double rho0, double w);
double gr_get_u_from_p(double rho0, double p);
/* --- */

#ifdef MAGNETIC
double Get_Gas_BField_RestFrame(int i_particle_id, double* Bprim, double *vphys, double *bvec);
#endif

/* --- Velocity-related functions --- */
static inline void Vphys(double vel[],double beta[],double alpha) {int k; for(k=0;k<3;k++) vel[k]=(vel[k]+beta[k])/alpha;}
static inline void Vcoord(double vel[],double beta[],double alpha) {int k; for(k=0;k<3;k++) vel[k]=alpha*vel[k]-beta[k];}
long double vector_sqnorm(double gamma[],double vec[]);
double lorentz_factor(double gamma[],double vec[], int gr_limiter);
/* --- */

/* --- Primitive/conservative functions --- */
void gr_get_conservative_variables(struct gr_primitives *prims,struct gr_conservatives *cons,double Wlorentz,double gamma[],double gdet);
double gr_get_primitive_variables(struct gr_conservatives *cons,struct gr_primitives *prims,struct gr_metric *gr, int gr_limiter);
/* --- */


/* --- Metric-related functions --- */
void set_metric_coefficients(double* dp,double potential,struct gr_metric *gr);
//
void get_covariant_metric(double* dp,double *gcov);
void get_metric_derivatives(int p,double Kij[][3],double d_gij[][4][4],double* dalpha, double dbeta[][3]);

#if (GR_METRIC>0)
#define EHradius(M,a) (M*(1+sqrt(1-a[0]*a[0]-a[1]*a[1]-a[2]*a[2])))
#endif
#if (GR_METRIC<0)
void load_metric();
#elif (GR_METRIC==2)
void load_cbwaves();
#ifdef GR_INCLUDE_STELLAR_HARDENING
void hardening_params(); /* defined in gravity/stellar_hardening.h */
#endif
void gr_update_binary_position();
#ifdef GR_REMOVE_GRAD
void remove_radial_gradient_at_outer_edge(double *grad, int target);
#endif
#endif
#ifdef PARTICLE_EXCISION
void smooth_conservatives_inside_horizon(int i, double r, double r_EH, double M);
void check_boundary(int target);
#endif
#ifdef GR_BETACOOLING
void implicit_cooling(void);
#endif
/* --- */

/* --- Tensor-vector operations --- */
void unpack_4tensor(double A[],double Aunpacked[][4]);
void unpack_3tensor(double A[],double Aunpacked[][3]);
void raise_lower_indices(double g[],double v_in[], double v_out[]);
void rotate_tensor(double A[],double cosphi,double sinphi,double costh,double sinth,double Aprime[]);
double spatial_metric_determinant(double gamma[]);
/* --- */

#endif // GR_GR_H
