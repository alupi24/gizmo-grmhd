/* gr/gr_functions.h -- device-callable GR physics: tensor algebra, metric, EOS, conserved variables, rest-frame field.
 *
 * Every function here is a KOKKOS_INLINE_FUNCTION (host and device), takes plain arrays or the gr_* structs of
 * gr/gr.h, and reads simulation parameters only through All (the device sees the synced mirror). No particle
 * indices, no printf, no endrun: a caller that needs to report a failure does it on the host.
 *
 * Included from declarations/cell_data.h (after particle_data.h and the declaration of All), so that the cell
 * methods can use the GR density and sound speed; include it directly in any other header that needs it.
 *
 * Port of gr/metric.cc, gr/gr_eos.cc and parts of gr/gr_utils.cc of the C/C++ GIZMO branches (A. Lupi). */
#ifndef GR_GR_FUNCTIONS_H
#define GR_GR_FUNCTIONS_H
#ifdef GENERAL_RELATIVITY

#if (GR_METRIC==2) || (GR_METRIC<0)
#error "GR_METRIC=2 (superposed binary) and GR_METRIC<0 (tabulated metric) are not ported to the Kokkos code yet: use GR_METRIC=0 (flat) or 1 (Kerr)"
#endif

/* ------------------------------------------------------------------ tensor algebra (spatial metric packed as xx,xy,xz,yy,yz,zz) */
KOKKOS_INLINE_FUNCTION void raise_lower_indices(const double g[], const double vin[], double vout[])
{
    vout[0] = g[0]*vin[0] + g[1]*vin[1] + g[2]*vin[2];
    vout[1] = g[1]*vin[0] + g[3]*vin[1] + g[4]*vin[2];
    vout[2] = g[2]*vin[0] + g[4]*vin[1] + g[5]*vin[2];
}

KOKKOS_INLINE_FUNCTION double vector_sqnorm(const double gamma[], const double v_in[]) /* g_ij v^i v^j */
{
    double v_rl[3]; raise_lower_indices(gamma, v_in, v_rl);
    return v_rl[0]*v_in[0] + v_rl[1]*v_in[1] + v_rl[2]*v_in[2];
}

/* Lorentz factor of the contravariant velocity v_cont; with gr_limiter the velocity is limited IN PLACE to the maximum Lorentz
   factor GR_GAMMA_LIMITER. (The C code used long double here; device code uses double.) */
KOKKOS_INLINE_FUNCTION double lorentz_factor(const double gamma[], double v_cont[], int gr_limiter)
{
    double denom, R = vector_sqnorm(gamma, v_cont);
    if(gr_limiter == 1)
    {
        double RMAX = 1.0 - 1./((double)GR_GAMMA_LIMITER*(double)GR_GAMMA_LIMITER);
        double Rstar = 0.5*(RMAX + R - fabs(R - RMAX)), fac = sqrt(Rstar/(R + 1e-70));
        v_cont[0] *= fac; v_cont[1] *= fac; v_cont[2] *= fac;
        denom = 1. - Rstar;
    }
    else
    {
        denom = 1. - R;
        if(denom < MIN_REAL_NUMBER) {denom = MIN_REAL_NUMBER;}
    }
    return 1./sqrt(denom);
}

KOKKOS_INLINE_FUNCTION double spatial_metric_determinant(const double gamma[])
{
    return (-gamma[2]*gamma[2]*gamma[3] + 2.0*gamma[1]*gamma[2]*gamma[4] - gamma[0]*gamma[4]*gamma[4]
            - gamma[1]*gamma[1]*gamma[5] + gamma[0]*gamma[3]*gamma[5]);
}

KOKKOS_INLINE_FUNCTION void unpack_4tensor(const double Gin[], double Gout[][4])
{
    Gout[0][0]=Gin[0]; Gout[1][1]=Gin[4]; Gout[2][2]=Gin[7]; Gout[3][3]=Gin[9];
    Gout[0][1]=Gout[1][0]=Gin[1]; Gout[0][2]=Gout[2][0]=Gin[2]; Gout[0][3]=Gout[3][0]=Gin[3];
    Gout[1][2]=Gout[2][1]=Gin[5]; Gout[1][3]=Gout[3][1]=Gin[6]; Gout[2][3]=Gout[3][2]=Gin[8];
}

KOKKOS_INLINE_FUNCTION void unpack_3tensor(const double Gin[], double Gout[][3])
{
    Gout[0][0]=Gin[0]; Gout[1][1]=Gin[3]; Gout[2][2]=Gin[5];
    Gout[0][1]=Gout[1][0]=Gin[1]; Gout[0][2]=Gout[2][0]=Gin[2]; Gout[1][2]=Gout[2][1]=Gin[4];
}

KOKKOS_INLINE_FUNCTION void rotate_tensor(const double A[], double cosphi, double sinphi, double costh, double sinth, double Aprime[])
{
    double R[3][3] = {{costh*cosphi,-costh*sinphi,-sinth},{sinphi,cosphi,0},{sinth*cosphi,-sinth*sinphi,costh}};
    double tmp[3][3] = {{A[0],A[1],A[2]},{A[1],A[3],A[4]},{A[2],A[4],A[5]}}, rtmp[3][3];
    int i,j;
    for(i=0;i<3;i++) {for(j=0;j<3;j++) {rtmp[i][j] = tmp[i][0]*R[j][0] + tmp[i][1]*R[j][1] + tmp[i][2]*R[j][2];}}
    for(i=0;i<3;i++) {for(j=0;j<3;j++) {tmp[i][j] = R[i][0]*rtmp[0][j] + R[i][1]*rtmp[1][j] + R[i][2]*rtmp[2][j];}}
    Aprime[0]=tmp[0][0]; Aprime[1]=tmp[0][1]; Aprime[2]=tmp[0][2]; Aprime[3]=tmp[1][1]; Aprime[4]=tmp[1][2]; Aprime[5]=tmp[2][2];
}

/* coordinate velocity dx/dt <-> velocity measured by the Eulerian observer: v = (dx/dt + beta)/alpha */
KOKKOS_INLINE_FUNCTION void Vphys(double vel[], const double beta[], double alpha) {int k; for(k=0;k<3;k++) {vel[k] = (vel[k] + beta[k])/alpha;}}
KOKKOS_INLINE_FUNCTION void Vcoord(double vel[], const double beta[], double alpha) {int k; for(k=0;k<3;k++) {vel[k] = alpha*vel[k] - beta[k];}}

/* ------------------------------------------------------------------ metric */
/* 1 inside GAMMA_LIMITER_RADIUS of the origin (where the Lorentz factor is limited), measured from the raw coordinates as in the original code */
KOKKOS_INLINE_FUNCTION int gr_limiter_flag(const double pos[3])
{
    return (pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2] <= (double)GAMMA_LIMITER_RADIUS*(double)GAMMA_LIMITER_RADIUS) ? 1 : 0;
}

/* covariant 4-metric at pos, packed as tt, tx, ty, tz, xx, xy, xz, yy, yz, zz; the source is at the box centre if the box is periodic */
KOKKOS_INLINE_FUNCTION void get_covariant_metric(const double pos[3], double gcov[10])
{
    double dp[3], dp2[3], R2 = 0; int i;
    for(i=0;i<3;i++) {dp[i] = pos[i];}
#ifdef BOX_PERIODIC
    dp[0] -= boxHalf_X; dp[1] -= boxHalf_Y; dp[2] -= boxHalf_Z;
    NEAREST_XYZ(dp[0],dp[1],dp[2],1);
#endif
    for(i=0;i<3;i++) {dp2[i] = dp[i]*dp[i]; R2 += dp2[i];}
#if (GR_METRIC==0) /* Minkowski */
    gcov[0] = -1; gcov[1]=gcov[2]=gcov[3]=gcov[5]=gcov[6]=gcov[8]=0.; gcov[4]=gcov[7]=gcov[9]=1.;
#elif (GR_METRIC==1) /* Kerr: Boyer-Lindquist (USE_BOYERLINDQUIST) or Kerr-Schild coordinates */
    double M = All.MetricParams[0], a = All.MetricParams[1]*M, a2 = a*a;
    double r2 = (R2 - a2 + sqrt((R2-a2)*(R2-a2) + 4*a2*dp2[2]))/2.;
#ifndef USE_KERRSCHILD_COORDINATES
    double sinth2 = 1. - dp2[2]/r2, rho2 = r2 + a2*dp2[2]/r2, delta = r2 - 2*M*sqrt(r2) + a2, Att = 2.0*M*sqrt(r2)/rho2;
    gcov[0] = Att - 1.0;
    double gtp = -Att*a*sinth2, gpp = (r2 + a2 + Att*a2*sinth2)*sinth2;
    gcov[1] = -dp[1]/(dp2[0]+dp2[1])*gtp;
    gcov[2] =  dp[0]/(dp2[0]+dp2[1])*gtp;
    gcov[3] = 0.;
    double fac1 = r2/rho2/delta, fac2 = gpp/((dp2[0]+dp2[1])*(dp2[0]+dp2[1])), fac3 = dp2[2]/rho2/(r2*sinth2), fac4 = (1 - dp2[2]*(r2+a2)/r2/rho2);
    gcov[4] = fac1*dp2[0] + fac2*dp2[1] + fac3*dp2[0];
    gcov[7] = fac1*dp2[1] + fac2*dp2[0] + fac3*dp2[1];
    gcov[9] = (1-fac4)*(a2+r2)/delta + rho2/(r2*sinth2)*fac4*fac4;
    gcov[5] = dp[0]*dp[1]*(fac1 - fac2 + fac3);
    gcov[6] = (a2+r2)*dp[0]*dp[2]/rho2/delta - dp[0]*dp[2]/(r2*sinth2)*fac4;
    gcov[8] = (a2+r2)*dp[1]*dp[2]/rho2/delta - dp[1]*dp[2]/(r2*sinth2)*fac4;
#else
    double rho2 = r2*r2 + a2*dp2[2], r = sqrt(r2), Mfac = 2*M*r*r*r/rho2, denom = r2 + a2;
    double l_cov[4] = {1, (r*dp[0] + a*dp[1])/denom, (r*dp[1] - a*dp[0])/denom, dp[2]/r};
    gcov[0] = -1 + Mfac*l_cov[0]*l_cov[0];
    gcov[1] = Mfac*l_cov[0]*l_cov[1]; gcov[2] = Mfac*l_cov[0]*l_cov[2]; gcov[3] = Mfac*l_cov[0]*l_cov[3];
    gcov[4] = 1 + Mfac*l_cov[1]*l_cov[1]; gcov[5] = Mfac*l_cov[1]*l_cov[2]; gcov[6] = Mfac*l_cov[1]*l_cov[3];
    gcov[7] = 1 + Mfac*l_cov[2]*l_cov[2]; gcov[8] = Mfac*l_cov[2]*l_cov[3];
    gcov[9] = 1 + Mfac*l_cov[3]*l_cov[3];
#endif
#endif
}

KOKKOS_INLINE_FUNCTION void get_upper_spatial_metric(const double gcov[6], double gcon[6])
{
    double invdet = 1.e0/spatial_metric_determinant(gcov);
    gcon[0] = (-gcov[4]*gcov[4] + gcov[3]*gcov[5])*invdet;
    gcon[1] = ( gcov[2]*gcov[4] - gcov[1]*gcov[5])*invdet;
    gcon[2] = (-gcov[2]*gcov[3] + gcov[1]*gcov[4])*invdet;
    gcon[3] = (-gcov[2]*gcov[2] + gcov[0]*gcov[5])*invdet;
    gcon[4] = ( gcov[1]*gcov[2] - gcov[0]*gcov[4])*invdet;
    gcon[5] = (-gcov[1]*gcov[1] + gcov[0]*gcov[3])*invdet;
}

/* contravariant spatial metric, shift and lapse from the covariant 4-metric */
KOKKOS_INLINE_FUNCTION void get_extra_metric_coefficients(const double *gcov, double *gcon, double beta[3], double *alpha)
{
    get_upper_spatial_metric(gcov+4, gcon);
    double betacov[3] = {gcov[1], gcov[2], gcov[3]};
    raise_lower_indices(gcon, betacov, beta);
    double beta2 = beta[0]*gcov[1] + beta[1]*gcov[2] + beta[2]*gcov[3];
    *alpha = sqrt(beta2 - gcov[0]);
}

/* fills the metric structure at pos (potential: Newtonian correction of the lapse, only with self-gravity) */
KOKKOS_INLINE_FUNCTION void set_metric_coefficients(const double pos[3], double potential, struct gr_metric *gr)
{
    double gcov[10]; int k;
    get_covariant_metric(pos, gcov);
    for(k=0;k<6;k++) {gr->gamma[k] = gcov[k+4];}
    get_extra_metric_coefficients(gcov, gr->ugamma, gr->beta, &(gr->alpha));
    gr->gdet = sqrt(spatial_metric_determinant(gr->gamma));
#ifndef SELFGRAVITY_OFF
    gr->alpha = sqrt(gr->alpha*gr->alpha + 2*potential);
#endif
    gr->limiter = gr_limiter_flag(pos);
}

/* spatial derivatives of the 4-metric, lapse and shift at dp, by centred finite differences (Cartesian directions) */
KOKKOS_INLINE_FUNCTION void compute_metric_jacobian(const double dp[3], double deriv[][4][4], double dalpha[], double dbeta[][3])
{
    int dir,j,k; double dpH[3], dpL[3], gcovH[10], gcovL[10], alphaH, alphaL, betaH[3], betaL[3], gtmpH[4][4], gtmpL[4][4], gcon[6];
    for(dir=0;dir<3;dir++) {dpH[dir]=dpL[dir]=dp[dir];}
    for(dir=0;dir<3;dir++)
    {
        dpH[dir] += GR_EPS_TOLERANCE; dpL[dir] -= GR_EPS_TOLERANCE;
        get_covariant_metric(dpH, gcovH); get_covariant_metric(dpL, gcovL);
        get_extra_metric_coefficients(gcovH, gcon, betaH, &alphaH); get_extra_metric_coefficients(gcovL, gcon, betaL, &alphaL);
        unpack_4tensor(gcovH, gtmpH); unpack_4tensor(gcovL, gtmpL);
        for(j=0;j<4;j++)
        {
            deriv[dir][j][j] = 0.5*(gtmpH[j][j]-gtmpL[j][j])/GR_EPS_TOLERANCE;
            for(k=j+1;k<4;k++) {deriv[dir][j][k] = deriv[dir][k][j] = 0.5*(gtmpH[j][k]-gtmpL[j][k])/GR_EPS_TOLERANCE;}
        }
        dalpha[dir] = 0.5*(alphaH-alphaL)/GR_EPS_TOLERANCE;
        for(j=0;j<3;j++) {dbeta[dir][j] = 0.5*(betaH[j]-betaL[j])/GR_EPS_TOLERANCE;}
        dpH[dir]=dpL[dir]=dp[dir];
    }
}

/* extrinsic curvature K_ij = -alpha Gamma^0_ij and the metric derivatives at pos, given the metric there */
KOKKOS_INLINE_FUNCTION void get_metric_derivatives(const double pos[3], const struct gr_metric *gr, double Kij[][3], double d_gij[][4][4], double dalpha[], double dbeta[][3])
{
    int i,j,w; double invalpha2 = 1./(gr->alpha*gr->alpha), gmunu[4][4];
    gmunu[0][0] = -invalpha2;
    gmunu[0][1]=gmunu[1][0]=gr->beta[0]*invalpha2; gmunu[0][2]=gmunu[2][0]=gr->beta[1]*invalpha2; gmunu[0][3]=gmunu[3][0]=gr->beta[2]*invalpha2;
    gmunu[1][1]=gr->ugamma[0] - gr->beta[0]*gr->beta[0]*invalpha2;
    gmunu[2][2]=gr->ugamma[3] - gr->beta[1]*gr->beta[1]*invalpha2;
    gmunu[3][3]=gr->ugamma[5] - gr->beta[2]*gr->beta[2]*invalpha2;
    gmunu[1][2]=gmunu[2][1]=gr->ugamma[1] - gr->beta[0]*gr->beta[1]*invalpha2;
    gmunu[1][3]=gmunu[3][1]=gr->ugamma[2] - gr->beta[0]*gr->beta[2]*invalpha2;
    gmunu[2][3]=gmunu[3][2]=gr->ugamma[4] - gr->beta[1]*gr->beta[2]*invalpha2;
    compute_metric_jacobian(pos, d_gij, dalpha, dbeta);
    for(i=0;i<3;i++) {for(j=0;j<3;j++)
    {
        Kij[i][j] = 0.5*gmunu[0][0]*(d_gij[j][0][i+1] + d_gij[i][0][j+1]); /* the time derivative of the metric is zero */
        for(w=1;w<4;w++) {Kij[i][j] += 0.5*gmunu[0][w]*(d_gij[j][w][i+1] + d_gij[i][w][j+1] - d_gij[w-1][i+1][j+1]);}
        Kij[i][j] *= -gr->alpha;
    }}
}

/* ------------------------------------------------------------------ equation of state (ideal gas with GAMMA_DEFAULT, as in the con2prim) */
/* rest-mass density from the conserved density D = sqrt(gamma) rho0 W and the coordinate velocity */
KOKKOS_INLINE_FUNCTION double gr_get_primitive_density(double density, const double velocity[], const struct gr_metric *coeffs, int gr_limiter)
{
    double vel_p[3]; int k; for(k=0;k<3;k++) {vel_p[k] = velocity[k];}
    Vphys(vel_p, coeffs->beta, coeffs->alpha);
    return density/lorentz_factor(coeffs->gamma, vel_p, gr_limiter)/coeffs->gdet;
}

KOKKOS_INLINE_FUNCTION double gr_pressure_rho0_u(double rho0, double u)
{
#ifdef GR_TEST_PARTICLE
    return 0;
#else
    return (GAMMA_DEFAULT - 1.)*rho0*u;
#endif
}

KOKKOS_INLINE_FUNCTION double gr_pressure_rho0_w(double rho0, double w) /* from w = rho0 + rho0 u + p */
{
#ifdef GR_TEST_PARTICLE
    return 0;
#else
    return (GAMMA_DEFAULT - 1.)*(w - rho0)/GAMMA_DEFAULT;
#endif
}

KOKKOS_INLINE_FUNCTION double gr_get_u_from_p(double rho0, double p)
{
#ifdef GR_TEST_PARTICLE
    return 0;
#else
    return p/rho0/(GAMMA_DEFAULT - 1.);
#endif
}

/* pressure and its derivatives with respect to W = rho h gamma^2 and v^2, for the con2prim iteration */
KOKKOS_INLINE_FUNCTION double gr_eos_info(double W, double vsq, double *dpdw, double *dpdvsq, const struct kernel_NR *kernel)
{
    double ftmp = 1. - vsq, gtmp = sqrt(ftmp), gam_m1_o_gam = (GAMMA_DEFAULT - 1.)/GAMMA_DEFAULT;
    *dpdw = gam_m1_o_gam*ftmp;
    *dpdvsq = gam_m1_o_gam*(0.5*kernel->D/gtmp - W);
    return gam_m1_o_gam*(W*ftmp - kernel->D*gtmp);
}

/* relativistic sound speed of the ideal gas: c_s^2 = Gamma P / (rho0 h), rho0 h = rho0 (1+u) + P */
KOKKOS_INLINE_FUNCTION double gr_sound_speed(double gamma_eos, double pressure, double rho0, double u)
{
    double rho_h = rho0*(1. + u) + pressure;
    if(rho_h > 0 && pressure > 0) {return sqrt(gamma_eos*pressure/rho_h);}
    return 0;
}

/* ------------------------------------------------------------------ conserved variables from the primitives */
KOKKOS_INLINE_FUNCTION void gr_get_conservative_variables(const struct gr_primitives *prims, struct gr_conservatives *cons, double Wlorentz, const double gamma[], double gdet)
{
    double hstar = 1.0 + prims->u + prims->p/prims->rho0, Pstar = prims->p, v_cov[3];
    raise_lower_indices(gamma, prims->v, v_cov);
    cons->D = gdet*prims->rho0*Wlorentz;
    cons->S[0] = cons->S[1] = cons->S[2] = 0; cons->tau = 0;
#ifdef MAGNETIC
    cons->B[0] = gdet*prims->B[0]; cons->B[1] = gdet*prims->B[1]; cons->B[2] = gdet*prims->B[2];
    Pstar += prims->bsq/2; hstar += prims->bsq/prims->rho0;
    cons->S[0] = -gdet*prims->b[0]*prims->b[1];
    cons->S[1] = -gdet*prims->b[0]*prims->b[2];
    cons->S[2] = -gdet*prims->b[0]*prims->b[3];
    cons->tau  = -gdet*prims->b[0]*prims->b[0];
#endif
    double DhW = cons->D*Wlorentz*hstar;
    cons->S[0] += DhW*v_cov[0]; cons->S[1] += DhW*v_cov[1]; cons->S[2] += DhW*v_cov[2];
    cons->tau += DhW - gdet*Pstar - cons->D;
    if(cons->tau < 0) {cons->tau = MIN_REAL_NUMBER;} /* (the C code printed a warning here) */
    cons->p = Pstar;
}

#ifdef MAGNETIC
/* b^2 = b^mu b_mu of the Eulerian field B (primitive) and the Eulerian velocity vphys (contravariant); fills bvec with (alpha b^0, b_i).
   vphys may be limited in place (gr_limiter). */
KOKKOS_INLINE_FUNCTION double gr_bfield_restframe(const struct gr_metric *gr, const double B[3], double vphys[3], double bvec[4], int gr_limiter)
{
    double Wl = lorentz_factor(gr->gamma, vphys, gr_limiter), tmpvec[3], beta_cov[3], bmag = 0; int k;
    raise_lower_indices(gr->gamma, vphys, tmpvec);
    for(k=0;k<3;k++) {bmag += B[k]*tmpvec[k];}
    bvec[0] = Wl*bmag; /* b^0 always appears multiplied by alpha */
    for(k=0;k<3;k++) {tmpvec[k] = B[k]/Wl + Wl*bmag*(vphys[k] - gr->beta[k]/gr->alpha);}
    raise_lower_indices(gr->gamma, gr->beta, beta_cov);
    raise_lower_indices(gr->gamma, tmpvec, bvec+1);
    for(k=1;k<4;k++) {bvec[k] += beta_cov[k-1]/gr->alpha*bvec[0];}
    return bmag*bmag + vector_sqnorm(gr->gamma, B)/Wl/Wl;
}

/* relativistic Alfven speed v_A = sqrt(b^2/(rho0 h + b^2)), with rho0 h = rho0 (1+u) + P */
KOKKOS_INLINE_FUNCTION double gr_alfven_speed(double b2, double rho0, double u, double pressure)
{
    double hstar = rho0*(1. + u) + pressure + b2;
    if(b2 > 0 && hstar > 0) {return sqrt(b2/hstar);}
    return 0;
}

/* coordinate speed of the fast magnetosonic wave along the direction n (positive root), from the Eulerian-frame speed V^2 = c_s^2 (1 - v_A^2) + v_A^2 */
KOKKOS_INLINE_FUNCTION double gr_fast_wavespeed_coord(double cs, double vA, double Wfac, double vcoord_n, double beta_n, double alpha, double gamma_uii)
{
    double vA2 = vA*vA, V2 = cs*cs*(1 - vA2) + vA2;
    double a = (V2 - 1)*(Wfac*Wfac) - V2;
    double b = -2*(Wfac*Wfac*vcoord_n*(V2 - 1) + V2*beta_n);
    double c = (Wfac*Wfac*vcoord_n*vcoord_n*(V2 - 1) + V2*(alpha*gamma_uii - beta_n*beta_n));
    double disc = sqrt(b*b - 4*a*c);
    return (-b + disc)/(2*a);
}
#endif

#ifdef DIVBCLEANING_GR
/* Route B cleaning speed from the local magnetosonic speed (1 in the legacy mode) */
KOKKOS_INLINE_FUNCTION double gr_cleaning_speed(double vms)
{
#ifdef GR_CLEANING_LEGACY
    return 1.;
#else
    double ch = vms*(1. + VEL_CLEANING)/(1. + VEL_CLEANING*vms*vms);
    if(ch > 1.) {ch = 1.;} if(ch < 0.) {ch = 0.;}
    return ch;
#endif
}
#endif

#include "gr_con2prim_functions.h"

#endif /* GENERAL_RELATIVITY */
#endif /* GR_GR_FUNCTIONS_H */
