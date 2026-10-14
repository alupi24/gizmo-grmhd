/* gr/gr_riemann_functions.h -- the GR(M)HD Riemann solvers, device-callable.
 *
 * Hydro (no MAGNETIC): HLL (GENERAL_RELATIVITY=0) or the HLLC of Mignone & Bodo 2005 (GENERAL_RELATIVITY=1), with the
 * star state of Duez et al. (Lupi 2023). MHD: HLL (0) or the GForce centred scheme (1) of Fedrigo & Lupi 2025, with the
 * Dedner (B_n, psi) subsystem of the Route-B cleaning (gr/ROUTE_B_cleaning_notes.pdf). The MHD HLLC (2) is not ported.
 *
 * gr_riemann_solver() is the entry point called from hydro_core_meshless_functions.h in place of the Newtonian
 * Riemann_solver(): the states are the primitives of each side (rest-mass density, Eulerian velocity v^i, pressure,
 * primitive field B^i, psi = c_h*psihat), the metric at the face is Riemann_vec.gr and the face velocity Riemann_vec.v_frame
 * (Eulerian). The fluxes are returned in the lab frame, per unit face area, already in the coordinate frame of the moving
 * face: Fluxes.rho = F(D), Fluxes.v = F(S_i), Fluxes.p = F(tau), Fluxes.B = F(sqrt(gamma) B^i). A failed solve leaves
 * P_M < 0, which triggers the lower-order retries of the caller; nothing here prints or stops.
 *
 * Port of gr/gr_riemann_solver.h and gr/gr_mhd_riemann_solver.h of the C/C++ GIZMO branches (A. Lupi, G. Fedrigo). */
#ifndef GR_RIEMANN_FUNCTIONS_H
#define GR_RIEMANN_FUNCTIONS_H
#ifdef GENERAL_RELATIVITY

#define GR_DEDNER_IMPLICIT_LIMITER (0.75) /* the psi jump may change the face-normal field by at most this fraction of it */

/* components of a symmetric tensor in the rotated basis whose rows are R: g' = R g R^T (R orthonormal) */
KOKKOS_INLINE_FUNCTION SymmetricTensor2<double> gr_rotate_tensor(const SymmetricTensor2<double> &g, const Mat3<double> &R)
{
    SymmetricTensor2<double> gp;
    for(int a = 0; a < 3; a++) {for(int b = a; b < 3; b++) {gp[a][b] = dot(R[a], g.matvec(R[b]));}}
    return gp;
}

/* the left/right state of a face: conserved variables and the quantities of the primitives the fluxes need */
struct gr_face_state
{
    struct gr_conservatives cons;
    Vec3<double> v;        /*!< Eulerian velocity (possibly limited) */
    Vec3<double> vcoord;   /*!< coordinate velocity alpha v - beta */
    double W;              /*!< Lorentz factor */
#ifdef MAGNETIC
    double ab0;             /*!< alpha b^0 */
    Vec3<double> b;        /*!< covariant rest-frame field b_i */
    double bsq;
    double vA;             /*!< relativistic Alfven speed */
#endif
};

/* ------------------------------------------------------------------ wave speeds */
#ifndef MAGNETIC
/* fastest coordinate signal speeds along n of a state with Eulerian velocity v (hydro): lam[0] = lambda_+, lam[1] = lambda_- */
KOKKOS_INLINE_FUNCTION void gr_lambda_hydro(const Vec3<double> &v, double cs, double g_nn, double beta_n, double alpha, const Vec3<double> &n_unit, double Wfac, double lam[2], int gr_limiter)
{
    double vtilde_n = dot(v, n_unit), cs2 = cs*cs;
    double v2max = 1. - MIN_REAL_NUMBER; if(gr_limiter) {v2max = 1. - 1./GR_GAMMA_LIMITER/GR_GAMMA_LIMITER;}
    double v2 = 1. - 1./Wfac/Wfac; if(v2 < 0) {v2 = 0.;}
    double discr = sqrt((1-v2)*(g_nn*(1-v2*cs2) - vtilde_n*vtilde_n*(1-cs2)));
    lam[0] = alpha*DMIN((vtilde_n*(1-cs2) + cs*discr)/(1-v2*cs2),  v2max) - beta_n;
    lam[1] = alpha*DMAX((vtilde_n*(1-cs2) - cs*discr)/(1-v2*cs2), -v2max) - beta_n;
}
#else
/* fastest coordinate signal speeds of a state in the face-aligned frame (x = normal), with the fast magnetosonic speed */
KOKKOS_INLINE_FUNCTION void gr_lambda_mhd(const Vec3<double> &v, double cs, double va, double g_nn, double beta_n, double alpha, double Wfac, double lam[2], int gr_limiter)
{
    double cs2 = cs*cs*(1-va*va) + va*va;
    double v2max = 1. - MIN_REAL_NUMBER; if(gr_limiter == 1) {v2max = 1. - 1./GR_GAMMA_LIMITER/GR_GAMMA_LIMITER;}
    double vtilde_n = alpha*v[0] - beta_n, W2 = Wfac*Wfac;
    double a = W2*(cs2-1) - cs2;
    double b = -2*(W2*(cs2-1)*vtilde_n + cs2*beta_n);
    double c = cs2*(alpha*alpha*g_nn - beta_n*beta_n) + W2*vtilde_n*vtilde_n*(cs2-1);
    double discr = sqrt(b*b - 4*a*c);
    lam[0] = DMIN((-b - discr)/(2*a),  v2max);
    lam[1] = DMAX((-b + discr)/(2*a), -v2max);
}
#endif

/* ------------------------------------------------------------------ face states */
/* conserved variables (and the rest-frame field) of one side from its primitives; v is limited in place where required */
KOKKOS_INLINE_FUNCTION void gr_face_state_from_primitives(const struct Conserved_var_Riemann &side, const struct gr_metric *gr, int gr_limiter, struct gr_face_state *st
#ifdef MAGNETIC
                                                          , double B_normal
#endif
                                                          )
{
    struct gr_primitives prims;
    prims.rho0 = side.rho; prims.v = side.v; prims.u = side.u; prims.p = side.p;
    st->W = lorentz_factor(gr->gamma, prims.v, gr_limiter);
#ifdef MAGNETIC
    prims.B = side.B; prims.B[0] = B_normal; /* the face-normal field is the common corrected value (x is the normal in the rotated frame) */
    prims.bsq = gr_bfield_restframe(gr, prims.B, prims.v, &prims.ab0, &prims.b, gr_limiter);
    st->ab0 = prims.ab0; st->b = prims.b; st->bsq = prims.bsq;
    st->vA = sqrt(prims.bsq/(prims.rho0*(1+prims.u) + (prims.p + prims.bsq)));
#endif
    gr_get_conservative_variables(&prims, &(st->cons), st->W, gr);
    st->v = prims.v;
    st->vcoord = Vcoord(prims.v, gr);
}

#ifdef MAGNETIC
/* physical flux of the conserved variables through a face with normal x (rotated frame): F = U vx + pressure and field terms */
KOKKOS_INLINE_FUNCTION void gr_mhd_physical_flux(const struct gr_face_state *st, const struct gr_metric *gr, double beta_n, struct gr_conservatives *flux)
{
    const struct gr_conservatives &U = st->cons;
    double vx = st->vcoord[0], Bx_over_W = U.B[0]/st->W, agp = gr->alpha*gr->gdet*U.p;
    flux->D   = U.D*vx;
    flux->tau = U.tau*vx + gr->gdet*U.p*(vx + beta_n) - gr->alpha*st->ab0*Bx_over_W;
    flux->S   = vx*U.S - (gr->alpha*Bx_over_W)*st->b;
    flux->S[0] += agp;
    flux->B   = vx*U.B - U.B[0]*st->vcoord;
    flux->B[0] = 0; /* the normal field is continuous: no flux */
    flux->p = U.p;
#ifdef DIVBCLEANING_GR
    flux->phi = 0; /* the psi flux is mass-based in the caller (Fluxes.rho*psihat): the solver's value is never read */
#endif
}

/* the face-normal field shared by both states, with the Dedner psi-jump correction of the Route-B cleaning; returns the corrected B_n
   and fills the cleaning outputs of Riemann_out */
KOKKOS_INLINE_FUNCTION double gr_mhd_face_normal_field(const struct Input_vec_Riemann &Riemann_vec, double gnn, struct Riemann_outputs *Riemann_out)
{
    double BxL = Riemann_vec.L.B[0], BxR = Riemann_vec.R.B[0];
    double B_normal_corrected = 0.5*(BxL + BxR);
#ifdef DIVBCLEANING_GR
    /* (B_perp, psi) subsystem, Route B (gr/ROUTE_B_cleaning_notes.pdf, eqs. 13-14): Rusanov states with characteristic speed alpha*c_h*sqrt(gamma^nn).
       L.phi/R.phi carry Dedner's psi = c_h*psihat of each side (multiplied in the meshless core); c_h is the face value (max over the pair).
       Powers of c_h: 1/c_h on the psi jump in B*, c_h on the B jump in psi*. For c_h = 1 this is Fedrigo & Lupi (2025) eq. (32). */
    double c_h = Riemann_vec.c_h, phiL = Riemann_vec.L.phi, phiR = Riemann_vec.R.phi, corr_norm = 1.0, corr_p = 0.;
    if(c_h > 0) {corr_p = 0.5*(phiL - phiR)*sqrt(gnn)/c_h;}
    double corr_p_abs = fabs(corr_p), corr_b_abs = GR_DEDNER_IMPLICIT_LIMITER*fabs(B_normal_corrected);
    if(corr_p_abs > corr_b_abs) {corr_norm *= corr_b_abs/corr_p_abs;}
    B_normal_corrected += corr_norm*corr_p;
    Riemann_out->phi_normal_mean = 0.5*corr_norm*(phiL + phiR);
    Riemann_out->phi_normal_db = corr_norm*0.5*(BxL - BxR)*c_h/sqrt(gnn);
#endif
    Riemann_out->B_normal_corrected = Riemann_vec.gr.gdet*B_normal_corrected;
    return B_normal_corrected;
}

/* writes the flux of the conserved variables in the frame of the face (moving at v_frame_line) into Riemann_out */
KOKKOS_INLINE_FUNCTION void gr_set_output_flux(const struct gr_conservatives *F, const struct gr_conservatives *U, double v_frame_line, struct Riemann_outputs *Riemann_out)
{
    Riemann_out->Fluxes.rho = F->D - U->D*v_frame_line;
    Riemann_out->Fluxes.p   = F->tau - U->tau*v_frame_line;
    Riemann_out->Fluxes.v   = F->S - v_frame_line*U->S;
    Riemann_out->Fluxes.B   = F->B - v_frame_line*U->B;
}

/* HLL solver (GENERAL_RELATIVITY=0), in the face-aligned frame */
KOKKOS_INLINE_FUNCTION void gr_mhd_riemann_hll(const struct Input_vec_Riemann &Riemann_vec, struct Riemann_outputs *Riemann_out, int gr_limiter)
{
    const struct gr_metric *gr = &(Riemann_vec.gr);
    double beta_n = gr->beta[0], gnn = gr->ugamma[0][0];
    double v_frame_line = gr->alpha*Riemann_vec.v_frame[0] - beta_n;
    double B_normal = gr_mhd_face_normal_field(Riemann_vec, gnn, Riemann_out);
    struct gr_face_state L, R;
    gr_face_state_from_primitives(Riemann_vec.L, gr, gr_limiter, &L, B_normal);
    gr_face_state_from_primitives(Riemann_vec.R, gr, gr_limiter, &R, B_normal);
    double lam_l[2], lam_r[2];
    gr_lambda_mhd(L.v, Riemann_vec.L.cs, L.vA, gnn, beta_n, gr->alpha, L.W, lam_l, gr_limiter);
    gr_lambda_mhd(R.v, Riemann_vec.R.cs, R.vA, gnn, beta_n, gr->alpha, R.W, lam_r, gr_limiter);
    double S_L = DMIN(lam_r[1], lam_l[1]), S_R = DMAX(lam_r[0], lam_l[0]);
    Riemann_out->S_M = 0.5*(S_L + S_R);
    Riemann_out->P_M = 0.5*(L.cons.p + R.cons.p);
    struct gr_conservatives fluxL, fluxR;
    gr_mhd_physical_flux(&L, gr, beta_n, &fluxL);
    gr_mhd_physical_flux(&R, gr, beta_n, &fluxR);
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    if(S_L > v_frame_line) {gr_set_output_flux(&fluxL, &(L.cons), v_frame_line, Riemann_out); return;}
    if(S_R < v_frame_line) {gr_set_output_flux(&fluxR, &(R.cons), v_frame_line, Riemann_out); return;}
#endif
    /* HLL intermediate state and flux */
    double facLR = 1./(S_R - S_L);
    struct gr_conservatives U_hll, F_hll;
    U_hll.D   = (S_R*R.cons.D   - S_L*L.cons.D   + fluxL.D   - fluxR.D)*facLR;
    U_hll.tau = (S_R*R.cons.tau - S_L*L.cons.tau + fluxL.tau - fluxR.tau)*facLR;
    U_hll.S   = facLR*(S_R*R.cons.S - S_L*L.cons.S + fluxL.S - fluxR.S);
    U_hll.B   = facLR*(S_R*R.cons.B - S_L*L.cons.B + fluxL.B - fluxR.B);
    F_hll.D   = (S_R*fluxL.D   - S_L*fluxR.D   + S_L*S_R*(R.cons.D   - L.cons.D))*facLR;
    F_hll.tau = (S_R*fluxL.tau - S_L*fluxR.tau + S_L*S_R*(R.cons.tau - L.cons.tau))*facLR;
    F_hll.S   = facLR*(S_R*fluxL.S - S_L*fluxR.S + (S_L*S_R)*(R.cons.S - L.cons.S));
    F_hll.B   = facLR*(S_R*fluxL.B - S_L*fluxR.B + (S_L*S_R)*(R.cons.B - L.cons.B));
#ifdef HYDRO_MESHLESS_FINITE_MASS
    /* the HLL equivalent of the zero-mass-flux face velocity */
    v_frame_line = F_hll.D/U_hll.D;
    Riemann_out->S_M = v_frame_line;
#endif
    gr_set_output_flux(&F_hll, &U_hll, v_frame_line, Riemann_out);
#ifdef HYDRO_MESHLESS_FINITE_MASS
    Riemann_out->Fluxes.rho = 0;
#endif
}

/* GForce solver (GENERAL_RELATIVITY=1): weighted Lax-Friedrichs / Lax-Wendroff fluxes, in the face-aligned frame */
KOKKOS_INLINE_FUNCTION void gr_mhd_riemann_gforce(const struct Input_vec_Riemann &Riemann_vec, struct Riemann_outputs *Riemann_out, int gr_limiter)
{
    const struct gr_metric *gr = &(Riemann_vec.gr);
    double beta_n = gr->beta[0], gnn = gr->ugamma[0][0];
    double v_frame_line = gr->alpha*Riemann_vec.v_frame[0] - beta_n;
    double omega_g = 1./(1. + All.CourantFac);
    double B_normal = gr_mhd_face_normal_field(Riemann_vec, gnn, Riemann_out);
    struct gr_face_state L, R;
    gr_face_state_from_primitives(Riemann_vec.L, gr, gr_limiter, &L, B_normal);
    gr_face_state_from_primitives(Riemann_vec.R, gr, gr_limiter, &R, B_normal);
    double lam_l[2], lam_r[2];
    gr_lambda_mhd(L.v, Riemann_vec.L.cs, L.vA, gnn, beta_n, gr->alpha, L.W, lam_l, gr_limiter);
    gr_lambda_mhd(R.v, Riemann_vec.R.cs, R.vA, gnn, beta_n, gr->alpha, R.W, lam_r, gr_limiter);
    double S_L = DMIN(lam_r[1], lam_l[1]), S_R = DMAX(lam_r[0], lam_l[0]);
    double S_MAX = 1./DMAX(fabs(S_L), fabs(S_R)); /* the (inverse) maximum signal speed: dt/dx of the centred schemes */
    Riemann_out->S_M = 0.5*(S_L + S_R);
    Riemann_out->P_M = 0.5*(L.cons.p + R.cons.p);
    struct gr_conservatives fluxL, fluxR;
    gr_mhd_physical_flux(&L, gr, beta_n, &fluxL);
    gr_mhd_physical_flux(&R, gr, beta_n, &fluxR);
    /* Lax-Friedrichs flux */
    struct gr_conservatives fluxLF;
    fluxLF.D   = 0.5*(fluxR.D + fluxL.D)     + 0.5*(L.cons.D - R.cons.D)/S_MAX;
    fluxLF.tau = 0.5*(fluxR.tau + fluxL.tau) + 0.5*(L.cons.tau - R.cons.tau)/S_MAX;
    fluxLF.S   = 0.5*(fluxR.S + fluxL.S)     + (0.5/S_MAX)*(L.cons.S - R.cons.S);
    fluxLF.B   = 0.5*(fluxR.B + fluxL.B)     + (0.5/S_MAX)*(L.cons.B - R.cons.B);
    fluxLF.B[0] = 0;
    /* Lax-Wendroff state, inverted to primitives for its flux */
    struct gr_face_state LW;
    LW.cons.D   = 0.5*(L.cons.D + R.cons.D)     + 0.5*S_MAX*(fluxL.D - fluxR.D);
    LW.cons.tau = 0.5*(L.cons.tau + R.cons.tau) + 0.5*S_MAX*(fluxL.tau - fluxR.tau);
    LW.cons.S   = 0.5*(L.cons.S + R.cons.S)     + (0.5*S_MAX)*(fluxL.S - fluxR.S);
    LW.cons.B   = 0.5*(L.cons.B + R.cons.B)     + (0.5*S_MAX)*(fluxL.B - fluxR.B);
    LW.cons.B[0] = gr->gdet*B_normal;
    LW.cons.p = 0;
#ifdef DIVBCLEANING_GR
    LW.cons.phi = 0;
#endif
    struct gr_primitives prims; /* the right state as the initial guess of the inversion */
    prims.rho0 = Riemann_vec.R.rho; prims.v = R.v; prims.u = Riemann_vec.R.u; prims.p = Riemann_vec.R.p;
    prims.B = R.cons.B/gr->gdet; prims.ab0 = R.ab0; prims.b = R.b; prims.bsq = R.bsq;
    double res_con2prim = gr_get_primitive_variables(&(LW.cons), &prims, gr, gr_limiter);
    if(res_con2prim < 0) {Riemann_out->P_M = -1; return;} /* the caller retries with lower-order states */
    LW.vcoord = Vcoord(prims.v, gr);
    LW.W = lorentz_factor(gr->gamma, prims.v, gr_limiter);
    LW.bsq = gr_bfield_restframe(gr, prims.B, prims.v, &(LW.ab0), &(LW.b), gr_limiter);
    LW.v = prims.v;
    LW.cons.p = prims.p + 0.5*LW.bsq; /* total pressure of the Lax-Wendroff state */
    struct gr_conservatives fluxLW;
    gr_mhd_physical_flux(&LW, gr, beta_n, &fluxLW);
    /* GForce = omega LW + (1-omega) LF, in the frame of the face */
    struct gr_conservatives F;
    F.D   = omega_g*fluxLW.D   + (1-omega_g)*fluxLF.D;
    F.tau = omega_g*fluxLW.tau + (1-omega_g)*fluxLF.tau;
    F.S   = omega_g*fluxLW.S   + (1-omega_g)*fluxLF.S;
    F.B   = omega_g*fluxLW.B   + (1-omega_g)*fluxLF.B;
#ifdef HYDRO_MESHLESS_FINITE_MASS
    v_frame_line = F.D/LW.cons.D; /* zero-mass-flux face velocity */
    Riemann_out->S_M = v_frame_line;
#endif
    gr_set_output_flux(&F, &(LW.cons), v_frame_line, Riemann_out);
#ifdef HYDRO_MESHLESS_FINITE_MASS
    Riemann_out->Fluxes.rho = 0;
#endif
}

#else /* ------------------------------------------------------------------ hydro solvers (no MAGNETIC) */

/* star state (contact wave speed and pressure) of Duez et al. / Mignone & Bodo 2005; returns 0 on success, 1 on failure */
KOKKOS_INLINE_FUNCTION int gr_starstate(const struct gr_conservatives &cons_L, const struct gr_conservatives &cons_R, const Vec3<double> &n_unit, double vline_L, double vline_R,
                                        const struct gr_metric *gr, double beta_n, double g_nn, double S_L, double S_R, double *S_M, double *P_M, int gr_limiter)
{
    double sXR = dot(gr->ugamma.matvec(cons_R.S), n_unit), sXL = dot(gr->ugamma.matvec(cons_L.S), n_unit); /* S^n */
    /* the coefficients of the quadratic equation for lambda* (Duez et al. eqs. 23-26) */
    double AR = (cons_R.tau + cons_R.D)*(S_R + beta_n) - gr->alpha*sXR;
    double AL = (cons_L.tau + cons_L.D)*(S_L + beta_n) - gr->alpha*sXL;
    double CR = gr->alpha*g_nn, CL = CR;
    double BR = CR*gr->gdet*cons_R.p + sXR*(vline_R - S_R);
    double BL = CL*gr->gdet*cons_L.p + sXL*(vline_L - S_L);
    double ER = -(beta_n + S_R), EL = -(beta_n + S_L);
    double a = (EL*AR - AL*ER), b = -(AL*CR + BL*ER - CL*AR - EL*BR), c = BR*CL - BL*CR, lamstar;
    if(a == 0)
    {
        if(b == 0) {return 1;}
        lamstar = -c/b;
    }
    else
    {
        double discr = b*b - 4*a*c;
        if(discr < 0) {return 1;}
        if(discr == 0) {lamstar = -0.5*b/a;}
        else
        {
            double tmp = (b > 0) ? -0.5*(b + sqrt(discr)) : -0.5*(b - sqrt(discr)); /* the negative root (Mignone & Bodo) */
            lamstar = c/tmp;
            if(fabs(lamstar) > 1) {if(gr_limiter) {lamstar = sqrt(1. - 1./GR_GAMMA_LIMITER/GR_GAMMA_LIMITER);} else {return 1;}}
        }
    }
    *P_M = (AL*lamstar + BL)/((CL + EL*lamstar)*gr->gdet); /* the pressure is continuous across the contact */
    *S_M = gr->alpha*lamstar - beta_n;
    return 0;
}

/* HLL fluxes of the hydro system (lab frame, direction n) */
KOKKOS_INLINE_FUNCTION void gr_hydro_hll_flux(const struct gr_conservatives &cons_L, const struct gr_conservatives &cons_R, double vline_L, double vline_R, double S_L, double S_R,
                                              double v_frame_line, const struct gr_metric *gr, double beta_n, const Vec3<double> &n_unit, struct Riemann_outputs *Riemann_out)
{
    double alpha = gr->alpha, gdet = gr->gdet;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    if(S_L > v_frame_line)
    {
        Riemann_out->Fluxes.rho = cons_L.D*(vline_L - v_frame_line);
        Riemann_out->Fluxes.p   = cons_L.tau*(vline_L - v_frame_line) + gdet*cons_L.p*(vline_L + beta_n);
        Riemann_out->Fluxes.v   = (vline_L - v_frame_line)*cons_L.S + (alpha*gdet*cons_L.p)*n_unit;
        return;
    }
    if(S_R < v_frame_line)
    {
        Riemann_out->Fluxes.rho = cons_R.D*(vline_R - v_frame_line);
        Riemann_out->Fluxes.p   = cons_R.tau*(vline_R - v_frame_line) + gdet*cons_R.p*(vline_R + beta_n);
        Riemann_out->Fluxes.v   = (vline_R - v_frame_line)*cons_R.S + (alpha*gdet*cons_R.p)*n_unit;
        return;
    }
#endif
    struct gr_conservatives U_hll;
    double facLR = 1./(S_R - S_L);
    U_hll.D   = (S_R*cons_R.D - S_L*cons_L.D + cons_L.D*vline_L - cons_R.D*vline_R)*facLR;
    U_hll.tau = (S_R*cons_R.tau - S_L*cons_L.tau + cons_L.tau*vline_L - cons_R.tau*vline_R + gdet*(cons_L.p*(vline_L + beta_n) - cons_R.p*(vline_R + beta_n)))*facLR;
    U_hll.S   = facLR*(S_R*cons_R.S - S_L*cons_L.S + vline_L*cons_L.S - vline_R*cons_R.S + (alpha*gdet*(cons_L.p - cons_R.p))*n_unit);
    double F_D = (S_R*cons_L.D*vline_L - S_L*cons_R.D*vline_R + S_L*S_R*(cons_R.D - cons_L.D))*facLR;
#ifdef HYDRO_MESHLESS_FINITE_MASS
    v_frame_line = F_D/U_hll.D; /* the HLL equivalent of the zero-mass-flux face velocity */
    Riemann_out->Fluxes.rho = 0;
#else
    Riemann_out->Fluxes.rho = F_D - v_frame_line*U_hll.D;
#endif
    Riemann_out->Fluxes.p = facLR*(S_R*(cons_L.tau*(vline_L - S_L) + gdet*cons_L.p*(vline_L + beta_n)) - S_L*(cons_R.tau*(vline_R - S_R) + gdet*cons_R.p*(vline_R + beta_n))) - v_frame_line*U_hll.tau;
    Riemann_out->Fluxes.v = facLR*(S_R*((vline_L - S_L)*cons_L.S + (alpha*gdet*cons_L.p)*n_unit) - S_L*((vline_R - S_R)*cons_R.S + (alpha*gdet*cons_R.p)*n_unit)) - v_frame_line*U_hll.S;
}

#if (GENERAL_RELATIVITY==1)
/* HLLC fluxes of the hydro system (Mignone & Bodo 2005), finite-volume form */
KOKKOS_INLINE_FUNCTION void gr_hydro_hllc_flux(const struct gr_conservatives &cons_L, const struct gr_conservatives &cons_R, double vline_L, double vline_R, double S_L, double S_R,
                                               double v_frame_line, const struct gr_metric *gr, double beta_n, const Vec3<double> &n_unit, struct Riemann_outputs *Riemann_out)
{
    double alpha = gr->alpha, gdet = gr->gdet, S_M = Riemann_out->S_M, P_M = Riemann_out->P_M;
    if((S_M == S_L) || (S_M == S_R))
    {
        Riemann_out->Fluxes.rho = 0;
        Riemann_out->Fluxes.p = gdet*P_M*S_M;
        Riemann_out->Fluxes.v = (alpha*gdet*P_M)*n_unit;
        return;
    }
    if(S_L > v_frame_line)
    {
        Riemann_out->Fluxes.rho = cons_L.D*(vline_L - v_frame_line);
        Riemann_out->Fluxes.p   = cons_L.tau*(vline_L - v_frame_line) + gdet*cons_L.p*(vline_L + beta_n);
        Riemann_out->Fluxes.v   = (vline_L - v_frame_line)*cons_L.S + (alpha*gdet*cons_L.p)*n_unit;
    }
    else if(v_frame_line <= S_M)
    {
        double v_fac = (S_L - vline_L)/(S_L - S_M); if(v_fac < 0) {v_fac = 0;}
        struct gr_conservatives U;
        U.D   = v_fac*cons_L.D;
        U.S   = v_fac*cons_L.S + (alpha*gdet*(cons_L.p - P_M)/(S_M - S_L))*n_unit;
        U.tau = v_fac*cons_L.tau + gdet*((vline_L + beta_n)*cons_L.p - (S_M + beta_n)*P_M)/(S_M - S_L);
        Riemann_out->Fluxes.rho = U.D*(S_M - v_frame_line);
        Riemann_out->Fluxes.p   = U.tau*(S_M - v_frame_line) + gdet*P_M*(S_M + beta_n);
        Riemann_out->Fluxes.v   = (S_M - v_frame_line)*U.S + (alpha*gdet*P_M)*n_unit;
    }
    else if(v_frame_line <= S_R)
    {
        double v_fac = (S_R - vline_R)/(S_R - S_M); if(v_fac < 0) {v_fac = 0;}
        struct gr_conservatives U;
        U.D   = v_fac*cons_R.D;
        U.S   = v_fac*cons_R.S + (alpha*gdet*(P_M - cons_R.p)/(S_R - S_M))*n_unit;
        U.tau = v_fac*cons_R.tau + gdet*((S_M + beta_n)*P_M - (vline_R + beta_n)*cons_R.p)/(S_R - S_M);
        Riemann_out->Fluxes.rho = U.D*(S_M - v_frame_line);
        Riemann_out->Fluxes.p   = U.tau*(S_M - v_frame_line) + gdet*P_M*(S_M + beta_n);
        Riemann_out->Fluxes.v   = (S_M - v_frame_line)*U.S + (alpha*gdet*P_M)*n_unit;
    }
    else
    {
        Riemann_out->Fluxes.rho = cons_R.D*(vline_R - v_frame_line);
        Riemann_out->Fluxes.p   = cons_R.tau*(vline_R - v_frame_line) + gdet*cons_R.p*(vline_R + beta_n);
        Riemann_out->Fluxes.v   = (vline_R - v_frame_line)*cons_R.S + (alpha*gdet*cons_R.p)*n_unit;
    }
}
#endif

/* the hydro solver, in the lab frame along n_unit */
KOKKOS_INLINE_FUNCTION void gr_hydro_riemann(const struct Input_vec_Riemann &Riemann_vec, struct Riemann_outputs *Riemann_out, const Vec3<double> &n_unit, int gr_limiter)
{
    const struct gr_metric *gr = &(Riemann_vec.gr);
    double beta_n = dot(gr->beta, n_unit), gnn = gr->ugamma.quad(n_unit); /* gamma^nn */
    double vline_L = gr->alpha*dot(Riemann_vec.L.v, n_unit) - beta_n, vline_R = gr->alpha*dot(Riemann_vec.R.v, n_unit) - beta_n;
    double v_frame_line = gr->alpha*dot(Riemann_vec.v_frame, n_unit) - beta_n;
    struct gr_face_state L, R;
    gr_face_state_from_primitives(Riemann_vec.L, gr, gr_limiter, &L);
    gr_face_state_from_primitives(Riemann_vec.R, gr, gr_limiter, &R);
    double lam_l[2], lam_r[2];
    gr_lambda_hydro(L.v, Riemann_vec.L.cs, gnn, beta_n, gr->alpha, n_unit, L.W, lam_l, gr_limiter);
    gr_lambda_hydro(R.v, Riemann_vec.R.cs, gnn, beta_n, gr->alpha, n_unit, R.W, lam_r, gr_limiter);
    double S_L = DMIN(lam_r[1], lam_l[1]), S_R = DMAX(lam_r[0], lam_l[0]);
    if(gr_starstate(L.cons, R.cons, n_unit, vline_L, vline_R, gr, beta_n, gnn, S_L, S_R, &(Riemann_out->S_M), &(Riemann_out->P_M), gr_limiter) != 0) {Riemann_out->P_M = -1; return;}
    if(Riemann_out->P_M < 0) {return;}
#if (GENERAL_RELATIVITY==1)
#ifdef HYDRO_MESHLESS_FINITE_MASS
    /* the face moves with the contact (v_frame = S_M): no mass flux, pressure work and force of the star state (the finite-volume
       form below with v_frame_line = S_M; the C/C++ GR code omitted the alpha*sqrt(gamma) factors here) */
    Riemann_out->Fluxes.rho = 0;
    Riemann_out->Fluxes.v = (gr->alpha*gr->gdet*Riemann_out->P_M)*n_unit;
    Riemann_out->Fluxes.p = gr->gdet*Riemann_out->P_M*(Riemann_out->S_M + beta_n);
#else
    gr_hydro_hllc_flux(L.cons, R.cons, vline_L, vline_R, S_L, S_R, v_frame_line, gr, beta_n, n_unit, Riemann_out);
#endif
#else
    gr_hydro_hll_flux(L.cons, R.cons, vline_L, vline_R, S_L, S_R, v_frame_line, gr, beta_n, n_unit, Riemann_out);
#endif
}
#endif /* MAGNETIC */

/* ------------------------------------------------------------------ entry point */
/* the GR replacement of Riemann_solver(): primitives of both sides in Riemann_vec (lab frame), fluxes back in the lab frame */
KOKKOS_INLINE_FUNCTION void gr_riemann_solver(struct Input_vec_Riemann Riemann_vec, struct Riemann_outputs *Riemann_out, const Vec3<double> &n_unit, int gr_limiter)
{
    if((Riemann_vec.L.p < 0 && Riemann_vec.R.p < 0) || (Riemann_vec.L.rho < 0) || (Riemann_vec.R.rho < 0)) {Riemann_out->P_M = 0; return;}
#ifndef EOS_GENERAL
    /* sound speeds and internal energies from pressure and rest-mass density (ideal gas, relativistic: c_s^2 = Gamma P/(rho0 h)) */
    Riemann_vec.L.u = Riemann_vec.L.p/((GAMMA_DEFAULT - 1.)*Riemann_vec.L.rho);
    Riemann_vec.R.u = Riemann_vec.R.p/((GAMMA_DEFAULT - 1.)*Riemann_vec.R.rho);
    Riemann_vec.L.cs = gr_sound_speed(GAMMA_DEFAULT, Riemann_vec.L.p, Riemann_vec.L.rho, Riemann_vec.L.u);
    Riemann_vec.R.cs = gr_sound_speed(GAMMA_DEFAULT, Riemann_vec.R.p, Riemann_vec.R.rho, Riemann_vec.R.u);
#endif
#ifdef MAGNETIC
    /* rotate the states, the face velocity, the shift and the metric to the face frame (x = normal), solve, rotate the fluxes back */
    struct rotation_matrix rot_matrix;
    rotate_states_to_face(&Riemann_vec, n_unit, &rot_matrix);
    Mat3<double> Rm; Rm[0] = rot_matrix.n; Rm[1] = rot_matrix.m; Rm[2] = rot_matrix.p;
    Riemann_vec.v_frame = Rm.matvec(Riemann_vec.v_frame);
    Riemann_vec.gr.beta = Rm.matvec(Riemann_vec.gr.beta);
    Riemann_vec.gr.gamma = gr_rotate_tensor(Riemann_vec.gr.gamma, Rm);
    Riemann_vec.gr.ugamma = gr_rotate_tensor(Riemann_vec.gr.ugamma, Rm);
#if (GENERAL_RELATIVITY==0)
    gr_mhd_riemann_hll(Riemann_vec, Riemann_out, gr_limiter);
#else
    gr_mhd_riemann_gforce(Riemann_vec, Riemann_out, gr_limiter);
#endif
    rotate_fluxes_back_to_lab(Riemann_out, rot_matrix);
#else
    gr_hydro_riemann(Riemann_vec, Riemann_out, n_unit, gr_limiter);
#endif
#ifdef SAVE_FACE_DENSITY
    Riemann_out->Face_Density = 0.5*(Riemann_vec.L.rho + Riemann_vec.R.rho);
#endif
}

#endif /* GENERAL_RELATIVITY */
#endif /* GR_RIEMANN_FUNCTIONS_H */
