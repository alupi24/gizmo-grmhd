/* gr/gr_con2prim_functions.h -- conserved-to-primitive inversion (Noble et al. 2006, 2D Newton-Raphson in W and v^2), device-callable.
 *
 * gr_get_primitive_variables() returns the Lorentz factor on success and a negative code on failure:
 *   -1 : rho0 <= 0 or u <= 0 after the inversion
 *   -3 : W out of range,  -4 : v^2 >= 1,  -5 : the initial Lorentz factor is not finite
 *   -(100*r + 1) : the Newton-Raphson iteration returned r (1 = not converged, 2 = non-finite residual)
 * It never prints: a host caller reports a failure with gr_con2prim_report_failure() (gr/gr_utils.cc).
 *
 * Port of gr/utoprim.cc of the C/C++ GIZMO branches (after the HARM utoprim_2d of Noble et al.), without the host diagnostics. */
#ifndef GR_CON2PRIM_FUNCTIONS_H
#define GR_CON2PRIM_FUNCTIONS_H
#ifdef GENERAL_RELATIVITY

#define GR_C2P_MAX_NEWT_ITER   (30)       /* maximum number of Newton-Raphson iterations */
#define GR_C2P_NEWT_TOL        (1.0e-10)  /* convergence tolerance */
#define GR_C2P_MIN_NEWT_TOL    (1.0e-10)
#define GR_C2P_EXTRA_NEWT_ITER (2)        /* extra iterations once converged */
#define GR_C2P_W_TOO_BIG       (1.e20)    /* gamma^2 (rho0 + u + p) above this is unphysical */
#define GR_C2P_MIN_DV          (1.0e-15)
#define GR_C2P_FAIL_VAL        (MAX_REAL_NUMBER)

KOKKOS_INLINE_FUNCTION double gr_c2p_vsq_calc(double W, const struct kernel_NR *kernel) /* v^2 from W = gamma^2 w */
{
    double Wsq = W*W, Bsq_W = kernel->Bsq + W, Xsq = Bsq_W*Bsq_W;
    return (Wsq*kernel->Qtsq + kernel->QdotBsq*(Bsq_W + W))/(Wsq*Xsq);
}

KOKKOS_INLINE_FUNCTION double gr_c2p_x1_of_x0(double x0, const struct kernel_NR *kernel) /* physical v^2 from W */
{
    double vsq = fabs(gr_c2p_vsq_calc(x0, kernel));
    return (vsq > 1.) ? (1.0 - GR_C2P_MIN_DV) : vsq;
}

/* residuals and Newton step for x = (W, v^2) */
KOKKOS_INLINE_FUNCTION void gr_c2p_func_vsq(const double x[], double dx[], double resid[], double jac[][2], double *f, double *df, const struct kernel_NR *kernel)
{
    double W = x[0], vsq = x[1], dPdvsq, dPdW;
    double p_tmp = gr_eos_info(W, vsq, &dPdW, &dPdvsq, kernel);
    double t2 = -kernel->half_Bsq + dPdvsq, t3 = kernel->Bsq + W, t4 = t3*t3, t23 = 1./W, t16 = kernel->QdotBsq*t23*t23;
    double t11 = kernel->Qtsq - vsq*t4 + t16*(kernel->Bsq + W + W);
    double t18 = -kernel->Qdotn - kernel->half_Bsq*(1.0 + vsq) + 0.5*t16 - W + p_tmp;
    double t24 = t16*t23, t25 = -1.0 + dPdW - t24;
    double t35 = t25*t3 + (kernel->Bsq - 2.0*dPdvsq)*(t16 + vsq*W)*t23;
    double t21 = 1/(t3*t35);
    dx[0] = -(t2*t11 + t4*t18)*t21;
    double t40 = -2*(vsq + t24)*t3;
    dx[1] = -(-t25*t11 + t40*t18)*t21;
    jac[0][0] = t40; jac[0][1] = -t4; jac[1][0] = t25; jac[1][1] = t2;
    resid[0] = t11; resid[1] = t18;
    *df = -resid[0]*resid[0] - resid[1]*resid[1];
    *f = -0.5*(*df);
}

/* 2D Newton-Raphson on (W, v^2): 0 converged, 1 not converged, 2 non-finite residual */
KOKKOS_INLINE_FUNCTION int gr_c2p_twod_newton_raphson(double x[], const struct kernel_NR *kernel)
{
    double f = 1., df = 1., dx[2], x_old[2], resid[2], jac[2][2], errx = 1.;
    int n_iter = 0, i_extra = 0, doing_extra = 0, keep_iterating = 1;
    const double dv = (1. - GR_C2P_MIN_DV);
    while(keep_iterating)
    {
        gr_c2p_func_vsq(x, dx, resid, jac, &f, &df, kernel);
        x_old[0] = x[0]; x_old[1] = x[1];
        x[0] += dx[0]; x[1] += dx[1];
        errx = (x[0] == 0.) ? fabs(dx[0]) : fabs(dx[0]/x[0]);
        if(x[0] < 0.) {x[0] = fabs(x[0]);} else {if(x[0] > GR_C2P_W_TOO_BIG) {x[0] = x_old[0];}}
        if(x[1] < 0.) {x[1] = 0.;} else {if(x[1] > 1.) {x[1] = dv;}}
        if((fabs(errx) <= GR_C2P_NEWT_TOL) && (doing_extra == 0) && (GR_C2P_EXTRA_NEWT_ITER > 0)) {doing_extra = 1;}
        if(doing_extra == 1) {i_extra++;}
        if(((fabs(errx) <= GR_C2P_NEWT_TOL) && (doing_extra == 0)) || (i_extra > GR_C2P_EXTRA_NEWT_ITER) || (n_iter >= (GR_C2P_MAX_NEWT_ITER - 1))) {keep_iterating = 0;}
        n_iter++;
    }
    if((!isfinite(f)) || (!isfinite(df))) {return 2;}
    if(fabs(errx) <= GR_C2P_NEWT_TOL) {return 0;}
    if((fabs(errx) <= GR_C2P_MIN_NEWT_TOL) && (fabs(errx) > GR_C2P_NEWT_TOL)) {return 0;}
    return 1;
}

/* primitives from the (densitized) conserved variables; prims->u and prims->v on input are the initial guess (the previous primitives) */
KOKKOS_INLINE_FUNCTION double gr_get_primitive_variables(const struct gr_conservatives *cons_pred, struct gr_primitives *prims, const struct gr_metric *gr, int gr_limiter)
{
    double x_2d[2], Svec[3], uSvec[3], tau, dens, gamma, inv_sqrtg = 1./gr->gdet, rho0, u, p, rhoh, W_last, W, vsq, QdotB, QdB_o_W;
    int i_increase, retval;
    struct kernel_NR kernel;
    dens = cons_pred->D*inv_sqrtg;
    Svec[0] = cons_pred->S[0]*inv_sqrtg; Svec[1] = cons_pred->S[1]*inv_sqrtg; Svec[2] = cons_pred->S[2]*inv_sqrtg;
    tau = cons_pred->tau*inv_sqrtg;
    raise_lower_indices(gr->ugamma, Svec, uSvec);
#ifdef GR_TEST_PARTICLE
    W = (tau + dens)/sqrt((tau + dens)*(tau + dens) - (uSvec[0]*Svec[0] + uSvec[1]*Svec[1] + uSvec[2]*Svec[2]));
    prims->rho0 = dens/W;
    prims->v[0] = uSvec[0]/(tau + dens); prims->v[1] = uSvec[1]/(tau + dens); prims->v[2] = uSvec[2]/(tau + dens);
    prims->u = 1e-20; prims->p = 0;
    return W;
#endif
#ifdef MAGNETIC
    double Bcov[3];
    prims->B[0] = cons_pred->B[0]*inv_sqrtg; prims->B[1] = cons_pred->B[1]*inv_sqrtg; prims->B[2] = cons_pred->B[2]*inv_sqrtg;
    raise_lower_indices(gr->gamma, prims->B, Bcov);
    kernel.Bsq = Bcov[0]*prims->B[0] + Bcov[1]*prims->B[1] + Bcov[2]*prims->B[2];
    QdotB = Svec[0]*prims->B[0] + Svec[1]*prims->B[1] + Svec[2]*prims->B[2];
    kernel.QdotBsq = QdotB*QdotB;
    kernel.half_Bsq = 0.5*kernel.Bsq;
#else
    kernel.Bsq = 0.0; QdotB = 0.0; kernel.QdotBsq = 0.0; kernel.half_Bsq = 0.0;
#endif
    kernel.Qdotn = -(tau + dens);
    kernel.Qtsq = (uSvec[0]*Svec[0] + uSvec[1]*Svec[1] + uSvec[2]*Svec[2]);
    kernel.D = dens;
    /* initial guess from the previous primitives */
    vsq = vector_sqnorm(gr->gamma, prims->v);
    if(vsq > 1) {gamma = lorentz_factor(gr->gamma, prims->v, gr_limiter); vsq = 1. - 1./gamma/gamma;}
    gamma = 1./sqrt(1.0 - vsq);
    if(!isfinite(gamma)) {return -5.;}
    rho0 = kernel.D/gamma;
    u = prims->u*rho0;
    p = gr_pressure_rho0_u(rho0, prims->u);
    rhoh = rho0 + u + p;
    W_last = rhoh*gamma*gamma;
    i_increase = 0;
    while(((W_last*W_last*W_last*(W_last + 2.*kernel.Bsq) - kernel.QdotBsq*(2.*W_last + kernel.Bsq)) <= W_last*W_last*(kernel.Qtsq - kernel.Bsq*kernel.Bsq)) && (i_increase < 10))
    {
        W_last *= 10.; i_increase++;
    }
    x_2d[0] = fabs(W_last);
    x_2d[1] = gr_c2p_x1_of_x0(W_last, &kernel);
    retval = gr_c2p_twod_newton_raphson(x_2d, &kernel);
    W = x_2d[0]; vsq = x_2d[1];
    if((retval != 0) || (W == GR_C2P_FAIL_VAL)) {return -(retval*100. + 1.);}
    if(W <= 0. || W > GR_C2P_W_TOO_BIG) {return -3.;}
    if(vsq >= 1.) {return -4.;}
    if(gr_limiter == 1) {double RMAX = 1.0 - 1./((double)GR_GAMMA_LIMITER*(double)GR_GAMMA_LIMITER); vsq = 0.5*(RMAX + vsq - fabs(vsq - RMAX));}
    gamma = 1./sqrt(1. - vsq);
    rho0 = kernel.D/gamma;
    rhoh = W*(1. - vsq);
    rhoh = fmax(rho0, rhoh);
    p = gr_pressure_rho0_w(rho0, rhoh);
    u = rhoh - (rho0 + p);
#ifndef EOS_GENERAL
    if(p <= 0 && fabs(p)/rhoh < 1e-10) {u = prims->u*rho0; p = gr_pressure_rho0_u(rho0, u/rho0);}
    if(u <= 0) {u = p/(GAMMA_DEFAULT - 1);}
#endif
    prims->rho0 = rho0;
    prims->u = u/rho0;
    prims->p = p;
    double W_inv = 1./(W + kernel.Bsq);
    prims->v[0] = W_inv*uSvec[0]; prims->v[1] = W_inv*uSvec[1]; prims->v[2] = W_inv*uSvec[2];
#ifdef MAGNETIC
    QdB_o_W = QdotB/W;
    prims->bsq = kernel.Bsq*(1. - vsq) + QdB_o_W*QdB_o_W;
    prims->v[0] += W_inv*QdB_o_W*prims->B[0]; prims->v[1] += W_inv*QdB_o_W*prims->B[1]; prims->v[2] += W_inv*QdB_o_W*prims->B[2];
#endif
    lorentz_factor(gr->gamma, prims->v, gr_limiter); /* limits the velocity in place where required */
    if((rho0 <= 0.) || (u <= 0.)) {return -1;}
    return gamma;
}

#endif /* GENERAL_RELATIVITY */
#endif /* GR_CON2PRIM_FUNCTIONS_H */
