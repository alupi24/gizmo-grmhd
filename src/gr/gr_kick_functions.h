/* gr/gr_kick_functions.h -- the GR integration of the conserved variables: metric source terms, kick, drift and predict.
 *
 * The hydro loop returns the rates of the conserved variables (DtEnergy = dtau/dt, DtMomentum = dS_j/dt, DtB, DtVolume);
 * the kick advances (D, S_j, tau, sqrt(gamma)B) with those rates and the metric source terms, inverts the result
 * (con2prim) and converts it into the rates of the primitive velocity and specific energy (HydroAccel, DtInternalEnergy)
 * that the generic kick/drift of the code apply. The drift re-inverts the kicked and predicted conserved variables with
 * the metric at the new position (implicit position correction along the path), and the predict advances the predicted
 * primitives from the conserved rates. Lupi 2023, MNRAS, 519, 1115; Fedrigo & Lupi 2025, A&A, 700, A112.
 *
 * Every function is a KOKKOS_INLINE_FUNCTION on the particle/cell pointers it is given (the drift runs on the device in
 * drift_particle_impl, the kick on the host in do_the_kick): no globals but All. A failed inversion prints and calls
 * endrun (the device variant traps, the host one requests a controlled stop).
 *
 * Port of gr/gr_kickdrift.cc and gr/gr_source.h of the C/C++ GIZMO branches (A. Lupi, G. Fedrigo). The superposed
 * binary metric (GR_METRIC=2) is not ported. */
#ifndef GR_KICK_FUNCTIONS_H
#define GR_KICK_FUNCTIONS_H
#ifdef GENERAL_RELATIVITY

#ifdef GR_DRIFT_HOLD_INTERNAL_ENERGY
/* Inversion with the energy slaved to the momentum: the conserved D, S_j (and sqrt(gamma) B^i) and the specific internal energy u are
   held fixed, and the Lorentz factor, the rest-mass density, the pressure and tau are recomputed with the metric gr. It is used by the
   drift (the metric changes along the path at fixed S_j and u), the kick removing the corresponding rate from the Valencia energy
   source (gr_drift_tau_path_rate). Carrying tau through the drift unchanged leaves (D, S_j, tau) inconsistent with the new metric by
   O(dt), which makes the velocity of the drifted state first-order accurate (circular geodesics in Kerr-Schild coordinates drift);
   for a test particle the modified scheme is the Hamiltonian splitting in the covariant momentum, second order. Ideal gas: h = 1 + Gamma u does not depend on rho0,
   so S_i S^i = (D h)^2 (W^2 - 1) gives W directly; with a field the rest-frame terms depend on v, so W is iterated.
   On input prims->v is the velocity before the drift (initial guess); returns the Lorentz factor, or -1 if the iteration failed. */
KOKKOS_INLINE_FUNCTION double gr_drift_reinvert(const struct gr_conservatives *cons, double u, double gamma_eos, const struct gr_metric *gr, int gr_limiter, struct gr_primitives *prims, struct gr_conservatives *cons_new)
{
    const double D = cons->D, h = 1. + gamma_eos*u, gdet = gr->gdet;
    Vec3<double> v = prims->v;
    double W = lorentz_factor(gr->gamma, v, gr_limiter);
#ifdef MAGNETIC
    Vec3<double> Bprim = cons->B/gdet;
    double ab0, bsq; Vec3<double> b;
    for(int it = 0; it < 50; it++)
    {
        Vec3<double> vtmp = v;
        bsq = gr_bfield_restframe(gr, Bprim, vtmp, &ab0, &b, gr_limiter);
        double rho0 = D/(gdet*W), hstar = h + bsq/rho0;
        Vec3<double> v_cov = (cons->S + (gdet*ab0)*b)/(D*hstar*W);
        Vec3<double> v_new = gr->ugamma.matvec(v_cov);
        double W_new = lorentz_factor(gr->gamma, v_new, gr_limiter);
        double err = (v_new - v).norm_sq();
        v = v_new; W = W_new;
        if(err < 1e-26) {break;}
        if(it == 49) {return -1;}
    }
    {Vec3<double> vtmp = v; bsq = gr_bfield_restframe(gr, Bprim, vtmp, &ab0, &b, gr_limiter);}
    prims->B = Bprim; prims->ab0 = ab0; prims->b = b; prims->bsq = bsq;
#else
    double S2 = gr->ugamma.quad(cons->S);
    W = sqrt(1. + S2/(D*h*D*h));
    v = gr->ugamma.matvec(cons->S)/(D*h*W);
    lorentz_factor(gr->gamma, v, gr_limiter);
#endif
    prims->v = v;
    prims->rho0 = D/(gdet*W);
    prims->u = u;
    prims->p = gr_pressure_rho0_u(prims->rho0, u);
    gr_get_conservative_variables(prims, cons_new, W, gr);
    return W;
}

/* rate of change of tau along the drift of cell i at fixed (D, S_j, sqrt(gamma) B^i, u), (d tau/d x^k) v^k with v the coordinate
   velocity of the state (the particle is drifted with it): the Valencia energy source is the total rate along the path, of which this part is now
   produced by the re-inversion of the drifted state, so it is removed from the source applied in the kick. Centred finite differences
   of the slaved inversion along the direction of motion. */
KOKKOS_INLINE_FUNCTION double gr_drift_tau_path_rate(int i, const struct gr_conservatives *cons, const struct gr_primitives *prims, struct particle_data *pp, struct gas_cell_data *cell)
{
    /* chain rule in the coordinates: the drift displaces the particle by dt vc^k, so the rate is (d tau/d x^k) vc^k with vc^k = dx^k/dt
       the coordinate velocity of the state; the finite-difference step is the coordinate displacement GR_EPS_TOLERANCE vc/|vc| (|vc| the
       Euclidean norm, a parametrisation of the step that cancels out: no spatial metric is involved) */
    const Vec3<double> vc = Vcoord(prims->v, &(cell[i].gr));
    const double vnorm = vc.norm();
    if(vnorm == 0) {return 0;}
    const Vec3<double> dx = (GR_EPS_TOLERANCE/vnorm)*vc;
    const double gamma_eos = cell[i].gamma_eos_value();
    double tau_pm[2];
    for(int k = 0; k < 2; k++)
    {
        struct gr_metric gr_k; struct gr_primitives prims_k = *prims; struct gr_conservatives cons_k;
        set_metric_coefficients((k == 0) ? pp[i].Pos + dx : pp[i].Pos - dx, pp[i].Potential, &gr_k);
        if(gr_drift_reinvert(cons, prims->u, gamma_eos, &gr_k, cell[i].gr.limiter, &prims_k, &cons_k) < 0) {return 0;}
        tau_pm[k] = cons_k.tau;
    }
    return (tau_pm[0] - tau_pm[1])*(0.5*vnorm/GR_EPS_TOLERANCE); /* [f(x+dx) - f(x-dx)] / (2 |dx|) * |vc| */
}
#endif

/* metric source terms of the momentum and energy equations at cell i, from the stress-energy tensor of the primitives
   (v may be limited in place where required): S_j: (alpha sqrt(gamma)/2) T^{mu nu} d_j g_{mu nu};
   tau: alpha sqrt(gamma) [T^{mu nu} K_ij-terms - (T^00 beta^j + T^0j) d_j alpha], with K_ij = -alpha Gamma^0_ij */
KOKKOS_INLINE_FUNCTION void gr_source_term(int i, struct gr_primitives *prims, struct gr_conservatives *Src, struct particle_data *pp, struct gas_cell_data *cell)
{
    const struct gr_metric *gr = &(cell[i].gr);
    const double alpha = gr->alpha, alphainv2 = 1./(alpha*alpha);
    double Pstar = prims->p, rhoh = prims->rho0*(1 + prims->u) + prims->p;
    double Wfac = lorentz_factor(gr->gamma, prims->v, gr->limiter);
    /* four-velocity u^mu = (W/alpha, W (v - beta/alpha)); v is the velocity seen by the Eulerian observer */
    double u0 = Wfac/alpha;
    Vec3<double> u_sp = Wfac*(prims->v - gr->beta/alpha);
#ifdef MAGNETIC
    /* contravariant rest-frame field b^mu = (W v.B/alpha, B/W + (v.B) u^i) [note: this b^0 is the true component; the ab0 stored in
       gr_primitives and the face states is alpha*b^0], total enthalpy and pressure */
    double vdotB = dot(gr->gamma.matvec(prims->v), prims->B);
    rhoh += prims->bsq; Pstar += 0.5*prims->bsq;
    double bcon0 = Wfac*vdotB/alpha;
    Vec3<double> bcon_sp = prims->B/Wfac + vdotB*u_sp;
#endif
    /* T^{mu nu} = rho h u u - b b + P* g^{mu nu}, with g^00 = -1/alpha^2, g^0i = beta^i/alpha^2, g^ij = gamma^ij - beta^i beta^j/alpha^2 */
    struct gr_sym4 T;
    T.tt = rhoh*u0*u0 - Pstar*alphainv2;
    T.ti = (rhoh*u0)*u_sp + (Pstar*alphainv2)*gr->beta;
    T.ij = rhoh*outer_product(u_sp) + Pstar*(gr->ugamma - alphainv2*outer_product(gr->beta));
#ifdef MAGNETIC
    T.tt -= bcon0*bcon0; T.ti -= bcon0*bcon_sp; T.ij -= outer_product(bcon_sp);
#endif
    /* metric derivatives at the cell */
    SymmetricTensor2<double> Kij; struct gr_sym4 dg[3]; Vec3<double> d_alpha; Mat3<double> d_beta;
    get_metric_derivatives(pp[i].Pos, gr, Kij, dg, d_alpha, d_beta);
    Src->D = 0.; /* no mass source */
    Src->S = (0.5*alpha*gr->gdet)*Vec3<double>{gr_sym4_contract(T, dg[0]), gr_sym4_contract(T, dg[1]), gr_sym4_contract(T, dg[2])};
    double t = T.tt*Kij.quad(gr->beta) + 2.*dot(T.ti, Kij.matvec(gr->beta)) + Kij.double_contract(T.ij) - dot(T.tt*gr->beta + T.ti, d_alpha);
#ifndef SELFGRAVITY_OFF
    /* Newtonian self-gravity as a metric perturbation */
    Src->S += (alpha*gr->gdet*T.tt)*pp[i].GravAccel;
    t += dot(T.tt*gr->beta + T.ti, pp[i].GravAccel)/alpha;
#endif
    Src->tau = alpha*gr->gdet*t;
    Src->p = 0;
#ifdef MAGNETIC
    Src->B = Vec3<double>{0,0,0};
#ifdef DIVBCLEANING_GR
    Src->phi = 0;
#endif
#endif
}

/* adds the source terms over dtime to the conserved variables (predictor, plus the implicit corrector when do_implicit) and
   returns the updated primitives in prims. Returns 0, or -1 after a failed inversion (reported, endrun requested) */
KOKKOS_INLINE_FUNCTION int gr_add_source_term(int i, struct gr_primitives *prims, struct gr_conservatives *cons, double dtime, int do_implicit, struct particle_data *pp, struct gas_cell_data *cell)
{
#ifdef BOX_BND_PARTICLES
    if(pp[i].ID == 0) {return 0;} /* no GR gravity for frozen boundary particles */
#endif
    const struct gr_metric *gr = &(cell[i].gr);
    struct gr_conservatives Src;
    gr_source_term(i, prims, &Src, pp, cell);
#ifdef GR_DRIFT_HOLD_INTERNAL_ENERGY
    Src.tau -= gr_drift_tau_path_rate(i, cons, prims, pp, cell); /* the part of the energy source the drift accounts for */
#endif
    /* predictor */
    cons->S += dtime*Src.S;
    cons->tau += dtime*Src.tau;
    double res_con2prim = gr_get_primitive_variables(cons, prims, gr, gr->limiter);
    if(res_con2prim < 0)
    {
        printf("GR: failed inversion (code %g) adding the sources: ID %llu at %g %g %g\n", res_con2prim, (unsigned long long) pp[i].ID, pp[i].Pos[0], pp[i].Pos[1], pp[i].Pos[2]);
        endrun(9993000); return -1;
    }
    if(do_implicit)
    {
        /* corrector: re-evaluate the sources with the updated primitives, iterated to convergence of the velocity */
        struct gr_conservatives Src_corr = Src;
        Vec3<double> p_err;
        double error = MAX_REAL_NUMBER; int iter = 0;
        do
        {
            p_err = cons->S/cons->D;
            Src = Src_corr;
            gr_source_term(i, prims, &Src_corr, pp, cell);
#ifdef GR_DRIFT_HOLD_INTERNAL_ENERGY
            Src_corr.tau -= gr_drift_tau_path_rate(i, cons, prims, pp, cell);
#endif
            cons->S += dtime*(Src_corr.S - Src.S);
            cons->tau += dtime*(Src_corr.tau - Src.tau);
            res_con2prim = gr_get_primitive_variables(cons, prims, gr, gr->limiter);
            if(res_con2prim < 0)
            {
                printf("GR: failed inversion (code %g) in the source corrector: ID %llu at %g %g %g\n", res_con2prim, (unsigned long long) pp[i].ID, pp[i].Pos[0], pp[i].Pos[1], pp[i].Pos[2]);
                endrun(9993001); return -1;
            }
            p_err = cons->S/cons->D - p_err;
            error = fmax(fmax(p_err[0]*p_err[0], p_err[1]*p_err[1]), p_err[2]*p_err[2]);
            error = sqrt(cons->D*cons->D*error/cons->S.norm_sq());
            iter++;
            if(iter > IMPLICIT_MAX_NITER)
            {
                printf("GR WARNING: %d iterations adding the sources: S/D %2.15g %2.15g %2.15g rho0 %g u %g error %g\n", IMPLICIT_MAX_NITER, cons->S[0]/cons->D, cons->S[1]/cons->D, cons->S[2]/cons->D, prims->rho0, prims->u, error);
                break;
            }
        } while(error > IMPLICIT_SOLVER);
    }
    return 0;
}

/* the primitives of cell i from its (kicked or predicted) state: Eulerian velocity of vcoord, mass and volume, specific energy, field B*V */
KOKKOS_INLINE_FUNCTION double gr_primitives_of_cell(int i, const Vec3<double> &vcoord, double mass, double Vi, double u, const Vec3<double> &BV, struct gr_primitives *Prims, struct particle_data *pp, struct gas_cell_data *cell)
{
    const struct gr_metric *gr = &(cell[i].gr);
    Prims->v = Vphys(vcoord, gr);
    double Wlorentz = lorentz_factor(gr->gamma, Prims->v, gr->limiter);
    Prims->rho0 = mass/Wlorentz/gr->gdet/Vi;
    Prims->u = u;
    Prims->p = gr_pressure_rho0_u(Prims->rho0, Prims->u);
#ifdef MAGNETIC
    Prims->B = BV/(Vi*gr->gdet);
    Prims->bsq = gr_bfield_restframe(gr, Prims->B, Prims->v, &(Prims->ab0), &(Prims->b), gr->limiter);
#else
    (void)BV;
#endif
    return Wlorentz;
}

/* the kick of the conserved variables over dt_kick (mode 0: first half-step, 1: second half-step), converted into the rates
   HydroAccel = dv/dt (coordinate velocity) and DtInternalEnergy = du/dt that the generic kick applies */
KOKKOS_INLINE_FUNCTION void gr_kick_impl(int i, double dt_kick, double mass_old, double mass_new, int mode, struct particle_data *pp, struct gas_cell_data *cell)
{
    const struct gr_metric *gr = &(cell[i].gr);
    struct gr_primitives Prims;
    struct gr_conservatives Cons;
    double Vi, Vnew; /* the volume before and after the kick, from the volume variation of the Riemann problems */
    if(mode == 1) {Vnew = mass_new/cell[i].Density; Vi = Vnew*exp(-cell[i].DtVolume*dt_kick/Vnew);}
    else          {Vi = mass_old/cell[i].Density; Vnew = Vi*exp(cell[i].DtVolume*dt_kick/Vi);}
#ifdef GR_TEST_PARTICLE
    mass_old = mass_new = Vi = Vnew = 1;
#endif
    Vec3<double> BV = {0,0,0};
#ifdef MAGNETIC
    BV = cell[i].B;
#endif
    double Wlorentz = gr_primitives_of_cell(i, pp[i].Vel, mass_old, Vi, cell[i].InternalEnergy, BV, &Prims, pp, cell);
    gr_get_conservative_variables(&Prims, &Cons, Wlorentz, gr);
    /* source terms, then the fluxes */
#if !(defined(HYDRO_MESHLESS_FINITE_MASS) && (GENERAL_RELATIVITY==1))
    if(mode == 0)
#endif
    {if(gr_add_source_term(i, &Prims, &Cons, dt_kick, 1-mode, pp, cell) < 0) {return;}}
    if(pp[i].Mass == 0) {return;}
    Cons.D = mass_new/Vnew;
    Cons.tau = (Cons.tau*Vi + cell[i].DtEnergy*dt_kick)/Vnew;
    Cons.S = (Vi*Cons.S + dt_kick*cell[i].DtMomentum)/Vnew;
#ifdef MAGNETIC
    Cons.B = (Vi*Cons.B + dt_kick*cell[i].DtB)/Vnew;
#endif
    double res_con2prim = gr_get_primitive_variables(&Cons, &Prims, gr, gr->limiter);
    if(res_con2prim < 0)
    {
        printf("GR: failed inversion (code %g) in the kick: ID %llu at %g %g %g, kernel radius %g\n", res_con2prim, (unsigned long long) pp[i].ID, pp[i].Pos[0], pp[i].Pos[1], pp[i].Pos[2], pp[i].KernelRadius);
        endrun(9992405); return;
    }
#if !(defined(HYDRO_MESHLESS_FINITE_MASS) && (GENERAL_RELATIVITY==1))
    if(mode == 1) {if(gr_add_source_term(i, &Prims, &Cons, dt_kick, 1-mode, pp, cell) < 0) {return;} if(pp[i].Mass == 0) {return;}}
#endif
    /* the rates from the old and new primitives */
    double rate_fac = 1./dt_kick;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    rate_fac *= cell[i].MassTrue/pp[i].Mass; /* account for the mass change */
#endif
    cell[i].DtInternalEnergy = (Prims.u - cell[i].InternalEnergy)*rate_fac;
#ifdef GR_ENTROPY_SWITCH
    if(cell[i].Entropy != 0)
    {
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        cell[i].Entropy = cell[i].Entropy + (cell[i].DtEntropy*dt_kick)/Vi;
#endif
        if(0.5*Prims.bsq >= 0.1*Prims.p) /* magnetically dominated: the energy from the entropy function */
        {
            double g = cell[i].gamma_eos_value();
            double EnergyDensFromEntropy = cell[i].Entropy*pow(Prims.rho0, g - 1)/(g - 1);
            cell[i].DtInternalEnergy = (EnergyDensFromEntropy - cell[i].InternalEnergy)*rate_fac;
        }
    }
#endif
    cell[i].HydroAccel = (Vcoord(Prims.v, gr) - pp[i].Vel)*rate_fac;
}

/* the drift of the conserved variables: the kicked and predicted states are re-inverted with the metric at the new position,
   with an implicit correction of the position along the path (the particle moves with the velocity of the drifted state) */
KOKKOS_INLINE_FUNCTION void gr_drift_impl(int i, double dt_drift, struct particle_data *pp, struct gas_cell_data *cell)
{
    struct gr_primitives Prims;
    struct gr_conservatives Cons, ConsPred;
    double Vi = pp[i].Mass/cell[i].Density, res_con2prim, Wlorentz;
    Vec3<double> v_old = pp[i].Vel, v_new;
#ifdef GR_TEST_PARTICLE
    pp[i].Mass = Vi = 1;
    if(pp[i].Pos.norm() <= 1.001) {pp[i].Mass = 0;}
#endif
    Vec3<double> BV = {0,0,0}, BVPred = {0,0,0};
#ifdef MAGNETIC
    BV = cell[i].B; BVPred = cell[i].BPred;
#endif
    /* conserved variables of the predicted and of the kicked state, with the metric of the old position */
    Wlorentz = gr_primitives_of_cell(i, cell[i].VelPred, pp[i].Mass, Vi, cell[i].InternalEnergyPred, BVPred, &Prims, pp, cell);
    gr_get_conservative_variables(&Prims, &ConsPred, Wlorentz, &(cell[i].gr));
    Wlorentz = gr_primitives_of_cell(i, pp[i].Vel, pp[i].Mass, Vi, cell[i].InternalEnergy, BV, &Prims, pp, cell);
    gr_get_conservative_variables(&Prims, &Cons, Wlorentz, &(cell[i].gr));
    /* implicit drift loop, accounting for the metric change along the path */
#ifdef IMPLICIT_SOLVER
    double error = MAX_REAL_NUMBER; int iter = 0;
    do {
#endif
        Vec3<double> x_err = pp[i].Pos;
        set_metric_coefficients(pp[i].Pos, pp[i].Potential, &(cell[i].gr)); /* metric at the predicted location */
#ifdef GR_DRIFT_HOLD_INTERNAL_ENERGY
        {struct gr_conservatives Ctmp; res_con2prim = gr_drift_reinvert(&Cons, cell[i].InternalEnergy, cell[i].gamma_eos_value(), &(cell[i].gr), cell[i].gr.limiter, &Prims, &Ctmp);}
#else
        res_con2prim = gr_get_primitive_variables(&Cons, &Prims, &(cell[i].gr), cell[i].gr.limiter);
#endif
        if(res_con2prim < -1)
        {
            printf("GR: failed inversion (code %g) in the drift: ID %llu at %g %g %g\n", res_con2prim, (unsigned long long) pp[i].ID, pp[i].Pos[0], pp[i].Pos[1], pp[i].Pos[2]);
            endrun(9992406); return;
        }
        v_new = Vcoord(Prims.v, &(cell[i].gr));
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION<5)
        break; /* no position correction except for (quasi-)Lagrangian motion */
#endif
        pp[i].Pos += (0.5*dt_drift)*(v_new - v_old); v_old = v_new; /* position with the corrected velocity */
#ifdef IMPLICIT_SOLVER
        x_err = pp[i].Pos - x_err;
        error = fmax(fmax(fabs(x_err[0]), fabs(x_err[1])), fabs(x_err[2]));
        iter++;
        if(iter > IMPLICIT_MAX_NITER)
        {
            printf("GR WARNING (drift): %d iterations: pos %2.15g %2.15g %2.15g, err %g %g %g / %g\n", IMPLICIT_MAX_NITER, pp[i].Pos[0], pp[i].Pos[1], pp[i].Pos[2], x_err[0], x_err[1], x_err[2], error);
            break;
        }
    } while(error > IMPLICIT_SOLVER);
#endif
    /* metric at the new location, kicked primitives (corrector) */
    set_metric_coefficients(pp[i].Pos, pp[i].Potential, &(cell[i].gr));
#ifdef GR_DRIFT_HOLD_INTERNAL_ENERGY
    res_con2prim = gr_drift_reinvert(&Cons, cell[i].InternalEnergy, cell[i].gamma_eos_value(), &(cell[i].gr), cell[i].gr.limiter, &Prims, &Cons);
#else
    res_con2prim = gr_get_primitive_variables(&Cons, &Prims, &(cell[i].gr), cell[i].gr.limiter);
#endif
    if(res_con2prim < 0) {printf("GR: failed inversion (code %g) in the drift corrector: ID %llu\n", res_con2prim, (unsigned long long) pp[i].ID); endrun(9992407); return;}
    pp[i].Vel = Vcoord(Prims.v, &(cell[i].gr)); /* velocity consistent with the metric change */
#ifdef MAGNETIC
    cell[i].B = (Vi*cell[i].gr.gdet)*Prims.B;
#endif
    cell[i].InternalEnergy = Prims.u;
#ifdef GR_ENTROPY_SWITCH
    if((cell[i].Entropy != 0) && (0.5*Prims.bsq >= 0.1*Prims.p))
    {
        double g = cell[i].gamma_eos_value();
        cell[i].InternalEnergy = cell[i].Entropy*pow(Prims.rho0, g - 1)/(g - 1);
    }
#endif
    /* predicted primitives */
#ifdef GR_DRIFT_HOLD_INTERNAL_ENERGY
    {Vec3<double> vguess = Vphys(cell[i].VelPred, &(cell[i].gr)); Prims.v = vguess;
     res_con2prim = gr_drift_reinvert(&ConsPred, cell[i].InternalEnergyPred, cell[i].gamma_eos_value(), &(cell[i].gr), cell[i].gr.limiter, &Prims, &ConsPred);}
#else
    res_con2prim = gr_get_primitive_variables(&ConsPred, &Prims, &(cell[i].gr), cell[i].gr.limiter);
#endif
    if(res_con2prim < 0) {printf("GR: failed inversion (code %g) in the drift (predicted state): ID %llu\n", res_con2prim, (unsigned long long) pp[i].ID); endrun(9992408); return;}
    Vec3<double> vpred_new = Vcoord(Prims.v, &(cell[i].gr));
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION>=5)
    cell[i].ParticleVel += vpred_new - cell[i].VelPred;
#endif
    cell[i].VelPred = vpred_new;
#ifdef MAGNETIC
    cell[i].BPred = (Vi*cell[i].gr.gdet)*Prims.B;
#endif
    cell[i].InternalEnergyPred = Prims.u;
#ifdef GR_ENTROPY_SWITCH
    if((cell[i].Entropy != 0) && (0.5*Prims.bsq >= 0.1*Prims.p))
    {
        double g = cell[i].gamma_eos_value();
        cell[i].InternalEnergyPred = cell[i].EntropyPred*pow(Prims.rho0, g - 1)/(g - 1);
    }
#endif
}

/* the prediction of the primitives over dt_hydrokick from the conserved rates (and the metric sources), replacing the generic
   VelPred/InternalEnergyPred prediction: HydroAccel and DtInternalEnergy are zeroed so that the generic drift adds nothing */
KOKKOS_INLINE_FUNCTION void gr_predict_impl(int i, double dt_hydrokick, struct particle_data *pp, struct gas_cell_data *cell)
{
    const struct gr_metric *gr = &(cell[i].gr);
    struct gr_primitives Prims;
    struct gr_conservatives ConsPred;
    double Vi = pp[i].Mass/cell[i].Density, mass_new = pp[i].Mass, res_con2prim;
#ifdef GR_TEST_PARTICLE
    mass_new = pp[i].Mass = Vi = 1.0;
#endif
#ifdef PARTICLE_EXCISION
    if(pp[i].Mass == 0) {return;}
#endif
    Vec3<double> BVPred = {0,0,0};
#ifdef MAGNETIC
    BVPred = cell[i].BPred;
#endif
    double Wlorentz = gr_primitives_of_cell(i, cell[i].VelPred, pp[i].Mass, Vi, cell[i].InternalEnergyPred, BVPred, &Prims, pp, cell);
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    mass_new = fmax(mass_new + cell[i].DtMass*dt_hydrokick, 0.5*cell[i].MassTrue);
#endif
    gr_get_conservative_variables(&Prims, &ConsPred, Wlorentz, gr);
    double Vnew = Vi*exp(cell[i].DtVolume*dt_hydrokick/Vi);
    if(gr_add_source_term(i, &Prims, &ConsPred, dt_hydrokick, 0, pp, cell) < 0) {return;}
    if(pp[i].Mass == 0) {return;}
    res_con2prim = gr_get_primitive_variables(&ConsPred, &Prims, gr, gr->limiter);
    if(res_con2prim < 0) {printf("GR: failed inversion (code %g) in the predict: ID %llu at %g %g %g\n", res_con2prim, (unsigned long long) pp[i].ID, pp[i].Pos[0], pp[i].Pos[1], pp[i].Pos[2]); endrun(99924091); return;}
    ConsPred.D = mass_new/Vnew;
    ConsPred.tau = (ConsPred.tau*Vi + cell[i].DtEnergy*dt_hydrokick)/Vnew;
    ConsPred.S = (Vi*ConsPred.S + dt_hydrokick*cell[i].DtMomentum)/Vnew;
#ifdef MAGNETIC
    ConsPred.B = (Vi*ConsPred.B + dt_hydrokick*cell[i].DtB)/Vnew;
#endif
    res_con2prim = gr_get_primitive_variables(&ConsPred, &Prims, gr, gr->limiter);
    if(res_con2prim < 0) {printf("GR: failed inversion (code %g) in the predict (fluxes): ID %llu\n", res_con2prim, (unsigned long long) pp[i].ID); endrun(9992409); return;}
    cell[i].VelPred = Vcoord(Prims.v, gr);
    cell[i].InternalEnergyPred = Prims.u;
#ifdef GR_ENTROPY_SWITCH
    if(cell[i].Entropy != 0)
    {
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        cell[i].EntropyPred = cell[i].EntropyPred + (cell[i].DtEntropy*dt_hydrokick)/Vi;
#endif
        if(0.5*Prims.bsq >= 0.1*Prims.p)
        {
            double g = cell[i].gamma_eos_value();
            cell[i].InternalEnergyPred = cell[i].EntropyPred*pow(Prims.rho0, g - 1)/(g - 1);
        }
    }
#endif
    /* zero the rates, to avoid double counting in the generic drift */
    cell[i].HydroAccel = Vec3<double>{0,0,0};
    cell[i].DtInternalEnergy = 0;
}

#endif /* GENERAL_RELATIVITY */
#endif /* GR_KICK_FUNCTIONS_H */
