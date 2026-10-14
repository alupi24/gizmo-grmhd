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

/* metric source terms of the momentum and energy equations at cell i, from the stress-energy tensor T^{mu nu} of the primitives
   (v may be limited in place where required) */
KOKKOS_INLINE_FUNCTION void gr_source_term(int i, struct gr_primitives *prims, struct gr_conservatives *Src, struct particle_data *pp, struct gas_cell_data *cell)
{
    const struct gr_metric *gr = &(cell[i].gr);
    const double alpha = gr->alpha, alphainv2 = 1./(alpha*alpha);
    double Tmunu[4][4], ucon[4];
    double Pstar = prims->p, rhoh = prims->rho0*(1 + prims->u) + prims->p;
    double Wfac = lorentz_factor(gr->gamma, prims->v, gr->limiter);
    ucon[0] = Wfac/alpha; /* four-velocity; v is the velocity seen by the Eulerian observer */
    for(int j = 0; j < 3; j++) {ucon[j+1] = Wfac*(prims->v[j] - gr->beta[j]/alpha);}
#ifdef MAGNETIC
    double b_con[4], vdotB = dot(gr->gamma.matvec(prims->v), prims->B);
    rhoh += prims->bsq; Pstar += 0.5*prims->bsq;
    b_con[0] = Wfac*vdotB/alpha;
    for(int j = 0; j < 3; j++) {b_con[j+1] = prims->B[j]/Wfac + vdotB*ucon[j+1];}
#endif
    for(int j = 0; j < 4; j++) {for(int k = 0; k < 4; k++)
    {
        Tmunu[j][k] = rhoh*ucon[j]*ucon[k];
#ifdef MAGNETIC
        Tmunu[j][k] -= b_con[j]*b_con[k];
#endif
    }}
    /* the pressure term P g^{mu nu}: g^00 = -1/alpha^2, g^0i = beta^i/alpha^2, g^ij = gamma^ij - beta^i beta^j/alpha^2 */
    Tmunu[0][0] += -Pstar*alphainv2;
    for(int j = 0; j < 3; j++) {Tmunu[0][j+1] = Tmunu[j+1][0] = Tmunu[0][j+1] + Pstar*gr->beta[j]*alphainv2;}
    for(int j = 0; j < 3; j++) {for(int k = j; k < 3; k++) {Tmunu[j+1][k+1] = Tmunu[k+1][j+1] = Tmunu[j+1][k+1] + Pstar*(gr->ugamma[j][k] - gr->beta[j]*gr->beta[k]*alphainv2);}}
    /* metric derivatives at the cell */
    SymmetricTensor2<double> Kij; struct gr_metric4 dg[3]; Vec3<double> d_alpha; Mat3<double> d_beta;
    get_metric_derivatives(pp[i].Pos, gr, Kij, dg, d_alpha, d_beta);
    Src->D = 0.; /* no mass source */
    for(int dir = 0; dir < 3; dir++) /* momentum: (alpha sqrt(gamma)/2) T^{mu nu} d_dir g_{mu nu} */
    {
        double s = 0;
        for(int j = 0; j < 4; j++) {for(int k = 0; k < 4; k++) {s += Tmunu[j][k]*gr_metric4_get(&dg[dir], j, k);}}
#ifndef SELFGRAVITY_OFF
        s += 2*Tmunu[0][0]*pp[i].GravAccel[dir]; /* Newtonian self-gravity as a metric perturbation */
#endif
        Src->S[dir] = 0.5*alpha*gr->gdet*s;
    }
    /* energy: alpha sqrt(gamma) [ K_jk (T^00 beta^j beta^k + 2 T^0j beta^k + T^jk) - (T^00 beta^j + T^0j) d_j alpha ], with K_ij = -alpha Gamma^0_ij */
    double t = 0;
    for(int j = 0; j < 3; j++) {for(int k = 0; k < 3; k++) {t += Kij[j][k]*(Tmunu[0][0]*gr->beta[j]*gr->beta[k] + 2.0*Tmunu[0][j+1]*gr->beta[k] + Tmunu[j+1][k+1]);}}
    for(int j = 0; j < 3; j++) {t += (Tmunu[0][0]*gr->beta[j] + Tmunu[0][j+1])*(-d_alpha[j]);}
#ifndef SELFGRAVITY_OFF
    for(int j = 0; j < 3; j++) {t += (Tmunu[0][0]*gr->beta[j] + Tmunu[0][j+1])*pp[i].GravAccel[j]/alpha;}
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
            cons->S += dtime*(Src_corr.S - Src.S);
            cons->tau += dtime*(Src_corr.tau - Src.tau);
            res_con2prim = gr_get_primitive_variables(cons, prims, gr, gr->limiter);
            if(res_con2prim < 0)
            {
                printf("GR: failed inversion (code %g) in the source corrector: ID %llu at %g %g %g\n", res_con2prim, (unsigned long long) pp[i].ID, pp[i].Pos[0], pp[i].Pos[1], pp[i].Pos[2]);
                endrun(9993001); return -1;
            }
            p_err = cons->S/cons->D - p_err;
            error = 0; for(int k = 0; k < 3; k++) {if(p_err[k]*p_err[k] > error) {error = p_err[k]*p_err[k];}}
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
    Prims->bsq = gr_bfield_restframe(gr, Prims->B, Prims->v, &(Prims->b0), &(Prims->b), gr->limiter);
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
        res_con2prim = gr_get_primitive_variables(&Cons, &Prims, &(cell[i].gr), cell[i].gr.limiter);
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
        error = 0; for(int k = 0; k < 3; k++) {if(fabs(x_err[k]) > error) {error = fabs(x_err[k]);}}
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
    res_con2prim = gr_get_primitive_variables(&Cons, &Prims, &(cell[i].gr), cell[i].gr.limiter);
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
    res_con2prim = gr_get_primitive_variables(&ConsPred, &Prims, &(cell[i].gr), cell[i].gr.limiter);
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
