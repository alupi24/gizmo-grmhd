/* gr/gr_utils.cc -- host-side utilities of the GR module: excision around the black hole(s), failure reports.
 *
 * The device-callable physics is in gr/gr_functions.h (included through declarations/cell_data.h). */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"

#ifdef GENERAL_RELATIVITY

/* prints the state of a failed conserved-to-primitive inversion (the con2prim itself never prints) */
void gr_con2prim_report_failure(const char *where, int code, const struct gr_conservatives *cons, const struct gr_primitives *prims, const struct gr_metric *gr)
{
    printf("GR con2prim failed in %s with code %d (task %d)\n", where, code, ThisTask);
    printf(" *dens        = %26.16e \n *sx          = %26.16e \n *sy          = %26.16e \n *sz          = %26.16e \n *tau         = %26.16e \n",
           cons->D, cons->S[0], cons->S[1], cons->S[2], cons->tau);
#ifdef MAGNETIC
    printf(" *Bconsx      = %26.16e \n *Bconsy      = %26.16e \n *Bconsz      = %26.16e \n", cons->B[0], cons->B[1], cons->B[2]);
#endif
    printf(" *rho (guess) = %26.16e \n *velx        = %26.16e \n *vely        = %26.16e \n *velz        = %26.16e \n *epsilon     = %26.16e \n *pressure    = %26.16e \n",
           prims->rho0, prims->v[0], prims->v[1], prims->v[2], prims->u, prims->p);
#ifdef MAGNETIC
    printf(" *Bx          = %26.16e \n *By          = %26.16e \n *Bz          = %26.16e \n *bsq         = %26.16e \n", prims->B[0], prims->B[1], prims->B[2], prims->bsq);
#endif
    printf(" metric: alpha %g gdet %g beta %g %g %g gamma %g %g %g %g %g %g\n", gr->alpha, gr->gdet, gr->beta[0], gr->beta[1], gr->beta[2],
           gr->gamma[0], gr->gamma[1], gr->gamma[2], gr->gamma[3], gr->gamma[4], gr->gamma[5]);
    fflush(stdout);
}

#ifdef PARTICLE_EXCISION
/* damps the velocity (and field) of a cell approaching the horizon, as in the original GR code */
void smooth_conservatives_inside_horizon(int i, double r, double r_EH, double M)
{
    double v_local[3], v_localPred[3]; int k;
    for(k=0;k<3;k++) {v_local[k] = P[i].Vel[k]; v_localPred[k] = CellP[i].VelPred[k];}
    Vphys(v_local, CellP[i].gr.beta, CellP[i].gr.alpha);
    Vphys(v_localPred, CellP[i].gr.beta, CellP[i].gr.alpha);
    double exponent = (0.5*sqrt(GR_EXCISION_RADIUS*M)*(r - r_EH)*(r - r_EH))/(r_EH - 0.4), Chi = exp(-exponent);
    for(k=0;k<3;k++) {v_local[k] *= Chi; v_localPred[k] *= Chi;}
    Vcoord(v_local, CellP[i].gr.beta, CellP[i].gr.alpha);
    Vcoord(v_localPred, CellP[i].gr.beta, CellP[i].gr.alpha);
    for(k=0;k<3;k++) {P[i].Vel[k] = v_local[k]; CellP[i].VelPred[k] = v_localPred[k];}
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION>=5)
    double par_v_local[3]; for(k=0;k<3;k++) {par_v_local[k] = CellP[i].ParticleVel[k];}
    Vphys(par_v_local, CellP[i].gr.beta, CellP[i].gr.alpha);
    for(k=0;k<3;k++) {par_v_local[k] *= Chi;}
    Vcoord(par_v_local, CellP[i].gr.beta, CellP[i].gr.alpha);
    for(k=0;k<3;k++) {CellP[i].ParticleVel[k] = par_v_local[k];}
#endif
#ifdef MAGNETIC
    CellP[i].B *= Chi; CellP[i].BPred *= Chi;
#endif
}

/* sets the mass to 0 inside GR_EXCISION_RADIUS horizon radii of the black hole (at the origin for the Kerr metric) */
void check_boundary(int target)
{
#if (GR_METRIC==1)
    double r_min = sqrt(P[target].Pos[0]*P[target].Pos[0] + P[target].Pos[1]*P[target].Pos[1] + P[target].Pos[2]*P[target].Pos[2]);
    double a[3] = {0, 0, All.MetricParams[1]};
    double horizon = EHradius(All.MetricParams[0], a);
    if(r_min < GR_EXCISION_RADIUS*horizon) {P[target].Mass = 0;}
#endif
}

#endif /* PARTICLE_EXCISION */

#endif /* GENERAL_RELATIVITY */
