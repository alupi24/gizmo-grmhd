/***********************************************************************************
    Copyright 2006 Scott C. Noble, Charles F. Gammie, Jonathan C. McKinney, 
                   and Luca Del Zanna.

                        PVS_GRMHD

    This file was derived from PVS_GRMHD.  The authors of PVS_GRMHD include 
    Scott C. Noble, Charles F. Gammie, Jonathan C. McKinney, and Luca Del Zanna.
    PVS_GRMHD is available under the GPL from:
    http://rainman.astro.uiuc.edu/codelib/  

    You are morally obligated to cite the following  paper in his/her 
    scientific literature that results from use of this file:

    [1] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006, 
        Astrophysical Journal, 641, 626.

    PVS_GRMHD is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    PVS_GRMHD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PVS_GRMHD; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

    If the user has any questions, please direct them to Scott C. Noble at 
    scn@astro.rit.edu  . 

***********************************************************************************/


#include "../declarations/allvars.h"
#include "../core/proto.h"

#ifdef GENERAL_RELATIVITY
/* Set this to be 1 if you want debug output */
#define DEBUG_CON2PRIMM (0)


/* Adiabatic index used for the state equation */

#define MAX_NEWT_ITER   (30)       /* Max. # of Newton-Raphson iterations for find_root_2D(); */
#define NEWT_TOL        (1.0e-10)  /* Min. of tolerance allowed for Newton-Raphson iterations */
#define MIN_NEWT_TOL    (1.0e-10)  /* Max. of tolerance allowed for Newton-Raphson iterations */
#define EXTRA_NEWT_ITER (2)

#define W_TOO_BIG       (1.e20)    /* \gamma^2 (\rho_0 + u + p) is assumed
                                      to always be smaller than this.  This
                                      is used to detect solver failures */
#define MIN_DV          (1.0e-15)

#define FAIL_VAL        (MAX_REAL_NUMBER) /* Generic value to which we set variables when a problem arises */

/**************************************************
  The following functions assume a Gamma-law EOS:
***************************************************/

// Declarations: 
static double vsq_calc(double W, struct kernel_NR *kernel);

static int twod_newton_raphson( double x[], struct kernel_NR *kernel);

static  void func_vsq( double [], double [], double [], double [][2], 
                       double *f, double *df, struct kernel_NR *kernel);

static double x1_of_x0(double x0, struct kernel_NR *kernel ) ;


/**********************************************************************/
/**********************************************************************************

  Con2PrimM_pt():
  -----------------------------

   -- Attempts an inversion from GRMHD conserved variables to primitive variables assuming a guess.

   -- Uses the 2D method of Noble et al. (2006): 
       -- Solves for two independent variables (W,v^2) via a 2D
          Newton-Raphson method 
       -- Can be used (in principle) with a general equation of state. 

   -- Minimizes two residual functions using a homemade Newton-Raphson routine. 
       -- It is homemade so that it can catch exceptions and handle them correctly, plus it is 
          optimized for this particular  problem. 

  -- Note that the notation used herein is that of Noble et al. (2006) except for the argument 
     list. 


INPUT:  (using GRHydro variable defintions)

   s[x,y,z]   =  scons[0,1,2]  = \alpha \sqrt(\gamma) T^0_i 
   dens, tau  =  as defined in GRHydro and are assumed to be densitized (i.e. with sqrt(\gamma))   
   dens       =  D = \sqrt(\gamma) W \rho    
   tau        =  \alpha^2 \sqrt(\gamma) T^{00} - D 
   g[x,y,z][x,y,x] = spatial metric corresponding to \gamma 
   u[x,y,z][x,y,z] = inverse of the spatial metric, g[x,y,z][x,y,x]
   det        =  sqrt(\gamma)
   B[x,y,z]   =  Bvec[0,1,2] 
   bsq        = b^\mu b_\mu  

   epsnegative = (integer)
               = 0  if rho and epsilon are positive
              != 0  otherwise 


  --  (users should set B[x,y,z] = 0  for hydrodynamic runs) 


OUTPUT:  (using GRHydro variable defintions)
   rho, eps   =  as defined in GRHydro, primitive variables
   vel[x,y,z] =  as defined in GRHydro, primitive variables


RETURN VALUE: of retval = (i*100 + j)  where 
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used) 
                   or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the 
                   given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence 
                   (occurrence of "nan" or "+/-inf" ;
             
         j = 0 -> success 
             1 -> failure: some sort of failure in Newton-Raphson; 
             2 -> failure: unphysical vsq = v^2  value at initial guess;
             3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: v^2 > 1 
             ( used to be  5 -> failure: rho,uu <= 0   but now sets epsnegative to non-zero )

**********************************************************************************/
double gr_get_primitive_variables(struct gr_conservatives *cons_pred,struct gr_primitives *prims,struct gr_metric *gr, int gr_limiter)
{
  double x_2d[2];
  //Un-densitised conserved variables
  double Svec[3],uSvec[3];
  double tau, dens;
  //local Lorentz factor and metric determinant
  double gamma;
  double inv_sqrtg = 1./gr->gdet;
  //Primitive variables
  double rho0,u,p,rhoh,Bcov[3],Bloc[3];
  //Temporary variables
  double gtmp,W_last,W,vsq;

  double QdotB;
  double g_o_WBsq, QdB_o_W;

  //Additional variables
  int  i_increase;
  double retval=0.;
 
  //Kernel of conserved quantities
  struct kernel_NR kernel; 

  /* Assume ok initially: */

#if(DEBUG_CON2PRIMM)
  fprintf(stdout," *dens        = %26.16e \n", cons_pred->D    );
  fprintf(stdout," *sx          = %26.16e \n", cons_pred->S[0] );    
  fprintf(stdout," *sy          = %26.16e \n", cons_pred->S[1] );    
  fprintf(stdout," *sz          = %26.16e \n", cons_pred->S[2] );    
  fprintf(stdout," *tau         = %26.16e \n", cons_pred->tau  );    
#ifdef MAGNETIC
  fprintf(stdout," *Pstar       = %26.16e \n", cons_pred->p    );    
  fprintf(stdout," *Bconsx      = %26.16e \n", cons_pred->B[0] );    
  fprintf(stdout," *Bconsy      = %26.16e \n", cons_pred->B[1] );    
  fprintf(stdout," *Bconsz      = %26.16e \n", cons_pred->B[2] );    
#endif
  fprintf(stdout," *rho         = %26.16e \n", prims->rho0    );    
  fprintf(stdout," *velx        = %26.16e \n", prims->v[0]    );    
  fprintf(stdout," *vely        = %26.16e \n", prims->v[1]    );    
  fprintf(stdout," *velz        = %26.16e \n", prims->v[2]    );    
  fprintf(stdout," *epsilon     = %26.16e \n", prims->u       );
  fprintf(stdout," *pressure    = %26.16e \n", prims->p       );
#ifdef MAGNETIC
  fprintf(stdout," *Bx          = %26.16e \n", prims->B[0]    );    
  fprintf(stdout," *By          = %26.16e \n", prims->B[1]    );    
  fprintf(stdout," *Bz          = %26.16e \n", prims->B[2]    );    
  fprintf(stdout," *bsq         = %26.16e \n", prims->bsq     );    
#endif
  fprintf(stdout," *gxx         = %26.16e \n", gr->gamma[0]   );    
  fprintf(stdout," *gxy         = %26.16e \n", gr->gamma[1]   );    
  fprintf(stdout," *gxz         = %26.16e \n", gr->gamma[2]   );    
  fprintf(stdout," *gyy         = %26.16e \n", gr->gamma[3]   );    
  fprintf(stdout," *gyz         = %26.16e \n", gr->gamma[4]   );    
  fprintf(stdout," *gzz         = %26.16e \n", gr->gamma[5]   );    
  fprintf(stdout," *uxx         = %26.16e \n", gr->ugamma[0]  );    
  fprintf(stdout," *uxy         = %26.16e \n", gr->ugamma[1]  );    
  fprintf(stdout," *uxz         = %26.16e \n", gr->ugamma[2]  );    
  fprintf(stdout," *uyy         = %26.16e \n", gr->ugamma[3]  );    
  fprintf(stdout," *uyz         = %26.16e \n", gr->ugamma[4]  );    
  fprintf(stdout," *uzz         = %26.16e \n", gr->ugamma[5]  );    
  fprintf(stdout," *sdet        = %26.16e \n", gr->gdet       );    
  fflush(stdout);
#endif


  /* First undensitize all conserved variables : */
  dens    = (cons_pred->D   )   * inv_sqrtg;
  Svec[0] = (cons_pred->S[0])   * inv_sqrtg;
  Svec[1] = (cons_pred->S[1])   * inv_sqrtg;
  Svec[2] = (cons_pred->S[2])   * inv_sqrtg;
  tau     = (cons_pred->tau )   * inv_sqrtg;

  /*Raise index of S_i*/
  raise_lower_indices(gr->ugamma,Svec,uSvec);
  //uSvec[0] = gr->ugamma[0]*Svec[0]+gr->ugamma[1]*Svec[1]+gr->ugamma[2]*Svec[2];
  //uSvec[1] = gr->ugamma[1]*Svec[0]+gr->ugamma[3]*Svec[1]+gr->ugamma[4]*Svec[2];
  //uSvec[2] = gr->ugamma[2]*Svec[0]+gr->ugamma[4]*Svec[1]+gr->ugamma[5]*Svec[2];

#ifdef GR_TEST_PARTICLE
  //W=tau/dens+1;
  //W=sqrt(1+(uSvec[0]*Svec[0] + uSvec[1]*Svec[1] + uSvec[2]*Svec[2])/dens/dens);
  //printf("%2.12g %2.12g\n",(uSvec[0]*Svec[0] + uSvec[1]*Svec[1] + uSvec[2]*Svec[2]),(tau+dens)*(tau+dens));
  W=(tau+dens)/sqrt((tau+dens)*(tau+dens)-(uSvec[0]*Svec[0] + uSvec[1]*Svec[1] + uSvec[2]*Svec[2]));
  //printf("My W is %2.12g %2.12g\n",W,tau/dens+1);
  prims->rho0=dens/W;
  /*prims->v[0]=uSvec[0]/dens/W;
  prims->v[1]=uSvec[1]/dens/W;
  prims->v[2]=uSvec[2]/dens/W;
  printf("My v is %2.12g %2.12g %2.12g %2.12g %2.12g %2.12g\n",prims->v[0],prims->v[1],prims->v[2],uSvec[0]/(ta$
  getchar();*/
  prims->v[0]=uSvec[0]/(tau+dens);
  prims->v[1]=uSvec[1]/(tau+dens);
  prims->v[2]=uSvec[2]/(tau+dens);
  prims->u = 1e-20;
  prims->p=0;
  return W;
#endif


#ifdef MAGNETIC
  prims->B[0] = cons_pred->B[0] * inv_sqrtg;
  prims->B[1] = cons_pred->B[1] * inv_sqrtg;
  prims->B[2] = cons_pred->B[2] * inv_sqrtg;
#endif

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
#ifdef MAGNETIC
  raise_lower_indices(gr->gamma,prims->B,Bcov);
  kernel.Bsq = Bcov[0]*prims->B[0]+Bcov[1]*prims->B[1]+Bcov[2]*prims->B[2];
  QdotB = Svec[0]*prims->B[0]+Svec[1]*prims->B[1]+Svec[2]*prims->B[2];
  kernel.QdotBsq = QdotB*QdotB;
  kernel.half_Bsq = 0.5*kernel.Bsq;
#else
  kernel.Bsq = 0.0;
  QdotB = 0.0;
  kernel.QdotBsq = 0.0;
  kernel.half_Bsq = 0.0;
#endif

  kernel.Qdotn = -(tau + dens) ;
  kernel.Qtsq  = (uSvec[0]*Svec[0] + uSvec[1]*Svec[1] + uSvec[2]*Svec[2]);
  kernel.D     = dens;

  /* calculate W from last timestep and use for guess */
  vsq = vector_sqnorm(gr->gamma,prims->v);
  if(vsq>1)
  {
	gamma = lorentz_factor(gr->gamma,prims->v, gr_limiter);
	vsq = 1.-1./gamma/gamma;
  }

  gamma = 1./sqrtl(1.0-vsq);
  if(isnan(gamma)){
   printf("I got NaN because vsq is %2.12g\n",vsq);
   endrun(9992404);
  }
#if(DEBUG_CON2PRIMM)
  printf("My vsq is %2.24g- gamma %2.12g\n",vsq,gamma);
#endif

  // Always calculate rho from D and gamma so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq) . 
  rho0 = kernel.D / gamma ;
  u = prims->u * rho0;
  p = gr_pressure_rho0_u(rho0,prims->u) ;  // EOS
  rhoh = rho0 + u + p ;

#if(DEBUG_CON2PRIMM)
  printf("Initial guess: local %2.12g %2.12g - input %2.12g %2.12g - w %2.12g\n",rho0,p,prims->rho0,prims->p,rhoh);
#endif

  W_last = rhoh*gamma*gamma;

  // Make sure that W is large enough so that v^2 < 1 : 
  i_increase = 0;
  while( ((W_last*W_last*W_last * (W_last + 2.*kernel.Bsq) - 
	kernel.QdotBsq * (2.*W_last + kernel.Bsq)) <= W_last*W_last*(kernel.Qtsq-kernel.Bsq*kernel.Bsq))
         && (i_increase < 10) ) {
    W_last *= 10.;
    i_increase++;
  }
  
  // Calculate W and vsq: 
  x_2d[0] =  fabs( W_last );
  x_2d[1] = x1_of_x0( W_last, &kernel ) ;
  retval = 1.0*twod_newton_raphson( x_2d, &kernel) ;

  W = x_2d[0];
  vsq = x_2d[1];
        
  /* Problem with solver, so return denoting error before doing anything further */
  if( ((retval) != 0.) || (W == FAIL_VAL) ) {
   fprintf(stdout," *dens        = %26.16e \n", cons_pred->D    );
   fprintf(stdout," *sx          = %26.16e \n", cons_pred->S[0] );    
   fprintf(stdout," *sy          = %26.16e \n", cons_pred->S[1] );    
   fprintf(stdout," *sz          = %26.16e \n", cons_pred->S[2] );    
   fprintf(stdout," *tau         = %26.16e \n", cons_pred->tau  );    
   fprintf(stdout," *rho         = %26.16e \n", prims->rho0    );    
   fprintf(stdout," *velx        = %26.16e \n", prims->v[0]    );    
   fprintf(stdout," *vely        = %26.16e \n", prims->v[1]    );    
   fprintf(stdout," *velz        = %26.16e \n", prims->v[2]    );    
   fprintf(stdout," *epsilon     = %26.16e \n", prims->u       );
   fprintf(stdout," *pressure    = %26.16e \n", prims->p      );
#ifdef MAGNETIC
  fprintf(stdout," *Bx          = %26.16e \n", prims->B[0]    );    
  fprintf(stdout," *By          = %26.16e \n", prims->B[1]    );    
  fprintf(stdout," *Bz          = %26.16e \n", prims->B[2]    );    
  fprintf(stdout," *bsq         = %26.16e \n", prims->bsq     );    
#endif

   printf("My vsq is %2.12g %2.12g - gamma %2.12g\n",vsq,W,gamma);
   printf("Initial guess: local %2.12g %2.12g - input %2.12g %2.12g - w %2.12g\n",rho0,p,prims->rho0,prims->p,rhoh);
   printf("My metric %g %g %g %g %g %g\n",gr->gamma[0],gr->gamma[1],gr->gamma[2],gr->gamma[3],gr->gamma[4],gr->gamma[5]);
   return -(retval*100.+1.);
  }
  else{
    if(W <= 0. || W > W_TOO_BIG) {
      return -3.;
    }
  }

  // Calculate v^2:
  if( vsq >= 1. ) {
    return -4.;
  }

  // Recover the primitive variables from the scalars and conserved variables:
  if(gr_limiter==1){
  	long double RMAX=1.0-1./GR_GAMMA_LIMITER/GR_GAMMA_LIMITER;
	vsq=0.5*(RMAX+vsq-fabs(vsq-RMAX));}

  gamma = 1./sqrt(1. - vsq);

  rho0 = kernel.D / gamma;

  rhoh = W * (1. - vsq) ;

  //Ensures that rhoh is at least rho0
  rhoh = fmax(rho0,rhoh);

  p = gr_pressure_rho0_w(rho0,rhoh) ;  // EOS
  //Let's trap this case (pressurless gas)

  u = rhoh - (rho0 + p) ;

  // User may want to handle this case differently, e.g. do NOT return upon 
  // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:
#ifndef EOS_GENERAL
  if(p<=0 && fabs(p)/rhoh<1e-10) {u = prims->u*rho0; p=gr_pressure_rho0_u(rho0,u/rho0);}
  if(u<=0) u = p/(GAMMA_DEFAULT-1);
#endif
  if( (rho0 <= 0.) || (u <= 0.) ) { 
   printf("My guess was %g %g %g %g %g\n",prims->rho0,prims->v[0],prims->v[1],prims->v[2],prims->u);
  }

  prims->rho0 = rho0;
  prims->u   = u / rho0;
  prims->p   = p ;

  double W_inv = 1./(W+kernel.Bsq);
#ifdef MAGNETIC
  QdB_o_W  = QdotB / W;
  prims->bsq = kernel.Bsq * (1.-vsq) + QdB_o_W*QdB_o_W;
#endif
  //W=rho*h*Wlorentz^2 in our notation
  prims->v[0] = W_inv * ( uSvec[0] );
  prims->v[1] = W_inv * ( uSvec[1] );
  prims->v[2] = W_inv * ( uSvec[2] );
#ifdef MAGNETIC
  prims->v[0] += W_inv * ( QdB_o_W*prims->B[0]);
  prims->v[1] += W_inv * ( QdB_o_W*prims->B[1]);
  prims->v[2] += W_inv * ( QdB_o_W*prims->B[2]);
#endif
  //Limit velocities to GR_GAMMA_LIMITER
  double WW=lorentz_factor(gr->gamma,prims->v, gr_limiter);
  // User may want to handle this case differently, e.g. do NOT return upon 
  // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:
  if( (rho0 <= 0.) || (u <= 0.) ) { 
    printf("Original %g %g %g sqrtg %g %g\n",cons_pred->D,cons_pred->S[0],cons_pred->tau,inv_sqrtg,tau);
    printf("My terms are rhoh %g rho0 %g p %g\n",rhoh,rho0,gr_pressure_rho0_w(rho0,rhoh));
    printf("res rho %g u %g p %g vsq %g- WL %g S %g %g %g tau %g\n",rho0,u,p,vsq,gamma,Svec[0],Svec[1],Svec[2],tau);
    return -1; //rho/u negative
  }

#if(DEBUG_CON2PRIMM)
  fprintf(stdout,"rho          = %26.16e \n",prims->rho0);
  fprintf(stdout,"epsilon      = %26.16e \n",prims->u   );
  fprintf(stdout,"pressure     = %26.16e \n",prims->p   );
  fprintf(stdout,"Wlorentz     = %26.16e \n",(double)gamma      );
  fprintf(stdout,"velx         = %26.16e \n",prims->v[0]);
  fprintf(stdout,"vely         = %26.16e \n",prims->v[1]);
  fprintf(stdout,"velz         = %26.16e \n",prims->v[2]);
#ifdef MAGNETIC
  fprintf(stdout,"bsq          = %26.16e \n",prims->bsq  );
#endif
  fflush(stdout);
  getchar();
#endif

  /* done! Return the Lorentz factor*/
  return gamma;

}


/**********************************************************************/ 
/****************************************************************************
   vsq_calc(): 
    
      -- evaluate v^2 (spatial, normalized velocity) from 
            W = \gamma^2 w 

****************************************************************************/
static double vsq_calc(double W, struct kernel_NR *kernel)
{
        double Wsq,Xsq,Bsq_W;
        
        Wsq = W*W ;
        Bsq_W = (kernel->Bsq + W);
        Xsq = Bsq_W * Bsq_W;

        return(  ( Wsq * kernel->Qtsq  + kernel->QdotBsq * (Bsq_W + W)) / (Wsq*Xsq) );
}


/********************************************************************

  x1_of_x0(): 
           
    -- calculates v^2 from W  with some physical bounds checking;
    -- asumes x0 is already physical
    -- makes v^2 physical  if not;

*********************************************************************/

static double x1_of_x0(double x0, struct kernel_NR *kernel ) 
{
  double vsq;
  double dv = MIN_DV;

  vsq = fabs(vsq_calc(x0,kernel)) ; // guaranteed to be positive 

  return( ( vsq > 1. ) ? (1.0 - dv) : vsq   ); 

}

/********************************************************************

  validate_x(): 
           
    -- makes sure that x[0,1] have physical values, based upon 
       their definitions:
    
*********************************************************************/

static void validate_x(double x[2], double x0[2] ) 
{
  
  const double dv = MIN_DV;

  /* Always take the absolute value of x[0] and check to see if it's too big:  */ 
  x[0] = fabs(x[0]);
  x[0] = (x[0] > W_TOO_BIG) ?  x0[0] : x[0];
  
  x[1] = (x[1] < 0.) ?   0.       : x[1];  /* if it's too small */
  x[1] = (x[1] > 1.) ?  (1. - dv) : x[1];  /* if it's too big   */

  return;

}

/************************************************************

  twod_newton_raphson(): 

    -- performs Newton-Rapshon method on an 2d system.

    -- inspired in part by Num. Rec.'s routine newt();

*****************************************************************/
static int twod_newton_raphson( double x[], struct kernel_NR *kernel)
{
  double f, df, dx[2], x_old[2];
  double resid[2], jac[2][2];
  double errx, x_orig[2];
  int    n_iter, i_extra, doing_extra;
  double vsq,W,W_old;
  const double dv = (1.-MIN_DV);

  int   keep_iterating;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;
  x_old[0] = x_orig[0] = x[0] ;
  x_old[1] = x_orig[1] = x[1] ;

  vsq = W = W_old = 0.;
  n_iter = 0;

  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    func_vsq(x, dx, resid, jac, &f, &df, kernel);  /* returns with new dx, f, df */
      
    if(isnan(dx[0]+dx[1])) printf("My x %g %g %g %g\n",x[0],x[1],dx[0],dx[1]);
    /* Save old values before calculating the new: */
    errx = 0.;
    x_old[0] = x[0] ; 
    x_old[1] = x[1] ; 

    /* Make the newton step: */
    x[0] += dx[0]  ; 
    x[1] += dx[1]  ; 

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    if( x[0] < 0. ) {  x[0] = fabs(x[0]);  } 
    else { 
     if(x[0] > W_TOO_BIG)  { x[0] = x_old[0] ; }
    }

    if( x[1] < 0. ) {  x[1] = 0.; } 
    else { 
      if( x[1] > 1. ) { x[1] = dv; }
    }
      
    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    
    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0)) 
        || (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)


    /*  Check for bad untrapped divergences : */
  if( (!isfinite(f)) ||  (!isfinite(df)) ) {
    return(2);
  }


  if( fabs(errx) <= NEWT_TOL ){
    return(0);
  }
  else if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    return(0);
  }
  else {
    return(1);
  } 

  

}



/**********************************************************************/
/*********************************************************************************
   func_vsq(): 

        -- calculates the residuals, and Newton step for general_newton_raphson();
        -- for this method, x=W,vsq here;

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
 *********************************************************************************/

static void func_vsq(double x[], double dx[], double resid[], 
                     double jac[][2], double *f, double *df, struct kernel_NR *kernel)
{

  
  double  W, vsq, p_tmp, dPdvsq, dPdW;
  double t11, t16,t18,t2,t21,   t23,   t24,   t25,   t3,   t35,  t4,   t40;


  W = x[0];
  vsq = x[1];
  
  p_tmp = gr_eos_info(W, vsq, &dPdW, &dPdvsq, kernel);

  // These expressions were calculated using Mathematica, but made into efficient 
  // code using Maple.  Since we know the analytic form of the equations, we can 
  // explicitly calculate the Newton-Raphson step:

  t2 = -kernel->half_Bsq+dPdvsq;
  t3 = kernel->Bsq+W;
  t4 = t3*t3;
  t23 = 1./W;
  t16 = kernel->QdotBsq*t23*t23;
  t11 = kernel->Qtsq-vsq*t4+t16*(kernel->Bsq+W+W);
  t18 = -kernel->Qdotn-kernel->half_Bsq*(1.0+vsq)+0.5*t16-W+p_tmp;
  t24 = t16*t23;
  t25 = -1.0+dPdW-t24;
  t35 = t25*t3+(kernel->Bsq-2.0*dPdvsq)*(t16+vsq*W)*t23;
  //  t21 = 1/t3;
  //  t36 = 1/t35;
  t21 = 1/(t3*t35);
  dx[0] = -(t2*t11+t4*t18)*t21;
  t40 = -2*(vsq+t24)*t3;
  dx[1] = -(-t25*t11+t40*t18)*t21;
  //  detJ = t3*t35;
  jac[0][0] = t40;
  jac[0][1] = -t4;
  jac[1][0] = t25;
  jac[1][1] = t2;
  resid[0] = t11;
  resid[1] = t18;

  *df = -resid[0]*resid[0] - resid[1]*resid[1];

  *f = -0.5 * ( *df );

}

/****************************************************************************** 
             END   
 ******************************************************************************/


#undef DEBUG_CON2PRIMM
#endif
