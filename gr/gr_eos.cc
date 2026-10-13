#include "../declarations/allvars.h"
#include "../core/proto.h"
#ifdef GENERAL_RELATIVITY
/********************************************************************** 
 ********************************************************************** 

 The following routines specify the equation of state. If the user
 wishes to use another equation of state, the below functions must be
 replaced by equivalent routines based upon the new EOS.

 **********************************************************************
**********************************************************************/
/* primitive density as a function of D and v */
double gr_get_primitive_density(double density,double velocity[],struct gr_metric *coeffs, int gr_limiter)
{
        double vel_p[3];
        int k;
	for(k=0;k<3;k++) vel_p[k]=velocity[k]; 
        Vphys(vel_p,coeffs->beta,coeffs->alpha);
        return density/lorentz_factor(coeffs->gamma,vel_p, gr_limiter)/coeffs->gdet;
}


/* pressure as a function of rho0 and u */
double gr_pressure_rho0_u(double rho0, double u)
{
#ifdef GR_TEST_PARTICLE
	return 0;
#else
	return (GAMMA_DEFAULT -1.)*rho0*u;
#endif
}

/* Pressure as a function of rho0 and w = rho0 + rho0*u + p */
double gr_pressure_rho0_w(double rho0, double w)
{
#ifdef GR_TEST_PARTICLE
	return 0;
#else
	//if P=(gamma-1)*rho*u, then h= rho + gamma*rho*u = rho + gamma/(gamma-1)*P
	return (GAMMA_DEFAULT-1.)*(w - rho0)/GAMMA_DEFAULT ;
#endif
}

/* Pressure as a function of rho0 and w = rho0 + rho0*u + p */
double gr_get_u_from_p(double rho0, double p)
{
#ifdef GR_TEST_PARTICLE
	return 0;
#else
	//if P=(gamma-1)*rho*u, then h= rho + gamma*rho*u = rho + gamma/(gamma-1)*P
	return p/rho0/(GAMMA_DEFAULT-1.);
#endif
}


/**********************************************************************/
/********************************************************************** 
  eos_info():
 
      -- returns with all the EOS-related values needed;
 **********************************************************************/
double gr_eos_info(double W, double vsq, double *dpdw, double *dpdvsq, struct kernel_NR* kernel)
{
  double ftmp,gtmp;

  ftmp = 1. - vsq;
  gtmp = sqrt(ftmp);

  double gam_m1_o_gam = ((GAMMA_DEFAULT-1.)/GAMMA_DEFAULT);

  *dpdw =  gam_m1_o_gam * ftmp ;
  *dpdvsq =  gam_m1_o_gam * ( 0.5 * kernel->D/gtmp  -  W ) ;

  return( gam_m1_o_gam * ( W * ftmp  -  kernel->D * gtmp )  );  // p 

}


#endif
