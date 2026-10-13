#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "string.h"
#ifdef GENERAL_RELATIVITY
//#define USE_BOYER_LINQUIST_COORDINATES

#if (GR_METRIC<0)
#include "custom_metric.h"
#elif (GR_METRIC==2)
#include "SuperposedBBH.h"
#endif

void get_covariant_metric(double* pos,double gcov[])
{
  double dp2[3],sinth2,dp[3];
  double R2=0.;
  double r2;
  int i;
  for(i=0;i<3;i++) dp[i]=pos[i];
#ifdef BOX_PERIODIC //We assume in this case that the source is at the centre of the box
  dp[0]-=boxHalf_X;
  dp[1]-=boxHalf_Y;
  dp[2]-=boxHalf_Z;
  NEAREST_XYZ(dp[0],dp[1],dp[2],1);
#endif

  for(i=0;i<3;i++) {dp2[i]=dp[i]*dp[i]; R2+=dp2[i]; }

#if (GR_METRIC<0) //User-defined metric

  get_custom_metric_coefficients(R2,dp,dp2,gcov);

#elif (GR_METRIC==0) //Minkowski flat space-time

  gcov[0]=-1;
  gcov[1]=gcov[2]=gcov[3]=gcov[5]=gcov[6]=gcov[8]=0.;
  gcov[4]=gcov[7]=gcov[9]=1.;

#elif (GR_METRIC==1) //Kerr-Schwarzchild metric
  double M=All.MetricParams[0];
  double a=All.MetricParams[1]*M;
  double a2=a*a;

  r2=(R2-a2+sqrt((R2-a2)*(R2-a2)+4*a2*dp2[2]))/2.;

#ifndef USE_KERRSCHILD_COORDINATES
  sinth2=1.-dp2[2]/r2;

  double rho2=r2+a2*dp2[2]/r2;
  double delta=r2-2*M*sqrt(r2)+a2;
  double Att=2.0*M*sqrt(r2)/rho2;
  //Now compute Gcov
  gcov[0]=Att-1.0; //t-t

  double gtp=-Att*a*sinth2; //t-phi
  double gpp=(r2+a2+Att*a2*sinth2)*sinth2;//phi-phi

  //Now the more complex terms
  //t-x,x-t
  gcov[1]=-dp[1]/(dp2[0]+dp2[1])*gtp;
  //t-y,y-t
  gcov[2]= dp[0]/(dp2[0]+dp2[1])*gtp;
  //t-z,z-t
  gcov[3]=0.;

  double fac1=r2/rho2/delta;
  double fac2=gpp/((dp2[0]+dp2[1])*(dp2[0]+dp2[1]));
  double fac3=dp2[2]/rho2/(r2*sinth2);
  double fac4=(1-dp2[2]*(r2+a2)/r2/rho2);
  //x-x
  gcov[4]=fac1*dp2[0] + fac2*dp2[1] + fac3*dp2[0];
  //y-y
  gcov[7]=fac1*dp2[1] + fac2*dp2[0] + fac3*dp2[1];
  //z-z
  gcov[9]=(1-fac4)*(a2+r2)/delta + rho2/(r2*sinth2)*fac4*fac4;
  //x-y,y-x
  gcov[5]=dp[0]*dp[1]*(fac1 - fac2 + fac3);
  //x-z,z-x
  gcov[6]=(a2+r2)*dp[0]*dp[2]/rho2/delta - dp[0]*dp[2]/(r2*sinth2)*fac4; //(1.-(dp2[2]*(a2+r2)/r2/rho2));
  //y-z,z-y
  gcov[8]=(a2+r2)*dp[1]*dp[2]/rho2/delta - dp[1]*dp[2]/(r2*sinth2)*fac4; //(1.-(dp2[2]*(a2+r2)/r2/rho2));
#else
  double rho2=r2*r2 + a2*dp2[2];
  double r = sqrt(r2);
  double Mfac = 2*M*r*r*r / rho2;
  double denom = r2 + a2;

  //Define the l_mu vector
  double l_cov[4];
  l_cov[0] = 1;
  l_cov[1] = (r*dp[0] + a*dp[1])/denom;
  l_cov[2] = (r*dp[1] - a*dp[0])/denom;
  l_cov[3] = dp[2]/r;

  //Define the Minkowski metric
  double eta[10];
  eta[0]=-1;
  eta[1]=eta[2]=eta[3]=eta[5]=eta[6]=eta[8]=0.;
  eta[4]=eta[7]=eta[9]=1.;

  //Finally compute the Kerr-Schild covariant metric
  gcov[0] = eta[0] + Mfac*l_cov[0]*l_cov[0]; //t-t
  gcov[1] = eta[1] + Mfac*l_cov[0]*l_cov[1]; //t-x, x-t
  gcov[2] = eta[2] + Mfac*l_cov[0]*l_cov[2]; //t-y, y-t
  gcov[3] = eta[3] + Mfac*l_cov[0]*l_cov[3]; //t-z, z-t

  gcov[4] = eta[4] + Mfac*l_cov[1]*l_cov[1]; //x-x
  gcov[5] = eta[5] + Mfac*l_cov[1]*l_cov[2]; //x-y, y-x
  gcov[6] = eta[6] + Mfac*l_cov[1]*l_cov[3]; //x-z, z-x

  gcov[7] = eta[7] + Mfac*l_cov[2]*l_cov[2]; //y-y
  gcov[8] = eta[8] + Mfac*l_cov[2]*l_cov[3]; //y-z, z-y

  gcov[9] = eta[9] + Mfac*l_cov[3]*l_cov[3]; //z-z

#endif
#elif (GR_METRIC==2)
  double traj_array[20]={grP[0].Pos[0],grP[0].Pos[1],grP[0].Pos[2],grP[1].Pos[0],grP[1].Pos[1],grP[1].Pos[2],
 			 grP[0].Vel[0],grP[0].Vel[1],grP[0].Vel[2],grP[1].Vel[0],grP[1].Vel[1],grP[1].Vel[2],
 			 grP[0].Spin[0],grP[0].Spin[1],grP[0].Spin[2],grP[1].Spin[0],grP[1].Spin[1],grP[1].Spin[2],
 			 grP[0].Mass,grP[1].Mass};
  SuperposedBBH(All.Time,dp, gcov, traj_array);
#endif

}

/* -------- DO NOT MODIFY HERE BELOW UNLESS YOU KNOW WHAT YOU ARE DOING -------- */
/*    Utility function to determine metric coefficients and their derivatives    */

double spatial_metric_determinant(double gamma[])
{
        //xx-xy-xz-yy-yz-zz
        return (-gamma[2]*gamma[2]*gamma[3] + 2.0*gamma[1]*gamma[2]*gamma[4] - gamma[0]*gamma[4]*gamma[4] \
                -gamma[1]*gamma[1]*gamma[5] + gamma[0]*gamma[3]*gamma[5]);
}


void get_upper_spatial_metric(double gcov[6],double gcon[6])
{
  double invdet = 1.e0 / spatial_metric_determinant(gcov);
  //xx
  gcon[0] =(-gcov[4]*gcov[4] + gcov[3]*gcov[5])*invdet;
  //xy
  gcon[1] =( gcov[2]*gcov[4] - gcov[1]*gcov[5])*invdet;
  //xz
  gcon[2] =(-gcov[2]*gcov[3] + gcov[1]*gcov[4])*invdet;
  //yy
  gcon[3] =(-gcov[2]*gcov[2] + gcov[0]*gcov[5])*invdet;
  //yz
  gcon[4] =( gcov[1]*gcov[2] - gcov[0]*gcov[4])*invdet;
  //zz
  gcon[5] =(-gcov[1]*gcov[1] + gcov[0]*gcov[3])*invdet;

}

void get_extra_metric_coefficients(double *gcov,double *gcon,double beta[3],double *alpha)
{

  get_upper_spatial_metric(gcov+4,gcon);

  //Now compute beta contravariant (beta_i = g0i)
  double beta2=0;
  double betacov[3]={gcov[1],gcov[2],gcov[3]};
  raise_lower_indices(gcon,betacov,beta);

  //Now compute alpha
  beta2=beta[0]*gcov[1]+beta[1]*gcov[2]+beta[2]*gcov[3];
  *alpha=sqrt(beta2-gcov[0]);
}


void set_metric_coefficients(double *pos, double potential, struct gr_metric *gr)
{

  double gcov[10];
  get_covariant_metric(pos,gcov);

  //Copy gamma
  int k;
  for(k=0;k<6;k++) gr->gamma[k]=gcov[k+4];

  //Get contravariant metric, lapse and shift
  get_extra_metric_coefficients(gcov,gr->ugamma,gr->beta,&(gr->alpha));

  //Finally get the sqrt(gamma)
  gr->gdet = sqrt(spatial_metric_determinant(gr->gamma));

  //If Newtonian self-gravity is included, we add it as a correction to the tt-term only (h_00)
#ifndef SELFGRAVITY_OFF
  gr->alpha = sqrt(gr->alpha*gr->alpha+2*potential);
#endif
}

#ifndef USE_BOYER_LINQUIST_COORDINATES
void compute_metric_jacobian(double* dp,double deriv[][4][4],double dalpha[],double dbeta[][3])
{
     int dir,j,k;
     //H = x+dx, L=x-dx//
     double dpH[3],dpL[3];
     //Metric coefficients at H and L
     double gcovH[10],gcovL[10];
     double alphaH,alphaL;
     double betaH[3],betaL[3];
     //Temporary variables to store the unpacked metric used for derivatives
     double gtmpH[4][4],gtmpL[4][4];
     //Dummy array to store the contravariant metric
     double gcon[6];


     for(dir=0;dir<3;dir++) dpH[dir]=dpL[dir]=dp[dir];

     for(dir=0;dir<3;dir++)
     {
	dpH[dir]+=GR_EPS_TOLERANCE;
	dpL[dir]-=GR_EPS_TOLERANCE;
	//Compute metric at L and H along direction 'dir'
	get_covariant_metric(dpH,gcovH);
	get_covariant_metric(dpL,gcovL);
	//Compute contravariant metric coefficients at H and L (lapse and shift)
	get_extra_metric_coefficients(gcovH,gcon,betaH,&alphaH);
	get_extra_metric_coefficients(gcovL,gcon,betaL,&alphaL);

	//Unpack metric tensor
	unpack_4tensor(gcovH,gtmpH);
	unpack_4tensor(gcovL,gtmpL);
	//Compute derivatives
	for(j=0;j<4;j++)
	{
	 deriv[dir][j][j]=(0.5*(gtmpH[j][j]-gtmpL[j][j])/GR_EPS_TOLERANCE);
	 for(k=j+1;k<4;k++)
	 {
	   deriv[dir][j][k]=deriv[dir][k][j]=(0.5*(gtmpH[j][k]-gtmpL[j][k])/GR_EPS_TOLERANCE);
	 }
	}
	dalpha[dir] = 0.5*(alphaH-alphaL)/GR_EPS_TOLERANCE;
//#ifdef DIVBCLEANING_GR
//	double gdetH,gdetL;
//	gdetH=sqrt(spatial_metric_determinant(gcovH+4));
//	gdetL=sqrt(spatial_metric_determinant(gcovL+4));
//	dalpha[3+dir] = 0.5*(alphaH/gdetH-alphaL/gdetL);
//#endif
	
	for(j=0;j<3;j++)  dbeta[dir][j] = 0.5*(betaH[j]-betaL[j])/GR_EPS_TOLERANCE;

	//Reset coordinates along this direction
	dpH[dir]=dpL[dir]=dp[dir];
     }
}

#else
/* --- This is the same as above, but for the 3 directions used for the derivatives, which correspond to BL coordinates
The projection in cartesian coordinates is done at the end of the routine
Note that this works only for single BHs, since the coordinates depend on the spin of the central object --- */
void compute_metric_jacobian_bl(double* dp,double deriv[][4][4],double dalpha[],double dbeta[][3])
{
     double a=All.MetricParams[1]; //Copy the spin into 'a'
     int dir,j,k;
     //H = x+dx, L=x-dx//
     double dpH[3],dpL[3];
     //Metric coefficients at H and L
     double gcovH[10],gcovL[10];
     double alphaH,alphaL;
     double betaH[3],betaL[3];

     //Temporary variables to store the unpacked metric used for derivatives
     double gtmpH[4][4],gtmpL[4][4];
     //Temporary variable to store the contravariant metric
     double gcon[6];
     //Temporary variables used for conversion into cartesian coordinates
     double d_tmp[3][4][4],d_betatmp[3][3];
//#ifdef DIVBCLEANING_GR
//     double d_alphatmp[6];
//#else
     double d_alphatmp[3];
//#endif

     for(dir=0;dir<3;dir++) dpH[dir]=dpL[dir]=dp[dir];
     //Instead of x,y,z, we use r,theta,phi (BL coordinates)
     double r3d=sqrt(dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2]);
     double r=sqrt((r3d*r3d-a*a+sqrt((r3d*r3d-a*a)*(r3d*r3d-a*a)+4*a*a*dp[2]*dp[2]))/2.);
     double sinth2=1-dp[2]*dp[2]/r/r;
     double theta=atan2(sqrt(sinth2),sqrt(1-sinth2));
     double phi=atan2(dp[1]/sqrt((r*r+a*a)*sinth2),dp[0]/sqrt((r*r+a*a)*sinth2));
     double rho2=r*r+a*a*(1-sinth2);

     /* --- Conversion factors to transform the derivatives into cartesian coordinates --- */
     double d_dx[3][3]={{r*dp[0]/rho2,dp[0]*dp[2]/rho2/(r*sin(theta)),-dp[1]/(dp[0]*dp[0]+dp[1]*dp[1])},
     	{r*dp[1]/rho2,dp[1]*dp[2]/rho2/(r*sin(theta)), dp[0]/(dp[0]*dp[0]+dp[1]*dp[1])},
	{dp[2]*(a*a+r*r)/r/rho2,-(1-dp[2]*dp[2]*(a*a+r*r)/r/r/rho2)/(r*sin(theta)),0.0}};
     /* --- */

     for(dir=0;dir<3;dir++)
     {
	//Reset coordinates along this direction
	for(j=0;j<3;j++) dpH[j]=dpL[j]=dp[j];
	if(dir==0)
	{
	  dpH[0]=sqrt((r+GR_EPS_TOLERANCE)*(r+GR_EPS_TOLERANCE)+a*a)*sin(theta)*cos(phi);
	  dpH[1]=sqrt((r+GR_EPS_TOLERANCE)*(r+GR_EPS_TOLERANCE)+a*a)*sin(theta)*sin(phi);
	  dpH[2]=(r+GR_EPS_TOLERANCE)*cos(theta);
	  dpL[0]=sqrt((r-GR_EPS_TOLERANCE)*(r-GR_EPS_TOLERANCE)+a*a)*sin(theta)*cos(phi);
	  dpL[1]=sqrt((r-GR_EPS_TOLERANCE)*(r-GR_EPS_TOLERANCE)+a*a)*sin(theta)*sin(phi);
	  dpL[2]=(r-GR_EPS_TOLERANCE)*cos(theta);
	}
	else if(dir==1)
	{
	  dpH[0]=sqrt(r*r+a*a)*sin(theta+GR_EPS_TOLERANCE)*cos(phi);
	  dpH[1]=sqrt(r*r+a*a)*sin(theta+GR_EPS_TOLERANCE)*sin(phi);
	  dpH[2]=r*cos(theta+GR_EPS_TOLERANCE);
	  dpL[0]=sqrt(r*r+a*a)*sin(theta-GR_EPS_TOLERANCE)*cos(phi);
	  dpL[1]=sqrt(r*r+a*a)*sin(theta-GR_EPS_TOLERANCE)*sin(phi);
	  dpL[2]=r*cos(theta-GR_EPS_TOLERANCE);
	}
	else
	{
	  dpH[0]=sqrt(r*r+a*a)*sin(theta)*cos(phi+GR_EPS_TOLERANCE);
	  dpH[1]=sqrt(r*r+a*a)*sin(theta)*sin(phi+GR_EPS_TOLERANCE);
	  dpL[0]=sqrt(r*r+a*a)*sin(theta)*cos(phi-GR_EPS_TOLERANCE);
	  dpL[1]=sqrt(r*r+a*a)*sin(theta)*sin(phi-GR_EPS_TOLERANCE);
	}
	//Compute metric at H and L alond direction 'dir'
	get_covariant_metric(dpH,gcovH);
	get_covariant_metric(dpL,gcovL);
	//Determine the upper metric coefficients
	compute_extra_metric_coefficients(gcovH,gcon,betaH,&alphaH);
	compute_extra_metric_coefficients(gcovL,gcon,betaL,&alphaL);

	unpack_4tensor(gcovH,gtmpH);
	unpack_4tensor(gcovL,gtmpL);
	for(j=0;j<4;j++)
	{
	 d_tmp[dir][j][j]=(0.5*(gtmpH[j][j]-gtmpL[j][j])/GR_EPS_TOLERANCE);
	 for(k=j+1;k<4;k++)
	 {
	   d_tmp[dir][j][k]=d_tmp[dir][k][j]=(0.5*(gtmpH[j][k]-gtmpL[j][k])/GR_EPS_TOLERANCE);
	 }
	}

	d_alphatmp[dir] = 0.5*(alphaH-alphaL)/GR_EPS_TOLERANCE;
//#ifdef DIVBCLEANING_GR
//	double gdetH,gdetL;
//	gdetH=sqrt(spatial_metric_determinant(gcovH+4));
//	gdetL=sqrt(spatial_metric_determinant(gcovL+4));
//	d_alphatmp[3+dir] = 0.5*(alphaH/gdetH-alphaL/gdetL);
//#endif

	for(j=0;j<3;j++) d_betatmp[dir][j] = 0.5*(betaH[j]-betaL[j])/GR_EPS_TOLERANCE;
     }
     for(dir=0;dir<3;dir++)
     {
	for(j=0;j<4;j++)
	{
	  for(k=j;k<4;k++){ 
	  deriv[dir][j][k]=deriv[dir][k][j] = d_tmp[0][j][k]*d_dx[dir][0]+d_tmp[1][j][k]*d_dx[dir][1]+d_tmp[2][j][k]*d_dx[dir][2];
	  }
	}
	dalpha[dir] = d_alphatmp[0]*d_dx[dir][0]+d_alphatmp[1]*d_dx[dir][1]+d_alphatmp[2]*d_dx[dir][2];
//#ifdef DIVBCLEANING_GR
//	dalpha[3+dir] = d_alphatmp[3]*d_dx[dir][0]+d_alphatmp[4]*d_dx[dir][1]+d_alphatmp[5]*d_dx[dir][2];
//#endif
	for(j=0;j<3;j++) dbeta[dir][j] = d_betatmp[0][j]*d_dx[dir][0]+d_betatmp[1][j]*d_dx[dir][1]+d_betatmp[2][j]*d_dx[dir][2];
     }
}
#endif

void get_metric_derivatives(int p,double Kij[][3],double d_gij[][4][4],double dalpha[],double dbeta[][3])
{
     int i,j,l,w;
     double invalpha2,gmunu[4][4];

     invalpha2=1./(CellP[p].gr.alpha*CellP[p].gr.alpha);
     gmunu[0][0]=-invalpha2;
     gmunu[0][1]=gmunu[1][0]=CellP[p].gr.beta[0]*invalpha2;
     gmunu[0][2]=gmunu[2][0]=CellP[p].gr.beta[1]*invalpha2;
     gmunu[0][3]=gmunu[3][0]=CellP[p].gr.beta[2]*invalpha2;

     gmunu[1][1]=CellP[p].gr.ugamma[0] - CellP[p].gr.beta[0]*CellP[p].gr.beta[0]*invalpha2;
     gmunu[2][2]=CellP[p].gr.ugamma[3] - CellP[p].gr.beta[1]*CellP[p].gr.beta[1]*invalpha2;
     gmunu[3][3]=CellP[p].gr.ugamma[5] - CellP[p].gr.beta[2]*CellP[p].gr.beta[2]*invalpha2;
     gmunu[1][2]=gmunu[2][1] = CellP[p].gr.ugamma[1] - CellP[p].gr.beta[0]*CellP[p].gr.beta[1]*invalpha2;
     gmunu[1][3]=gmunu[3][1] = CellP[p].gr.ugamma[2] - CellP[p].gr.beta[0]*CellP[p].gr.beta[2]*invalpha2;
     gmunu[2][3]=gmunu[3][2] = CellP[p].gr.ugamma[4] - CellP[p].gr.beta[1]*CellP[p].gr.beta[2]*invalpha2;

#ifdef USE_BOYER_LINQUIST_COORDINATES
     compute_metric_jacobian_bl(P[p].Pos,d_gij,dalpha,dbeta);
#else
     compute_metric_jacobian(P[p].Pos, d_gij,dalpha,dbeta);
#endif

     //Now rearrange terms to get the Christoffel symbol Gamma^0_ij
     //0.5*g^0mu*(d_j g_mui + d_i g_muj - d_mu g_ij)
     for(i=0;i<3;i++)
        for(j=0;j<3;j++)
	{
	   //First the 0 term
	   Kij[i][j]= 0.5*gmunu[0][0]*(d_gij[j][0][i+1]+d_gij[i][0][j+1]); //the time derivative is null d_gij[w][i][j]);
	   for(w=1;w<4;w++){ Kij[i][j]+= 0.5*gmunu[0][w]*(d_gij[j][w][i+1]+d_gij[i][w][j+1]-d_gij[w-1][i+1][j+1]);
	   }
	   //Before concluding, apply -alpha
	   Kij[i][j]*=-CellP[p].gr.alpha;
	}
}


#endif
