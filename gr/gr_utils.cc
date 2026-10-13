#include "../declarations/allvars.h"
#include "../core/proto.h"
#include <sys/stat.h>
#ifdef GENERAL_RELATIVITY
#ifdef GR_USECBWAVES
double **cbwaves;
double dt_cbwaves;

void load_cbwaves()
{
  int Ncbwaves=0;
  int i,j;
  FILE* fin=fopen("cbwaves.txt","r");
  if(fin==NULL) endrun(57373);
  fscanf(fin,"%d",&Ncbwaves);
  cbwaves=(double**)malloc(sizeof(double*)*Ncbwaves);
  for(i=0;i<Ncbwaves;i++)
  {
    cbwaves[i]=(double*)malloc(sizeof(double)*13);
    for(j=0;j<13;j++)
    fscanf(fin,"%lf",cbwaves[i]+j);
  }
  fclose(fin);
  dt_cbwaves=cbwaves[1][0]-cbwaves[0][0];
  if(ThisTask==0) printf("CBwaves trajectory loaded: dt %g\n",dt_cbwaves);
}
#endif

#if (GR_METRIC==2)
void gr_update_binary_position()
{
  int i,ibh,ispin;
#ifdef GR_USECBWAVES
  int itime=(int)(All.Time/dt_cbwaves);
  double t_frac=(All.Time-cbwaves[itime][0])/(cbwaves[itime+1][0]-cbwaves[itime][0]);
  grP[0].Pos[0]=cbwaves[itime+1][1]*t_frac + cbwaves[itime][1]*(1-t_frac);
  grP[0].Pos[1]=cbwaves[itime+1][2]*t_frac + cbwaves[itime][2]*(1-t_frac);
  grP[0].Pos[2]=0;//cbwaves[itime][3];
  grP[1].Pos[0]=cbwaves[itime+1][4]*t_frac + cbwaves[itime][4]*(1-t_frac);
  grP[1].Pos[1]=cbwaves[itime+1][5]*t_frac + cbwaves[itime][5]*(1-t_frac);
  grP[1].Pos[2]=0;//cbwaves[itime][6];

  grP[0].Vel[0]=cbwaves[itime+1][7]*t_frac + cbwaves[itime][7]*(1-t_frac);
  grP[0].Vel[1]=cbwaves[itime+1][8]*t_frac + cbwaves[itime][8]*(1-t_frac);
  grP[0].Vel[2]=0;//cbwaves[itime][9];
  grP[1].Vel[0]=cbwaves[itime+1][10]*t_frac + cbwaves[itime][10]*(1-t_frac);
  grP[1].Vel[1]=cbwaves[itime+1][11]*t_frac + cbwaves[itime][11]*(1-t_frac);
  grP[1].Vel[2]=0;//cbwaves[itime][12];

  grP[0].Mass   =0.5;
  grP[1].Mass   =0.5;
  grP[0].Spin[0]=grP[0].Spin[1]=grP[0].Spin[2]=0;
  grP[1].Spin[0]=grP[1].Spin[1]=grP[1].Spin[2]=0;
  if(ThisTask==0) printf("My BHs are in (%g %g) and (%g %g)\n",grP[0].Pos[0],grP[0].Pos[1],grP[1].Pos[0],grP[1].Pos[1]);
#else
  double BHprops[20]={0.},*BHprops_all;
  BHprops_all=(double*)malloc(20*NTask*sizeof(double));
  int found=0;
  if(All.TimeStep<=0)
  {
    int BHIds_loc[2]={0},*BHIds;
    BHIds=(int*)malloc(2*NTask*sizeof(int));
    for(i=0;i<2*NTask;i++) BHIds[i]=0;
    for(i=0;i<NumPart && found<2;i++) {if(P[i].Type!=5) continue; BHIds_loc[found++]=P[i].ID;}
    MPI_Allgather(BHIds_loc,2,MPI_INT,BHIds,2,MPI_INT,MPI_COMM_WORLD);
    //ID-ordered
    found=0;
    for(i=0;i<2*NTask && found<2;i++){
	 if(BHIds[i]<=0) continue;
 	 if(found==0) {grP[0].pID=BHIds[i];}
	 else
	 {
		if(grP[0].pID<BHIds[i]) grP[1].pID=BHIds[i]; else {grP[1].pID=grP[0].pID; grP[0].pID=BHIds[i];}
	 }
	 found++;
    }
    free(BHIds);
  }

  
  found=0; //Reset found
  for(i=0;i<NumPart && found<2;i++)
  {
    if(P[i].Type!=5) continue;
    if(P[i].ID==grP[0].pID) {ibh=0;ispin=0;} else {ibh=10;ispin=3;}
    BHprops[ibh]=P[i].Pos[0];
    BHprops[ibh+1]=P[i].Pos[1];
    BHprops[ibh+2]=P[i].Pos[2];
    BHprops[ibh+3]=P[i].Vel[0];
    BHprops[ibh+4]=P[i].Vel[1];
    BHprops[ibh+5]=P[i].Vel[2];
    BHprops[ibh+6]=P[i].Mass;
    BHprops[ibh+7]=All.MetricSpins[ispin]; //BPP(i).BH_Specific_AngMom[0];
    BHprops[ibh+8]=All.MetricSpins[ispin+1]; //BPP(i).BH_Specific_AngMom[1];
    BHprops[ibh+9]=All.MetricSpins[ispin+2]; //BPP(i).BH_Specific_AngMom[2];
    found++;
  }
  MPI_Allgather(BHprops,20,MPI_DOUBLE,BHprops_all,20,MPI_DOUBLE,MPI_COMM_WORLD);
  found=0;
  for(i=0;i<2*NTask && found<2;i++)
  {
    if(BHprops_all[i*10+6]>0)
    {
	ibh=i%2;
	grP[ibh].Pos[0] =BHprops_all[i*10];
	grP[ibh].Pos[1] =BHprops_all[i*10+1];
	grP[ibh].Pos[2] =BHprops_all[i*10+2];
	grP[ibh].Vel[0] =BHprops_all[i*10+3];
	grP[ibh].Vel[1] =BHprops_all[i*10+4];
	grP[ibh].Vel[2] =BHprops_all[i*10+5];
	grP[ibh].Mass   =BHprops_all[i*10+6];
	grP[ibh].Spin[0]=BHprops_all[i*10+7];
	grP[ibh].Spin[1]=BHprops_all[i*10+8];
	grP[ibh].Spin[2]=BHprops_all[i*10+9];
	found++;
    }
  }
#endif

  if(ThisTask == 0) {mkdir(All.OutputDir, 02755);}
  char path[300];
  sprintf(path,"%s/%s",All.OutputDir,"mdot_s1.dat");
  FILE* fsink1=fopen(path,"a");
  sprintf(path,"%s/%s",All.OutputDir,"mdot_s2.dat");
  FILE* fsink2=fopen(path,"a");
  sprintf(path,"%s/%s",All.OutputDir,"outflow_z.dat");
  FILE* foutflow=fopen(path,"a");
  double r2_0=0,r2_1=0;
  int k;
  //Now assign the distance of each gas resolution element from the BHs
  for(i=0;i<NumPart;i++){
        if(P[i].Type!=0 || P[i].Mass<=0) continue;
        r2_0=sqrt((P[i].Pos[0]-grP[0].Pos[0])*(P[i].Pos[0]-grP[0].Pos[0])+(P[i].Pos[1]-grP[0].Pos[1])*(P[i].Pos[1]-grP[0].Pos[1])+
                (P[i].Pos[2]-grP[0].Pos[2])*(P[i].Pos[2]-grP[0].Pos[2]));
        r2_1=sqrt((P[i].Pos[0]-grP[1].Pos[0])*(P[i].Pos[0]-grP[1].Pos[0])+(P[i].Pos[1]-grP[1].Pos[1])*(P[i].Pos[1]-grP[1].Pos[1])+
                (P[i].Pos[2]-grP[1].Pos[2])*(P[i].Pos[2]-grP[1].Pos[2]));
        if(r2_0/grP[0].Mass<r2_1/grP[1].Mass) ibh=0; else ibh=1;

        for(k=0;k<3;k++) {P[i].Min_xyz_to_Sink[k]=(grP[ibh].Pos[k]-P[i].Pos[k])/grP[ibh].Mass;} //We store the relative position in units of M (so that the excision is applied consistently to any BH)
        P[i].Min_Distance_to_Sink = sqrt(P[i].Min_xyz_to_Sink[0]*P[i].Min_xyz_to_Sink[0]+P[i].Min_xyz_to_Sink[1]*P[i].Min_xyz_to_Sink[1]+P[i].Min_xyz_to_Sink[2]*P[i].Min_xyz_to_Sink[2])*grP[ibh].Mass;

        if(All.TimeStep>0)
        {
        double v_local[3]={CellP[i].VelPred[0],CellP[i].VelPred[1],CellP[i].VelPred[2]};
        Vphys(v_local,CellP[i].gr.beta,CellP[i].gr.alpha);
        int gr_limiter=0;
	double dpl[3]; //for(k=0;k<3;k++) {dpl[k]=P[i].Pos[k];}
	for(k=0;k<3;k++) {dpl[k] = -P[i].Min_xyz_to_Sink[k];}
	if(dpl[0]*dpl[0]+dpl[1]*dpl[1]+dpl[2]*dpl[2] <= pow(GAMMA_LIMITER_RADIUS,2)) gr_limiter=1;
//      if(P[i].Pos[0]<LEFT_GAMMA_BOUNDARY||P[i].Pos[0]>RIGHT_GAMMA_BOUNDARY) gr_limiter=1;
        double Wlorentz=lorentz_factor(CellP[i].gr.gamma,v_local,gr_limiter);
        double rest_mass=P[i].Mass/Wlorentz/CellP[i].gr.gdet;
//        if(ibh==0 && r2_0<2*grP[0].Mass){
        if(ibh==0 && r2_0<EHradius(grP[0].Mass,grP[0].Spin)){
                fprintf(fsink1,"%1.10g\t %d %d %g %g %g\n",All.Time,P[i].ID,P[i].ID_child_number,P[i].Mass,rest_mass,r2_0);
		fflush(fsink1);}
//        else if(ibh==1 && r2_1<2*grP[1].Mass){
        else if(ibh==1 && r2_1<EHradius(grP[1].Mass,grP[1].Spin)){
                fprintf(fsink2,"%1.10g\t %d %d %g %g %g\n",All.Time,P[i].ID,P[i].ID_child_number,P[i].Mass,rest_mass,r2_1);
		fflush(fsink2);}
#ifdef BOX_PERIODIC
	if(fabs(P[i].Pos[2]-CellP[i].VelPred[2]*All.TimeStep-boxHalf)<(PRINT_OUTFLOW_Z) && fabs(P[i].Pos[2]-boxHalf)>PRINT_OUTFLOW_Z){
#else
	if(fabs(P[i].Pos[2]-CellP[i].VelPred[2]*All.TimeStep)<(PRINT_OUTFLOW_Z) && fabs(P[i].Pos[2])>PRINT_OUTFLOW_Z){
#endif
		fprintf(foutflow,"%1.10g\t %d %d %g %g %g %g %g\n",All.Time,P[i].ID,P[i].ID_child_number,P[i].Mass,rest_mass,P[i].Pos[0],P[i].Pos[1],P[i].Pos[2]);
		fflush(foutflow);}
        }
  }
  fclose(fsink1);
  fclose(fsink2);
  fclose(foutflow);
#ifndef GR_USECBWAVES
  free(BHprops_all);
#endif
}
#endif

#if defined(PARTICLE_EXCISION) && defined(USE_KERRSCHILD_COORDINATES)
void smooth_conservatives_inside_horizon(int i, double r, double r_EH, double M)
{
	double Chi, exponent;
	double v_local[3], v_localPred[3];
	for(int k=0; k<3; k++) v_local[k] = P[i].Vel[k];
	for(int k=0; k<3; k++) v_localPred[k] = CellP[i].VelPred[k];
	Vphys(v_local,CellP[i].gr.beta,CellP[i].gr.alpha);
	Vphys(v_localPred,CellP[i].gr.beta,CellP[i].gr.alpha);

	exponent = (0.5*sqrt(GR_EXCISION_RADIUS*M)*(r-r_EH)*(r-r_EH))/(r_EH - 0.4);
	Chi = exp(-exponent);

	v_local[0] *= Chi;
	v_local[1] *= Chi;
	v_local[2] *= Chi;
	v_localPred[0] *= Chi;
	v_localPred[1] *= Chi;
	v_localPred[2] *= Chi;

	Vcoord(v_local,CellP[i].gr.beta,CellP[i].gr.alpha);
	Vcoord(v_localPred,CellP[i].gr.beta,CellP[i].gr.alpha);

	P[i].Vel[0] = v_local[0];
	P[i].Vel[1] = v_local[1];
	P[i].Vel[2] = v_local[2];
	CellP[i].VelPred[0] = v_localPred[0];
	CellP[i].VelPred[1] = v_localPred[1];
	CellP[i].VelPred[2] = v_localPred[2];

#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION>=5)
	double par_v_local[3];
	for(int k=0; k<3; k++) par_v_local[k] = CellP[i].ParticleVel[k];
	Vphys(par_v_local,CellP[i].gr.beta,CellP[i].gr.alpha);
	par_v_local[0] *= Chi;
	par_v_local[1] *= Chi;
	par_v_local[2] *= Chi;
	Vcoord(par_v_local,CellP[i].gr.beta,CellP[i].gr.alpha);
	CellP[i].ParticleVel[0] = par_v_local[0];
	CellP[i].ParticleVel[1] = par_v_local[1];
	CellP[i].ParticleVel[2] = par_v_local[2];
#endif
#ifdef MAGNETIC
	CellP[i].B[0] *= Chi;
	CellP[i].B[1] *= Chi;
	CellP[i].B[2] *= Chi;
	CellP[i].BPred[0] *= Chi;
	CellP[i].BPred[1] *= Chi;
	CellP[i].BPred[2] *= Chi;
#endif
}
#endif



#ifdef PARTICLE_EXCISION
//This sets the mass to 0 if outside the physical boundaries (e.g. the horizon radius)
void check_boundary(int target)
{
   double horizon;
   double r_min;
#if (GR_METRIC==1)
   r_min=sqrt(P[target].Pos[0]*P[target].Pos[0]+P[target].Pos[1]*P[target].Pos[1]+P[target].Pos[2]*P[target].Pos[2]);
   double a[3] = {0,0,All.MetricParams[1]};
   horizon=EHradius(All.MetricParams[0],a);
#elif (GR_METRIC==2)
   double r[2],dp[3];
   int i;
   for(i=0;i<2;i++)
   {
     dp[0]=P[target].Pos[0]-grP[i].Pos[0];
     dp[1]=P[target].Pos[1]-grP[i].Pos[1];
     dp[2]=P[target].Pos[2]-grP[i].Pos[2];
     r[i]=dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2];
   }
   horizon = r[0]<r[1] ? EHradius(grP[0].Mass,grP[0].Spin) : EHradius(grP[1].Mass,grP[1].Spin);
   r_min = sqrt(fmin(r[0],r[1]));
#else
   return;
#endif
   if(r_min<GR_EXCISION_RADIUS*horizon){P[target].Mass=0;} // printf("Excised particle\n");}
}
#endif

void gr_get_conservative_variables(struct gr_primitives *prims,struct gr_conservatives *cons,double Wlorentz,double gamma[],double gdet)
{
	double hstar=1.0 + prims->u + prims->p/prims->rho0;
	double Pstar=prims->p;
        double v_cov[3];
        //Get covariant velocity (by contracting with the covariant metric)
        raise_lower_indices(gamma,prims->v,v_cov);
        cons->D=gdet*prims->rho0*Wlorentz;
	cons->S[0] = cons->S[1] = cons->S[2] = 0;
	cons->tau=0;
#ifdef MAGNETIC
	cons->B[0]=gdet*prims->B[0];
	cons->B[1]=gdet*prims->B[1];
	cons->B[2]=gdet*prims->B[2];

        //Now modify the pressure to compute the momentum-energy tensor
        Pstar+=prims->bsq/2;
        hstar+=prims->bsq/prims->rho0;

	cons->S[0] = - gdet*prims->b[0]*prims->b[1];
	cons->S[1] = - gdet*prims->b[0]*prims->b[2];
	cons->S[2] = - gdet*prims->b[0]*prims->b[3];
        cons->tau  = - gdet*prims->b[0]*prims->b[0];
#endif
        double DhW=cons->D * Wlorentz * hstar;
        cons->S[0] += DhW*v_cov[0];
        cons->S[1] += DhW*v_cov[1];
        cons->S[2] += DhW*v_cov[2];
        cons->tau  += DhW - gdet*Pstar - cons->D;

        if(cons->tau<0){printf("WARNING: tau<0 in Conservative evaluation. Tau %g \n", cons->tau); cons->tau = MIN_REAL_NUMBER;}
        cons->p = Pstar;
}


#ifdef GR_BETACOOLING
//Beta cooling like in Noble+(2012). Using target sound speed (Farris+(2014)) instead of target entropy
void implicit_cooling(void)
{
    PRINT_STATUS("...cooling for GR quantities");
    /* Determine indices of active particles. */
    int N_active=0, i, j, target, *active_indices; active_indices = (int *) malloc(N_gas * sizeof(int));
    for (i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type != 0) {continue;}
        if(P[i].Mass <= 0) {continue;}
        active_indices[N_active] = i;
        N_active++;
    }

#ifdef _OPENMP
#pragma omp parallel private(i, j)
#endif
    { /* open parallel block */
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for(j=0;j<N_active;j++)
    {
        target=active_indices[j]; /* actual particle index */
        /* do the actual cooling */
	int k;
	double res_con2prim=0;
	struct gr_primitives Prims;
	struct gr_conservatives Cons;
	double Wlorentz;
	double v_local[3];
	double v_cov[3];
	int iter=0;
        int gr_limiter=0;
	double Vi=P[target].Mass/CellP[target].Density;
	double dtime = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(target);
	//Load the primitives
	for(k=0;k<3;k++) v_local[k]=CellP[target].VelPred[k];
	Vphys(v_local,CellP[target].gr.beta,CellP[target].gr.alpha);
        double dpl[3]; for(k=0;k<3;k++) {dpl[k]=P[target].Pos[k];}
#if (GR_METRIC==2)
        for(k=0;k<3;k++) {dpl[k] = -P[target].Min_xyz_to_Sink[k];}
#endif
        if(dpl[0]*dpl[0]+dpl[1]*dpl[1]+dpl[2]*dpl[2] <= pow(GAMMA_LIMITER_RADIUS,2)) gr_limiter=1;
//      if(P[target].Pos[0]<LEFT_GAMMA_BOUNDARY||P[target].Pos[0]>RIGHT_GAMMA_BOUNDARY) gr_limiter=1;
        Wlorentz=lorentz_factor(CellP[target].gr.gamma,v_local,gr_limiter);
        Prims.rho0=P[target].Mass/Wlorentz/CellP[target].gr.gdet/Vi;
        Prims.u=CellP[target].InternalEnergyPred;
        Prims.p=gr_pressure_rho0_u(Prims.rho0,Prims.u);
        for(k=0;k<3;k++) Prims.v[k]=v_local[k];
	raise_lower_indices(CellP[target].gr.gamma,Prims.v,v_cov);
#ifdef MAGNETIC
        for(k=0;k<3;k++) Prims.B[k]=CellP[target].BPred[k]/Vi/CellP[target].gr.gdet;
        Prims.bsq=Get_Gas_BField_RestFrame(target,Prims.B, Prims.v,Prims.b);
#endif
        //Determine conservative variables
        gr_get_conservative_variables(&Prims,&Cons,Wlorentz,CellP[target].gr.gamma,CellP[target].gr.gdet);

	//Calculate cooling parameters
        double Omega_cool, target_SoundSpeed, eff_distance;
        double dp[3]; for(k=0;k<3;k++) {dp[k]=P[target].Pos[k];}
#ifdef BOX_PERIODIC
        dp[0] -= boxHalf_X;
        dp[1] -= boxHalf_Y;
        dp[2] -= boxHalf_Z;
#endif
	double r = sqrt(dp[0]*dp[0]+dp[1]*dp[1]); //Cylindrical radius
#if (GR_METRIC==2)
        double r1,r2,r1_EH,r2_EH,r_m;
        r1 = sqrt((P[target].Pos[0]-grP[0].Pos[0])*(P[target].Pos[0]-grP[0].Pos[0])+(P[target].Pos[1]-grP[0].Pos[1])*(P[target].Pos[1]-grP[0].Pos[1]));
        r2 = sqrt((P[target].Pos[0]-grP[1].Pos[0])*(P[target].Pos[0]-grP[1].Pos[0])+(P[target].Pos[1]-grP[1].Pos[1])*(P[target].Pos[1]-grP[1].Pos[1]));
	r_m = 0.5*(r1+r2);
	eff_distance = (1-exp(-r_m))*(r1*exp(-r1)+r2*exp(-r2))/(exp(-r1)+exp(-r2)) + exp(-r_m)*r_m;
        r1_EH = 2*grP[0].Mass;
        r2_EH = 2*grP[1].Mass;
        Omega_cool = sqrt(grP[0].Mass/(r1+r1_EH)/(r1+r1_EH)/(r1+r1_EH) + grP[1].Mass/(r2+r2_EH)/(r2+r2_EH)/(r2+r2_EH));
        target_SoundSpeed = HR_COOLING * sqrt(grP[0].Mass/(r1+r1_EH) + grP[1].Mass/(r2+r2_EH));
#elif (GR_METRIC==1)
        double r_EH = 2*All.MetricParams[0];
        Omega_cool = sqrt(All.MetricParams[0]/(r+r_EH)/(r+r_EH)/(r+r_EH));
        target_SoundSpeed = HR_COOLING * sqrt(All.MetricParams[0]/(r+r_EH));
	eff_distance = r;
#endif
        double cooling_factor = (Get_Gas_effective_soundspeed_i(target)-target_SoundSpeed)/DMAX(target_SoundSpeed,Get_Gas_effective_soundspeed_i(target));
	if(fabs(P[target].Pos[2])>HR_COOLING*eff_distance){cooling_factor *= exp(1 - fabs(P[target].Pos[2])/HR_COOLING/eff_distance);}
	double cool_rate_S[3];
	double cool_rate_tau;
	//First predictor step
        cool_rate_tau = -(Prims.rho0*Prims.u/2/M_PI/GR_BETA*Omega_cool*cooling_factor*Wlorentz)*CellP[target].gr.alpha*CellP[target].gr.gdet;
        for(k=0;k<3;k++) {cool_rate_S[k] = -(Prims.rho0*Prims.u/2/M_PI/GR_BETA*Omega_cool*cooling_factor*Wlorentz*v_cov[k])*CellP[target].gr.alpha*CellP[target].gr.gdet;}

//	if(P[target].ID==PARTICLE || P[target].ID==PARTICLE2 || P[target].ID==PARTICLE3){printf("Inizio Prims.u %g - Prims.v %g %g %g \nCooling rate tau %g - S %g %g %g \n",Prims.u,Prims.v[0],Prims.v[1],Prims.v[2],cool_rate_tau,cool_rate_S[0],cool_rate_S[1],cool_rate_S[2]);}
//	if(P[target].ID==PARTICLE || P[target].ID==PARTICLE2 || P[target].ID==PARTICLE3){printf("BH1 Pos %g %g %g - BH2 Pos %g %g %g\nPosizione %g %g %g\ncs %g - target %g \n",grP[0].Pos[0],grP[0].Pos[1],grP[0].Pos[2],grP[1].Pos[0],grP[1].Pos[1],grP[1].Pos[2],P[target].Pos[0],P[target].Pos[1],P[target].Pos[2],Get_Gas_effective_soundspeed_i(target),target_SoundSpeed);}

	Cons.tau += cool_rate_tau*dtime;
	for(k=0;k<3;k++) Cons.S[k] += cool_rate_S[k]*dtime;
        if(Cons.tau <= 0) {printf("WARNING: Tau < 0 after the cooling procedure (prediction)\n"); Cons.tau -= cool_rate_tau*dtime; for(k=0;k<3;k++){Cons.S[k] -= cool_rate_S[k]*dtime;} continue;}

	//Recover the primitives
	res_con2prim=gr_get_primitive_variables(&Cons,&Prims,&(CellP[target].gr), gr_limiter);
	if(res_con2prim<0) {printf("Warning: failed inversion with code %g - particle %d at %g %g %g - kernelsize %g\n",res_con2prim, P[target].ID,P[target].Pos[0],P[target].Pos[1],P[target].Pos[2],P[target].KernelRadius);endrun(9993010);}

//	if(P[target].ID==PARTICLE || P[target].ID==PARTICLE2 || P[target].ID==PARTICLE3){printf("iter 0 Prims.u %g - Prims.v %g %g %g \n",Prims.u,Prims.v[0],Prims.v[1],Prims.v[2]);}

	double cool_rate_S_corr[3];
	for(k=0;k<3;k++) cool_rate_S_corr[k] = cool_rate_S[k];
	double cool_rate_tau_corr=cool_rate_tau;
	double u_error=MAX_REAL_NUMBER;
	do{
		u_error = Prims.u;
		cool_rate_tau = cool_rate_tau_corr;
		for(k=0;k<3;k++) cool_rate_S[k] = cool_rate_S_corr[k];
		//Corrected conservatives
		Wlorentz=lorentz_factor(CellP[target].gr.gamma,Prims.v,gr_limiter);

		cool_rate_tau_corr = -(Prims.rho0*Prims.u/2/M_PI/GR_BETA*Omega_cool*cooling_factor*Wlorentz)*CellP[target].gr.alpha*CellP[target].gr.gdet;
		for(k=0;k<3;k++) {cool_rate_S_corr[k] = -(Prims.rho0*Prims.u/2/M_PI/GR_BETA*Omega_cool*cooling_factor*Wlorentz*v_cov[k])*CellP[target].gr.alpha*CellP[target].gr.gdet;}

		//Correct prediction
		Cons.tau += (cool_rate_tau_corr - cool_rate_tau)*dtime;
		for(k=0;k<3;k++) Cons.S[k] += (cool_rate_S_corr[k] - cool_rate_S[k])*dtime;
                if(Cons.tau <= 0) {printf("WARNING: Tau < 0 after the cooling procedure, iteration %d\n",iter); Cons.tau = MIN_REAL_NUMBER; break;}

		//Recover the primitives
		res_con2prim=gr_get_primitive_variables(&Cons,&Prims,&(CellP[target].gr), gr_limiter);
		if(res_con2prim<0) {printf("Warning: failed inversion with code %g - particle %d at %g %g %g - kernelsize %g\n",res_con2prim, P[target].ID,P[target].Pos[0],P[target].Pos[1],P[target].Pos[2],P[target].KernelRadius);endrun(9993010);}

//		if(P[target].ID==PARTICLE || P[target].ID==PARTICLE2 || P[target].ID==PARTICLE3){printf("iter %d Prims.u %g - Prims.v %g %g %g \nCooling rate tau %g - S %g %g %g \n\n",iter, Prims.u,Prims.v[0],Prims.v[1],Prims.v[2],cool_rate_tau,cool_rate_S[0],cool_rate_S[1],cool_rate_S[2]);}


		//Compute error
		u_error = Prims.u - u_error;
		iter++;
		if(iter>IMPLICIT_MAX_NITER)
		{
			printf("WARNING: Iter >100 in cooling procedure - u %g - error %g\n", Prims.u, u_error);
			break;
		}
	}while(u_error>IMPLICIT_SOLVER);

        //Convert vel back to dx/dt
        Vcoord(Prims.v,CellP[target].gr.beta,CellP[target].gr.alpha);
        for(k=0;k<3;k++) CellP[target].VelPred[k] += (Prims.v[k]-CellP[target].VelPred[k]);
        CellP[target].InternalEnergyPred += (Prims.u-CellP[target].InternalEnergyPred);
	//Keeping kicked and predicted quanitites synchronized
	CellP[target].InternalEnergy = CellP[target].InternalEnergyPred;
	for(k=0;k<3;k++) P[target].Vel[k] = CellP[target].VelPred[k];
	CellP[target].Pressure = get_pressure(target);
    }
    } /* close parallel block */
    free(active_indices); /* free memory */
}
#endif

#if (GR_METRIC==2) && defined(GR_REMOVE_GRAD)
void remove_radial_gradient_at_outer_edge(double *grad, int target)
{
	int k;
	double r;
	double tmp_vec[2]; //We are dealing with the radial gradient in cylindrical coordinates
	double dp[3]; for(k=0;k<2;k++) {dp[k]=P[target].Pos[k]; tmp_vec[k] = grad[k];}
#ifdef BOX_PERIODIC
        dp[0] -= boxHalf_X;
        dp[1] -= boxHalf_Y;
#endif
	r = sqrt(dp[0]*dp[0]+dp[1]*dp[1]);
	if(r < GR_REMOVE_GRAD*(1-8/100))
	{
		return;
	}else{
		if(dp[0]*tmp_vec[0]/r + dp[1]*tmp_vec[1]/r > 0){
			grad[0] -= dp[0]/r * (dp[0]*tmp_vec[0]/r + dp[1]*tmp_vec[1]/r);
			grad[1] -= dp[1]/r * (dp[0]*tmp_vec[0]/r + dp[1]*tmp_vec[1]/r);
		}
	}
}
#endif

long double vector_sqnorm(double gamma[],double v_in[])
{
        double v_rl[3];
        long double vnorm=0;
        int i;
	raise_lower_indices(gamma,v_in,v_rl); //g_ij * v^i
        for(i=0;i<3;i++) {vnorm+= v_rl[i]*v_in[i];}
        return vnorm;
}

double lorentz_factor(double gamma[],double v_cont[], int gr_limiter)
{
	//To avoid too coarse approximation of the velocity near the maximum gamma value, we use quadruple precision here
	long double denom;
        long double R=vector_sqnorm(gamma,v_cont);
	if(gr_limiter==1)
	{
        	long double RMAX=1.0-1./GR_GAMMA_LIMITER/GR_GAMMA_LIMITER;
        	long double Rstar=0.5*(RMAX+R-fabsl(R-RMAX));

        	v_cont[0]*=(double)sqrtl(Rstar/(R+1e-70));
        	v_cont[1]*=(double)sqrtl(Rstar/(R+1e-70));
        	v_cont[2]*=(double)sqrtl(Rstar/(R+1e-70));
		denom = 1. - Rstar;
	}else{
	denom = 1. - R;
	if(denom<MIN_REAL_NUMBER) denom=MIN_REAL_NUMBER;
	}
        return (double)(1./sqrtl(denom));

}

void raise_lower_indices(double g[],double vin[], double vout[])
{
	vout[0]=g[0]*vin[0]+g[1]*vin[1]+g[2]*vin[2];
	vout[1]=g[1]*vin[0]+g[3]*vin[1]+g[4]*vin[2];
	vout[2]=g[2]*vin[0]+g[4]*vin[1]+g[5]*vin[2];
}

void rotate_tensor(double A[],double cosphi,double sinphi,double costh,double sinth,double Aprime[])
{
	long double R[3][3]={{costh*cosphi,-costh*sinphi,-sinth},{sinphi,cosphi,0},{sinth*cosphi,-sinth*sinphi,costh}};
	long double tmp[3][3]={{A[0],A[1],A[2]},{A[1],A[3],A[4]},{A[2],A[4],A[5]}};
	long double rtmp[3][3];
	int i,j;

	//Now multiply on the right
	for(i=0;i<3;i++)
	 for(j=0;j<3;j++)
	 {
	   rtmp[i][j]=tmp[i][0]*R[j][0]+tmp[i][1]*R[j][1]+tmp[i][2]*R[j][2];
	 }
	//Now multiply on the left
	for(i=0;i<3;i++)
	 for(j=0;j<3;j++)
	 {
	   tmp[i][j]=R[i][0]*rtmp[0][j]+R[i][1]*rtmp[1][j]+R[i][2]*rtmp[2][j];
	 }
	Aprime[0]=tmp[0][0];
	Aprime[1]=tmp[0][1];
	Aprime[2]=tmp[0][2];
	Aprime[3]=tmp[1][1];
	Aprime[4]=tmp[1][2];
	Aprime[5]=tmp[2][2];
}

void unpack_4tensor(double Gin[],double Gout[][4])
{
	Gout[0][0]=Gin[0]; //tt
	Gout[1][1]=Gin[4]; //xx
	Gout[2][2]=Gin[7]; //yy
	Gout[3][3]=Gin[9]; //zz

	Gout[0][1]=Gout[1][0]=Gin[1]; //tx
	Gout[0][2]=Gout[2][0]=Gin[2]; //ty
	Gout[0][3]=Gout[3][0]=Gin[3]; //tz

	Gout[1][2]=Gout[2][1]=Gin[5]; //xy
	Gout[1][3]=Gout[3][1]=Gin[6]; //xz
	Gout[2][3]=Gout[3][2]=Gin[8]; //yz
}

void unpack_3tensor(double Gin[],double Gout[][3])
{
	Gout[0][0]=Gin[0]; //xx
	Gout[1][1]=Gin[3]; //yy
	Gout[2][2]=Gin[5]; //z

	Gout[0][1]=Gout[1][0]=Gin[1]; //xy
	Gout[0][2]=Gout[2][0]=Gin[2]; //xz
	Gout[1][2]=Gout[2][1]=Gin[4]; //yz
}

#endif
