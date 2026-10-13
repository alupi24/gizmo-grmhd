#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "gr_source.h"

#ifdef GENERAL_RELATIVITY

void compute_initial_accelerations(void)
{
	int i;
	double dt_0=DMIN(All.MinSizeTimestep*1.1,All.MaxSizeTimestep);
        for(i=0;i<NumPart;i++)
	{
		if(P[i].Type!=0) continue;
		do_kick_general_relativity(i,dt_0,P[i].Mass,P[i].Mass,0);
	}
}

void do_kick_general_relativity(int target,double dt_kick, double mass_old,double mass_new,int mode)
{
	int k;
	double res_con2prim=0;

	struct gr_primitives Prims;
	struct gr_conservatives Cons,Src;
	double Wlorentz;
	//Get particle volume
	double Vi, Vnew;
	if(mode==1) {Vnew = mass_new/CellP[target].Density; Vi = Vnew * exp(-CellP[target].DtVolume*dt_kick/Vnew);}
	if(mode==0) {Vi = mass_old/CellP[target].Density; Vnew = Vi * exp(CellP[target].DtVolume*dt_kick/Vi);}
	//Convert vel into physical velocity
        double v_local[3];
        for(k=0;k<3;k++) v_local[k]=P[target].Vel[k];
	Vphys(v_local,CellP[target].gr.beta,CellP[target].gr.alpha);

#ifdef GR_TEST_PARTICLE
	mass_old=mass_new=Vi=Vnew=1;
#endif
	//Now fill the primitive variable structure
	int gr_limiter=0;
	double dpl[3]; for(k=0;k<3;k++) {dpl[k]=P[target].Pos[k];}
#if (GR_METRIC==2)
	for(k=0;k<3;k++) {dpl[k] = -P[target].Min_xyz_to_Sink[k];}
#endif
	if(dpl[0]*dpl[0]+dpl[1]*dpl[1]+dpl[2]*dpl[2] <= pow(GAMMA_LIMITER_RADIUS,2)) gr_limiter=1;
//	if(P[target].Pos[0]<LEFT_GAMMA_BOUNDARY||P[target].Pos[0]>RIGHT_GAMMA_BOUNDARY) gr_limiter=1;
	Wlorentz=lorentz_factor(CellP[target].gr.gamma,v_local, gr_limiter);
	Prims.rho0=mass_old/Wlorentz/CellP[target].gr.gdet/Vi;
	Prims.u=CellP[target].InternalEnergy;
	Prims.p=gr_pressure_rho0_u(Prims.rho0,Prims.u);
	for(k=0;k<3;k++) Prims.v[k]=v_local[k];
#ifdef MAGNETIC
        for(k=0;k<3;k++) Prims.B[k]=CellP[target].B[k]/Vi/CellP[target].gr.gdet;
        Prims.bsq=Get_Gas_BField_RestFrame(target,Prims.B, Prims.v,Prims.b);
#endif
	//Determine conservative variables
	gr_get_conservative_variables(&Prims,&Cons,Wlorentz,CellP[target].gr.gamma,CellP[target].gr.gdet);
	//Update quantities accounting for fluxes and source terms

	//Add source terms
#if !(defined(HYDRO_MESHLESS_FINITE_MASS) && (GENERAL_RELATIVITY==1))
	if(mode==0)
#endif
	{gr_add_source_term(target,&Prims,&Cons,dt_kick,1-mode);}
	if(P[target].Mass==0){return;}
	Cons.D = mass_new/Vnew;
	Cons.tau = (Cons.tau*Vi + (CellP[target].DtEnergy*dt_kick))/Vnew;
	for(k=0;k<3;k++) Cons.S[k] = (Cons.S[k]*Vi + (CellP[target].DtMomentum[k]*dt_kick))/Vnew;
#ifdef MAGNETIC
	for(k=0;k<3;k++) Cons.B[k] = (Cons.B[k]*Vi + (CellP[target].DtB[k]*dt_kick))/Vnew;
#endif

	//Find predicted primitives
	res_con2prim=gr_get_primitive_variables(&Cons,&Prims,&(CellP[target].gr), gr_limiter);
	if(res_con2prim<0) {printf("Warning: failed inversion with code %g - particle %d at %g %g %g - kernelsize %g\n",res_con2prim, P[target].ID,P[target].Pos[0],P[target].Pos[1],P[target].Pos[2],P[target].KernelRadius); endrun(9992405);}

#if !(defined(HYDRO_MESHLESS_FINITE_MASS) && (GENERAL_RELATIVITY==1))
	//Now add source terms
	if(mode==1){gr_add_source_term(target,&Prims,&Cons,dt_kick,1-mode);if(P[target].Mass==0){return;}}
#endif


#if(GR_METRIC==2)
	if(Prims.u>=1.0){
		if(sqrt(dpl[0]*dpl[0]+dpl[1]*dpl[1]+dpl[2]*dpl[2])>5){
			printf("WARNING: u>1.0 for particle %d at %g from BH. u = %g\n", P[target].ID,sqrt(dpl[0]*dpl[0]+dpl[1]*dpl[1]+dpl[2]*dpl[2]),Prims.u);
		}
		P[target].Mass=0;
	}
#endif

	//Finally determine the update rate using old and new primitives
	double rate_fac=1./dt_kick;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME //Account for the mass change
        rate_fac*=CellP[target].MassTrue/P[target].Mass;
#endif

	//Convert vel back to dx/dt
	Vcoord(Prims.v,CellP[target].gr.beta,CellP[target].gr.alpha);

	//Finally update accelerations/udot
	CellP[target].DtInternalEnergy = (Prims.u-CellP[target].InternalEnergy)*rate_fac;
#ifdef GR_ENTROPY_SWITCH
	//Entropy Switch
	if(CellP[target].Entropy != 0){
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
		CellP[target].Entropy = CellP[target].Entropy + (CellP[target].DtEntropy*dt_kick)/Vi;
#endif
		if(0.5*Prims.bsq >= 0.1*Prims.p){
			double EnergyDensFromEntropy;
			EnergyDensFromEntropy = CellP[target].Entropy * pow(Prims.rho0,GAMMA(target)-1) / (GAMMA(target)-1);
			CellP[target].DtInternalEnergy = (EnergyDensFromEntropy - CellP[target].InternalEnergy)*rate_fac; ////Entropy Switch
		}
	}
#endif

	for(k=0;k<3;k++) CellP[target].HydroAccel[k] = (Prims.v[k]-P[target].Vel[k])*rate_fac;
}


void do_drift_general_relativity(int target,double dt_drift)
{
	int k;
	struct gr_primitives Prims;
	struct gr_conservatives Cons,ConsPred;
        double x_err[3];
	double res_con2prim,Wlorentz;
	double v_old[3]={P[target].Vel[0],P[target].Vel[1],P[target].Vel[2]};
	double v_new[3];
	//Get particle volume
	double Vi = P[target].Mass/CellP[target].Density;
	int gr_limiter=0;
	double dpl[3]; for(k=0;k<3;k++) {dpl[k]=P[target].Pos[k];}
#if (GR_METRIC==2)
	for(k=0;k<3;k++) {dpl[k] = -P[target].Min_xyz_to_Sink[k];}
#endif
	if(dpl[0]*dpl[0]+dpl[1]*dpl[1]+dpl[2]*dpl[2] <= pow(GAMMA_LIMITER_RADIUS,2)) gr_limiter=1;
//	if(P[target].Pos[0]<LEFT_GAMMA_BOUNDARY||P[target].Pos[0]>RIGHT_GAMMA_BOUNDARY) gr_limiter=1;


#ifdef GR_TEST_PARTICLE
	P[target].Mass=Vi=1;
	double r=sqrt(P[target].Pos[0]*P[target].Pos[0]+P[target].Pos[1]*P[target].Pos[1]+P[target].Pos[2]*P[target].Pos[2]);
	if(r<=1.001) P[target].Mass=0;
#endif
#if defined(PARTICLE_EXCISION) && (GR_METRIC==2)
	double mass_tmp = P[target].Mass;
	check_boundary(target);
	if(P[target].Mass==0){
		Vphys(v_old,CellP[target].gr.beta,CellP[target].gr.alpha);
		Wlorentz=lorentz_factor(CellP[target].gr.gamma,v_old, gr_limiter);
 	 	char path[300];
	  	sprintf(path,"%s/%s",All.OutputDir,"mdot_s1.dat");
		FILE* fsink1=fopen(path,"a");
	 	sprintf(path,"%s/%s",All.OutputDir,"mdot_s2.dat");
		FILE* fsink2=fopen(path,"a");
		double r2_0=0,r2_1=0;
		int ibh;
	        r2_0=sqrt((P[target].Pos[0]-grP[0].Pos[0])*(P[target].Pos[0]-grP[0].Pos[0])+(P[target].Pos[1]-grP[0].Pos[1])*(P[target].Pos[1]-grP[0].Pos[1])+
			(P[target].Pos[2]-grP[0].Pos[2])*(P[target].Pos[2]-grP[0].Pos[2]));
        	r2_1=sqrt((P[target].Pos[0]-grP[1].Pos[0])*(P[target].Pos[0]-grP[1].Pos[0])+(P[target].Pos[1]-grP[1].Pos[1])*(P[target].Pos[1]-grP[1].Pos[1])+
			(P[target].Pos[2]-grP[1].Pos[2])*(P[target].Pos[2]-grP[1].Pos[2]));
	        if(r2_0/grP[0].Mass<r2_1/grP[1].Mass) ibh=0; else ibh=1;
	        double rest_mass=mass_tmp/Wlorentz/CellP[target].gr.gdet;
	        if(ibh==0){
	                fprintf(fsink1,"%1.10g\t %d %d %g %g %g\n",All.Time,P[target].ID,P[target].ID_child_number,mass_tmp,rest_mass,r2_0);
        	        fflush(fsink1);}
        	else if(ibh==1){
                	fprintf(fsink2,"%1.10g\t %d %d %g %g %g\n",All.Time,P[target].ID,P[target].ID_child_number,mass_tmp,rest_mass,r2_1);
                	fflush(fsink2);}
		fclose(fsink1);
		fclose(fsink2);
       	 	return;
	}
#endif
	/* --- Determine conservative variables for predicted quantities (to be updated at the end of the implicit loop) --- */

	//Fill the primitive variable structure first
        for(k=0;k<3;k++) Prims.v[k]=CellP[target].VelPred[k];
        Vphys(Prims.v,CellP[target].gr.beta,CellP[target].gr.alpha);
	Wlorentz=lorentz_factor(CellP[target].gr.gamma,Prims.v, gr_limiter);
	Prims.rho0 = P[target].Mass/Wlorentz/CellP[target].gr.gdet/Vi;
        Prims.u    = CellP[target].InternalEnergyPred;
        Prims.p    = gr_pressure_rho0_u(Prims.rho0,Prims.u);
#ifdef MAGNETIC
        for(k=0;k<3;k++){ Prims.B[k]=CellP[target].BPred[k]/Vi/CellP[target].gr.gdet;}
        Prims.bsq=Get_Gas_BField_RestFrame(target,Prims.B,Prims.v,Prims.b);
#endif

	//Determine conservative variables
        gr_get_conservative_variables(&Prims,&ConsPred,Wlorentz,CellP[target].gr.gamma,CellP[target].gr.gdet);

	/* --- Determine conservative variables for kicked quantities (to be updated at the end of the implicit loop) --- */

	//Fill the primitive variable structure first
	for(k=0;k<3;k++) Prims.v[k] = P[target].Vel[k];
	Vphys(Prims.v,CellP[target].gr.beta,CellP[target].gr.alpha);
	Wlorentz=lorentz_factor(CellP[target].gr.gamma,Prims.v, gr_limiter);
	Prims.rho0 = P[target].Mass/Wlorentz/CellP[target].gr.gdet/Vi;
	Prims.u    = CellP[target].InternalEnergy;
	Prims.p    = gr_pressure_rho0_u(Prims.rho0,Prims.u);
#ifdef MAGNETIC
        for(k=0;k<3;k++){ Prims.B[k]=CellP[target].B[k]/Vi/CellP[target].gr.gdet;}
        Prims.bsq=Get_Gas_BField_RestFrame(target,Prims.B,Prims.v,Prims.b);
#endif

	//Determine conserved variables
	gr_get_conservative_variables(&Prims,&Cons,Wlorentz,CellP[target].gr.gamma,CellP[target].gr.gdet);

	/* --- Start implicit drift loop, accounting for metric change --- */
#ifdef IMPLICIT_SOLVER
	double error=MAX_REAL_NUMBER;
	int iter=0;
	do{
#endif
	 for(k=0;k<3;k++) x_err[k]=P[target].Pos[k];

	 //Update metric at predicted location
	 set_metric_coefficients(P[target].Pos,P[target].Potential,&(CellP[target].gr));

      	 res_con2prim=gr_get_primitive_variables(&Cons,&Prims,&(CellP[target].gr), gr_limiter);
         if(res_con2prim<-1) {printf("Warning: failed inversion with code in predict %g - ID %d - pos %g %g %g\n", res_con2prim,P[target].ID, P[target].Pos[0],P[target].Pos[1],P[target].Pos[2]);
#if (GR_METRIC==2)
		double r2_0=sqrt((P[target].Pos[0]-grP[0].Pos[0])*(P[target].Pos[0]-grP[0].Pos[0])+(P[target].Pos[1]-grP[0].Pos[1])*(P[target].Pos[1]-grP[0].Pos[1])+
                        (P[target].Pos[2]-grP[0].Pos[2])*(P[target].Pos[2]-grP[0].Pos[2]));
		double r2_1=sqrt((P[target].Pos[0]-grP[1].Pos[0])*(P[target].Pos[0]-grP[1].Pos[0])+(P[target].Pos[1]-grP[1].Pos[1])*(P[target].Pos[1]-grP[1].Pos[1])+
                        (P[target].Pos[2]-grP[1].Pos[2])*(P[target].Pos[2]-grP[1].Pos[2]));
	 	if(r2_0<2*EHradius(grP[0].Mass,grP[0].Spin) || r2_1<2*EHradius(grP[1].Mass,grP[1].Spin)){P[target].Mass=0; return;}
		else{endrun(9992406);}
#else
		endrun(9992406);
#endif
	 }

	 //Copy to vnew for conversion to dx/dt
	 for(k=0;k<3;k++) v_new[k]=Prims.v[k];
	 //Convert vel back to dx/dt
	 Vcoord(v_new,CellP[target].gr.beta,CellP[target].gr.alpha);
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION<5)
	 //No need to correct position except for (quasi-)lagrangian approaches
	 break;
#endif

	 //Finally update position using corrected values
	 for(k=0;k<3;k++) {P[target].Pos[k] += 0.5*(v_new[k]-v_old[k])*dt_drift; v_old[k]=v_new[k];}

#if defined(PARTICLE_EXCISION) && (GR_METRIC==2)
	 mass_tmp = P[target].Mass;
	 check_boundary(target);
	 if(P[target].Mass==0){
 	 	 char path[300];
	  	 sprintf(path,"%s/%s",All.OutputDir,"mdot_s1.dat");
		 FILE* fsink1=fopen(path,"a");
	 	 sprintf(path,"%s/%s",All.OutputDir,"mdot_s2.dat");
		 FILE* fsink2=fopen(path,"a");
		 double r2_0=0,r2_1=0;
		 int ibh;
	         r2_0=sqrt((P[target].Pos[0]-grP[0].Pos[0])*(P[target].Pos[0]-grP[0].Pos[0])+(P[target].Pos[1]-grP[0].Pos[1])*(P[target].Pos[1]-grP[0].Pos[1])+(P[target].Pos[2]-grP[0].Pos[2])*(P[target].Pos[2]-grP[0].Pos[2]));
        	 r2_1=sqrt((P[target].Pos[0]-grP[1].Pos[0])*(P[target].Pos[0]-grP[1].Pos[0])+(P[target].Pos[1]-grP[1].Pos[1])*(P[target].Pos[1]-grP[1].Pos[1])+(P[target].Pos[2]-grP[1].Pos[2])*(P[target].Pos[2]-grP[1].Pos[2]));
	         if(r2_0/grP[0].Mass<r2_1/grP[1].Mass) ibh=0; else ibh=1;
	         double rest_mass=mass_tmp/Wlorentz/CellP[target].gr.gdet;
	         if(ibh==0){
	                 fprintf(fsink1,"%1.10g\t %d %d %g %g %g\n",All.Time,P[target].ID,P[target].ID_child_number,mass_tmp,rest_mass,r2_0);
        	         fflush(fsink1);}
        	 else if(ibh==1){
                 	fprintf(fsink2,"%1.10g\t %d %d %g %g %g\n",All.Time,P[target].ID,P[target].ID_child_number,mass_tmp,rest_mass,r2_1);
                 	fflush(fsink2);}
		 fclose(fsink1);
		 fclose(fsink2);
       	 	 return;
	 }
#endif

#ifdef IMPLICIT_SOLVER
	 //Compute error
         error=0;
         for(k=0;k<3;k++) {x_err[k]=P[target].Pos[k]-x_err[k]; if(fabs(x_err[k])>error) error=fabs(x_err[k]); }
//	 double r=sqrt(P[target].Pos[0]*P[target].Pos[0]+P[target].Pos[1]*P[target].Pos[1]+
//                P[target].Pos[2]*P[target].Pos[2]);

//	 if(r>MIN_REAL_NUMBER)
//	   error/=sqrt(P[target].Pos[0]*P[target].Pos[0]+P[target].Pos[1]*P[target].Pos[1]+
//		P[target].Pos[2]*P[target].Pos[2]);
	 iter++;
	 if(iter>IMPLICIT_MAX_NITER)
	 {
		printf("WARNING (drift): Iter >%d => Pos %2.15g %2.15g %2.15g; err %g %g %g / %g\n",IMPLICIT_MAX_NITER,
 			P[target].Pos[0],P[target].Pos[1],P[target].Pos[2],x_err[0],x_err[1],x_err[2],error);
		break;
	 }
	}while(error>IMPLICIT_SOLVER);
#endif

	//Update metric at new location
	set_metric_coefficients(P[target].Pos,P[target].Potential,&(CellP[target].gr));

        //Get primitive quantities (corrector)
        res_con2prim=gr_get_primitive_variables(&Cons,&Prims,&(CellP[target].gr), gr_limiter);
        if(res_con2prim<0) {printf("Warning: failed inversion with code %g - particle %d \n",res_con2prim, P[target].ID);endrun(9992407);}

	for(k=0;k<3;k++) v_new[k]=Prims.v[k];
	Vcoord(v_new,CellP[target].gr.beta,CellP[target].gr.alpha);

	//Now update particle velocity to be consistent with the metric change
	for(k=0;k<3;k++)
	{
	 P[target].Vel[k] += (v_new[k]-P[target].Vel[k]);
#ifdef MAGNETIC
         CellP[target].B[k] = Prims.B[k]*Vi*CellP[target].gr.gdet;
#endif
	}
	CellP[target].InternalEnergy += (Prims.u-CellP[target].InternalEnergy);
#ifdef GR_ENTROPY_SWITCH
	//Entropy Switch
	if(CellP[target].Entropy != 0){
		if(0.5*Prims.bsq >= 0.1*Prims.p){
			double EnergyDensFromEntropy;
			EnergyDensFromEntropy = CellP[target].Entropy * pow(Prims.rho0,GAMMA(target)-1) / (GAMMA(target)-1);
			CellP[target].InternalEnergy += EnergyDensFromEntropy - CellP[target].InternalEnergy; //Entropy Switch
		}
	}
#endif
	//Finally update predicted quantities
	res_con2prim=gr_get_primitive_variables(&ConsPred,&Prims,&(CellP[target].gr), gr_limiter);
        if(res_con2prim<0) {printf("Warning: failed inversion with code %g - particle %d \n",res_con2prim, P[target].ID);endrun(9992408);}
	Vcoord(Prims.v,CellP[target].gr.beta,CellP[target].gr.alpha);
	for(k=0;k<3;k++)
        {
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION>=5)
         CellP[target].ParticleVel[k] += Prims.v[k]-CellP[target].VelPred[k];
#endif
	 CellP[target].VelPred[k] += Prims.v[k] - CellP[target].VelPred[k];
#ifdef MAGNETIC
	 CellP[target].BPred[k] = Prims.B[k]*Vi*CellP[target].gr.gdet;
#endif
	}

	CellP[target].InternalEnergyPred+=Prims.u - CellP[target].InternalEnergyPred;
#ifdef GR_ENTROPY_SWITCH
	//Entropy Switch
	if(CellP[target].Entropy != 0){
		if(0.5*Prims.bsq >= 0.1*Prims.p){
			double EnergyDensFromEntropyPred;
			EnergyDensFromEntropyPred = CellP[target].EntropyPred * pow(Prims.rho0,GAMMA(target)-1) / (GAMMA(target)-1);
			CellP[target].InternalEnergyPred += EnergyDensFromEntropyPred - CellP[target].InternalEnergyPred; //Entropy switch
		}
	}
#endif

}

void do_predict_general_relativity(int target,double dt_hydrokick)
{
        int k;
	struct gr_primitives Prims;
        struct gr_conservatives ConsPred;
        double res_con2prim;
	//Get particle volume
	double Vi=P[target].Mass/CellP[target].Density;
	double mass_new = P[target].Mass;

#ifdef GR_TEST_PARTICLE
	mass_new=P[target].Mass=Vi=1.0;
#endif
#ifdef PARTICLE_EXCISION
//        check_boundary(target);
        if(P[target].Mass==0) return;
#endif

	//Fill the primitive variable structure first
        for(k=0;k<3;k++) Prims.v[k]=CellP[target].VelPred[k];
        Vphys(Prims.v,CellP[target].gr.beta,CellP[target].gr.alpha);
	int gr_limiter=0;
	double dpl[3]; for(k=0;k<3;k++) {dpl[k]=P[target].Pos[k];}
#if (GR_METRIC==2)
	for(k=0;k<3;k++) {dpl[k] = -P[target].Min_xyz_to_Sink[k];}
#endif
	if(dpl[0]*dpl[0]+dpl[1]*dpl[1]+dpl[2]*dpl[2] <= pow(GAMMA_LIMITER_RADIUS,2)) gr_limiter=1;
//	if(P[target].Pos[0]<LEFT_GAMMA_BOUNDARY||P[target].Pos[0]>RIGHT_GAMMA_BOUNDARY) gr_limiter=1;
	double Wlorentz=lorentz_factor(CellP[target].gr.gamma,Prims.v, gr_limiter);

        Prims.rho0 = P[target].Mass/Wlorentz/CellP[target].gr.gdet/Vi;
        Prims.u    = CellP[target].InternalEnergyPred;
        Prims.p    = gr_pressure_rho0_u(Prims.rho0,Prims.u);
#ifdef MAGNETIC
        for(k=0;k<3;k++) Prims.B[k]=CellP[target].BPred[k]/Vi/CellP[target].gr.gdet;
        Prims.bsq=Get_Gas_BField_RestFrame(target,Prims.B,Prims.v,Prims.b);
#endif

#ifdef HYDRO_MESHLESS_FINITE_VOLUME
	mass_new = fmax(mass_new+CellP[target].DtMass*dt_hydrokick,0.5*CellP[target].MassTrue);
#endif

	//Determine conservative variables
        gr_get_conservative_variables(&Prims,&ConsPred,Wlorentz,CellP[target].gr.gamma,CellP[target].gr.gdet);

        //We can now update predicted quantities
	double Vnew = Vi*exp(CellP[target].DtVolume*dt_hydrokick/Vi);

	//Add source terms
	gr_add_source_term(target,&Prims,&ConsPred,dt_hydrokick,0);
	if(P[target].Mass==0){return;}

	res_con2prim=gr_get_primitive_variables(&ConsPred,&Prims,&(CellP[target].gr), gr_limiter);
 	if(res_con2prim<0) {printf("Warning: failed inversion with code %g - particle %d - %g %g %g \n",res_con2prim, P[target].ID,P[target].Pos[0],P[target].Pos[1],P[target].Pos[2]);endrun(99924091);}


	ConsPred.D = mass_new/Vnew;
	ConsPred.tau = (ConsPred.tau*Vi + (CellP[target].DtEnergy*dt_hydrokick))/Vnew;
	for(k=0;k<3;k++) ConsPred.S[k] = (ConsPred.S[k]*Vi + (CellP[target].DtMomentum[k]*dt_hydrokick))/Vnew;

#ifdef MAGNETIC
	for(k=0;k<3;k++) ConsPred.B[k] = (ConsPred.B[k]*Vi + (CellP[target].DtB[k]*dt_hydrokick))/Vnew;
#endif

	//Find predicted primitives
	res_con2prim=gr_get_primitive_variables(&ConsPred,&Prims,&(CellP[target].gr), gr_limiter);
 	if(res_con2prim<0) {printf("Warning: failed inversion with code %g - particle %d \n",res_con2prim, P[target].ID);endrun(9992409);}

	//Convert vel back to dx/dt
        Vcoord(Prims.v,CellP[target].gr.beta,CellP[target].gr.alpha);
	for(k=0;k<3;k++) CellP[target].VelPred[k] += Prims.v[k]-CellP[target].VelPred[k];
	CellP[target].InternalEnergyPred += (Prims.u-CellP[target].InternalEnergyPred);

#if(GR_METRIC==2)
	if(CellP[target].InternalEnergyPred>=1.0){
		if(sqrt(dpl[0]*dpl[0]+dpl[1]*dpl[1]+dpl[2]*dpl[2])>5){
			printf("WARNING: u>1.0 for particle %d at %g from BH. uPred = %g\n", P[target].ID,sqrt(dpl[0]*dpl[0]+dpl[1]*dpl[1]+dpl[2]*dpl[2]),CellP[target].InternalEnergyPred);
		}
		P[target].Mass=0;
	}
#endif

#ifdef GR_ENTROPY_SWITCH
	//Entropy Switch
	if(CellP[target].Entropy != 0){
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
		CellP[target].EntropyPred = CellP[target].EntropyPred + (CellP[target].DtEntropy*dt_hydrokick)/Vi;
#endif
		if(0.5*Prims.bsq >= 0.1*Prims.p){
			double EnergyDensFromEntropy;
			EnergyDensFromEntropy = CellP[target].EntropyPred * pow(Prims.rho0,GAMMA(target)-1) / (GAMMA(target)-1);
			CellP[target].InternalEnergyPred += (EnergyDensFromEntropy-CellP[target].InternalEnergyPred); //Entropy Switch
		}
	}
#endif

	//Now zero out the rates of change, to avoid double counting accelerations
	for(k=0;k<3;k++) CellP[target].HydroAccel[k] = 0;
	CellP[target].DtInternalEnergy = 0;
}


#endif

