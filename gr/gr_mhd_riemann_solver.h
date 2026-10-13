#define DEDNER_IMPLICIT_LIMITER_GR 0.75
#define C_G All.CourantFac //should be between 0 and 1
#define BX_MIN 1e-4
//#define DEBUG

#ifdef MAGNETIC
void compute_lambda_hll(double v[],double cs,double va,double g_nn,double beta_n, double alpha,double Wfac,double lam[2], int gr_limiter)
{
	//Get u_n (Euclidean projection)
	double cs2=cs*cs*(1-va*va)+va*va;
        double v2max = 1.-MIN_REAL_NUMBER;
	if(gr_limiter==1)
        	v2max = 1.-1./GR_GAMMA_LIMITER/GR_GAMMA_LIMITER;
	long double v2=1.-1./Wfac/Wfac;
	if(v2<0) v2=0.;
	double vtilde_n =alpha*v[0]-beta_n;
	double a = Wfac*Wfac*(cs2-1) - cs2;
	double b = -2*(Wfac*Wfac*(cs2-1)*vtilde_n + (cs2*beta_n));
	double c = cs2*(alpha*alpha*g_nn - beta_n*beta_n) + (Wfac*Wfac*vtilde_n*vtilde_n*(cs2-1)); 
	long double discr = sqrtl(b*b - 4*a*c);
	lam[0] = DMIN((-b - discr)/(2*a),v2max);
	lam[1] = DMAX((-b + discr)/(2*a),-v2max);
}

void GR_MHD_Riemann_solver(struct Input_vec_Riemann Riemann_vec,struct Riemann_outputs *Riemann_out, int gr_limiter)
{
	int k;
	//First rotate beta along n and the metric tensot
        double beta_n=Riemann_vec.gr.beta[0];

	//We are already aligned with x, hence we take the xx component
	double gnn=Riemann_vec.gr.ugamma[0];

	//Compute the frame line velocity
	double v_frame_line = Riemann_vec.gr.alpha*Riemann_vec.v_frame[0] - beta_n;

	//Get conservative
	struct gr_conservatives cons_L,cons_R;
	struct gr_primitives prims;
	double W_l,W_r,vL[3],vR[3],tmpvec[3];
	double vAL,vAR;
	double bL[4],bR[4];
	double lam_l[2];
	double lam_r[2];
	double S_L,S_R;
	double beta_cov[3];
	double BxL,BxR;
	raise_lower_indices(Riemann_vec.gr.gamma,Riemann_vec.gr.beta,beta_cov);

	//Bx continuity
	BxL = Riemann_vec.L.B[0];
	BxR = Riemann_vec.R.B[0];
	double B_normal_corrected = 0.5*(BxL + BxR);

#ifdef DIVBCLEANING_GR
	/* (B_perp, psi) subsystem, Route B (gr/ROUTE_B_cleaning_notes.pdf, eqs. 13-14): Rusanov states with characteristic speed alpha*c_h*sqrt(gamma^nn).
	   L.phi/R.phi carry Dedner's psi = c_h*psihat of each side (multiplied in hydro_core_meshless.h); c_h is the face value (max over the pair).
	   Powers of c_h: 1/c_h on the psi jump in B*, c_h on the B jump in psi*. For c_h = 1 this is Fedrigo & Lupi (2025) eq. (32). */
	double c_h = Riemann_vec.c_h;
	double ch2 = c_h*c_h;
	double phiL = Riemann_vec.L.phi, phiR = Riemann_vec.R.phi;
	double corr_norm = 1.0, corr_p = 0.;
	if(c_h > 0) {corr_p = 0.5*(phiL - phiR)*sqrt(gnn)/c_h;}
	double corr_p_abs = fabs(corr_p);
	double corr_b_abs = DEDNER_IMPLICIT_LIMITER_GR * fabs(B_normal_corrected);
	if(corr_p_abs > corr_b_abs) {corr_norm *= corr_b_abs/corr_p_abs; if(stamp==1){ printf("Correzione nel Riemann Solver AVVENUTA, fattore corr_norm %g \n", corr_norm);}}
	B_normal_corrected += corr_norm * corr_p;
	Riemann_out->phi_normal_mean = 0.5 * corr_norm * (phiL + phiR);
	Riemann_out->phi_normal_db = corr_norm * 0.5*(Riemann_vec.L.B[0] - Riemann_vec.R.B[0]) * c_h / sqrt(gnn);
#endif
	BxL = B_normal_corrected;
	BxR = B_normal_corrected;
	Riemann_out->B_normal_corrected = Riemann_vec.gr.gdet*B_normal_corrected;

        //Left state
        prims.rho0=Riemann_vec.L.rho;
        for(k=0;k<3;k++) prims.v[k]=vL[k]=Riemann_vec.L.v[k];
        prims.u=Riemann_vec.L.u;
        prims.p=Riemann_vec.L.p;
        W_l=lorentz_factor(Riemann_vec.gr.gamma,prims.v, gr_limiter);

        //Now get Rest Frame b field
        prims.bsq=0;
        raise_lower_indices(Riemann_vec.gr.gamma,prims.v,vL); //Get covariant
        for(k=1;k<3;k++) {prims.bsq += Riemann_vec.L.B[k]*vL[k]; prims.B[k]=Riemann_vec.L.B[k];}
        prims.bsq += BxL*vL[0];
        prims.B[0] = BxL;

        prims.b[0] = W_l*prims.bsq;
        for(k=0;k<3;k++) tmpvec[k] = prims.B[k]/W_l + W_l*prims.bsq*(prims.v[k]-(Riemann_vec.gr.beta[k]/Riemann_vec.gr.alpha));
        //Now lower indices for j components
        raise_lower_indices(Riemann_vec.gr.gamma,tmpvec,prims.b+1);
        for(k=1;k<=3; k++) {prims.b[k] += beta_cov[k-1]/Riemann_vec.gr.alpha*prims.b[0];}
        prims.bsq = prims.bsq*prims.bsq + vector_sqnorm(Riemann_vec.gr.gamma,prims.B)/W_l/W_l;

        vAL = sqrt(prims.bsq/(prims.rho0*(1+prims.u)+(prims.p+prims.bsq)));

        for(k=0;k<4;k++) bL[k]=prims.b[k];
        gr_get_conservative_variables(&prims,&cons_L,W_l,Riemann_vec.gr.gamma,Riemann_vec.gr.gdet);

        //Right state
        prims.rho0=Riemann_vec.R.rho;
        for(k=0;k<3;k++) prims.v[k]=vR[k]=Riemann_vec.R.v[k];
        prims.u=Riemann_vec.R.u;
        prims.p=Riemann_vec.R.p;
        W_r=lorentz_factor(Riemann_vec.gr.gamma,prims.v, gr_limiter);

        //Now get Rest Frame b field
        prims.bsq=0;
        raise_lower_indices(Riemann_vec.gr.gamma,prims.v,vR); //Get covariant
        for(k=1;k<3;k++) {prims.bsq += Riemann_vec.R.B[k]*vR[k]; prims.B[k]=Riemann_vec.R.B[k];}
        prims.bsq += BxR*vR[0];
        prims.B[0] = BxR;

        prims.b[0] = W_r*prims.bsq;
        for(k=0;k<3;k++) tmpvec[k] = prims.B[k]/W_r + W_r*prims.bsq*(prims.v[k]-(Riemann_vec.gr.beta[k]/Riemann_vec.gr.alpha));
        //Now lower indices for j components
        raise_lower_indices(Riemann_vec.gr.gamma,tmpvec,prims.b+1);
        for(k=1;k<=3; k++) {prims.b[k] += beta_cov[k-1]/Riemann_vec.gr.alpha*prims.b[0];}
        prims.bsq = prims.bsq*prims.bsq + vector_sqnorm(Riemann_vec.gr.gamma,prims.B)/W_r/W_r;

        vAR = sqrt(prims.bsq/(prims.rho0*(1+prims.u)+(prims.p+prims.bsq)));

        for(k=0;k<4;k++) bR[k]=prims.b[k];
        gr_get_conservative_variables(&prims,&cons_R,W_r,Riemann_vec.gr.gamma,Riemann_vec.gr.gdet);

        for(k=0;k<3;k++) {vL[k] = Riemann_vec.gr.alpha*Riemann_vec.L.v[k]-Riemann_vec.gr.beta[k]; vR[k]=Riemann_vec.gr.alpha*Riemann_vec.R.v[k]-Riemann_vec.gr.beta[k];}

        //Compute lambda_hll with the corrected Bx
        lam_l[0]=lam_r[0]=lam_l[1]=lam_r[1]=S_L=S_R=0;

        compute_lambda_hll(Riemann_vec.L.v,Riemann_vec.L.cs,vAL,gnn,beta_n,Riemann_vec.gr.alpha,W_l,lam_l, gr_limiter);
        compute_lambda_hll(Riemann_vec.R.v,Riemann_vec.R.cs,vAR,gnn,beta_n,Riemann_vec.gr.alpha,W_r,lam_r, gr_limiter);

        S_L = DMIN(lam_r[1],lam_l[1]); //Lambda minus
        S_R = DMAX(lam_r[0],lam_l[0]); //Lambda plus


	Riemann_out->S_M = 0.5*(S_L+S_R);
	Riemann_out->P_M = 0.5*(cons_L.p+cons_R.p);

	struct gr_conservatives fluxL,fluxR;

	//Left flux
	fluxL.D    = cons_L.D  *vL[0];
        fluxL.tau  = cons_L.tau*vL[0] + Riemann_vec.gr.gdet*cons_L.p*(vL[0]+beta_n) - Riemann_vec.gr.alpha*bL[0]*cons_L.B[0]/W_l;

        for(k=0;k<3;k++)
	{
                fluxL.S[k] = cons_L.S[k]*vL[0] - Riemann_vec.gr.alpha*bL[k+1]*cons_L.B[0]/W_l;
		fluxL.B[k] = cons_L.B[k]*vL[0] - cons_L.B[0]*vL[k];
	}
	//Null Bx Flux
	fluxL.B[0] = 0;
	fluxL.S[0] += Riemann_vec.gr.alpha*Riemann_vec.gr.gdet*cons_L.p;
#ifdef DIVBCLEANING_GR
	fluxL.phi =  cons_L.D * vL[0] * phiL;
#endif
	//Right flux
	fluxR.D    = cons_R.D  *vR[0];
 	fluxR.tau  = cons_R.tau*vR[0] + Riemann_vec.gr.gdet*cons_R.p*(vR[0]+beta_n) - Riemann_vec.gr.alpha*bR[0]*cons_R.B[0]/W_r;

        for(k=0;k<3;k++)
	{
                fluxR.S[k] = cons_R.S[k]*vR[0] - Riemann_vec.gr.alpha*bR[k+1]*cons_R.B[0]/W_r;
		fluxR.B[k] = cons_R.B[k]*vR[0] - cons_R.B[0]*vR[k];
	}
	//Null Bx Flux
	fluxR.B[0]=0;
	fluxR.S[0] += Riemann_vec.gr.alpha*Riemann_vec.gr.gdet*cons_R.p;
#ifdef DIVBCLEANING_GR
	fluxR.phi =  cons_R.D * vR[0] * phiR;
#endif

	//OUTPUT

#ifdef HYDRO_MESHLESS_FINITE_VOLUME
	if(S_L>v_frame_line)
	{
	  Riemann_out->Fluxes.rho = fluxL.D - cons_L.D*v_frame_line;
	  Riemann_out->Fluxes.p = fluxL.tau - cons_L.tau*v_frame_line;
	  for(k=0;k<3;k++){
		Riemann_out->Fluxes.v[k] = fluxL.S[k] - cons_L.S[k]*v_frame_line;
		Riemann_out->Fluxes.B[k] = fluxL.B[k] - cons_L.B[k]*v_frame_line;
	  }
#ifdef DIVBCLEANING_GR
	  Riemann_out->Fluxes.phi = fluxL.phi - cons_L.D*phiL*v_frame_line;
#endif

	}
	else if(S_R<v_frame_line)
	{
	  Riemann_out->Fluxes.rho = fluxR.D - cons_R.D*v_frame_line;
	  Riemann_out->Fluxes.p = fluxR.tau - cons_R.tau*v_frame_line;
	  for(k=0;k<3;k++){
		Riemann_out->Fluxes.v[k] = fluxR.S[k] - cons_R.S[k]*v_frame_line;
		Riemann_out->Fluxes.B[k] = fluxR.B[k] - cons_R.B[k]*v_frame_line;
	  }
#ifdef DIVBCLEANING_GR
	  Riemann_out->Fluxes.phi = fluxR.phi - cons_R.D*phiR*v_frame_line;
#endif
	}
	else if(S_L<=v_frame_line && v_frame_line<=S_R)
	{
#endif
	 //Get HLL intermediate state and fluxes
	 struct gr_conservatives U_hll;
         double facLR=1./(S_R-S_L);
	 U_hll.D = (S_R*cons_R.D-S_L*cons_L.D + fluxL.D - fluxR.D)*facLR;
	 U_hll.tau = (S_R*cons_R.tau-S_L*cons_L.tau + fluxL.tau - fluxR.tau)*facLR;
	 double Bx_hll;
	 for(k=0;k<3;k++)
	 {
	 	U_hll.S[k] = (S_R*cons_R.S[k]-S_L*cons_L.S[k] + fluxL.S[k] - fluxR.S[k])*facLR;
	 	U_hll.B[k] = (S_R*cons_R.B[k]-S_L*cons_L.B[k] + fluxL.B[k] - fluxR.B[k])*facLR;
	 }
#ifdef DIVBCLEANING_GR
	 U_hll.phi = (sqrt(ch2)*cons_R.D*phiR-sqrt(ch2)*cons_L.D*phiL + fluxL.phi - fluxR.phi)*facLR;
#endif
#ifdef HYDRO_MESHLESS_FINITE_MASS
	 //This is the HLL equivalent of the zero-mass-flux face velocity
	 double w_bck=v_frame_line;
	 //Fluxrho = (SR*FL - SL*FR + SL*SR*(UR-UL))/(SR-SL)
	 //Uhll    = (SR*UR - SL*UL + FL - FR)/(SR-SL)
	 v_frame_line = (S_R*cons_L.D*vL[0]-S_L*cons_R.D*vR[0] + S_L*S_R*(cons_R.D - cons_L.D))*facLR/U_hll.D;
	 Riemann_out->S_M=v_frame_line;
	 Riemann_out->Fluxes.rho=0;
	 //We store the density in the intermediate state, for face correction
	 //Riemann_out->Fluxes.rho=U_hll.D;
#else
	 Riemann_out->Fluxes.rho = (S_R*fluxL.D-S_L*fluxR.D + S_L*S_R*(cons_R.D - cons_L.D))*facLR - v_frame_line*U_hll.D;
#ifdef DIVBCLEANING_GR
	 Riemann_out->Fluxes.phi = (sqrt(ch2)*fluxL.phi-sqrt(ch2)*fluxR.phi + ch2*(cons_R.D*phiR - cons_L.D*phiL))*facLR - v_frame_line*U_hll.phi;
#endif
#endif
	 Riemann_out->Fluxes.p = (S_R*fluxL.tau-S_L*fluxR.tau + S_L*S_R*(cons_R.tau - cons_L.tau))*facLR - v_frame_line*U_hll.tau;

	 for(k=0;k<3;k++)
	 {
	 	Riemann_out->Fluxes.v[k] = (S_R*fluxL.S[k]-S_L*fluxR.S[k] + S_L*S_R*(cons_R.S[k] - cons_L.S[k]))*facLR - v_frame_line*U_hll.S[k];
	 	Riemann_out->Fluxes.B[k] = (S_R*fluxL.B[k]-S_L*fluxR.B[k] + S_L*S_R*(cons_R.B[k] - cons_L.B[k]))*facLR - v_frame_line*U_hll.B[k];
	 }
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
	}
	else
	{
		printf("Error here v_frame_line %g - SL %g - SR %g\n", v_frame_line,S_L,S_R);
		printf("Right: rho %g - v %g %g - u %g - p %g - B %g %g \n", Riemann_vec.R.rho,Riemann_vec.R.v[0],Riemann_vec.R.v[1],Riemann_vec.R.u, Riemann_vec.R.p, Riemann_vec.R.B[0],Riemann_vec.R.B[1]);
		printf("Left: rho %g - v %g %g - u %g - p %g - B %g %g \n", Riemann_vec.L.rho,Riemann_vec.L.v[0],Riemann_vec.L.v[1],Riemann_vec.L.u, Riemann_vec.L.p, Riemann_vec.L.B[0],Riemann_vec.L.B[1]);
		endrun(47372);
	}
#endif

/*-----------------------------------------*/

#ifdef DEBUG
	printf("My flux is %2.12g %12.5e %2.12g %2.12g %2.12g %2.12g %2.12g %2.12g\n",Riemann_out->Fluxes.rho,Riemann_out->Fluxes.v[0],Riemann_out->Fluxes.v[1],
			Riemann_out->Fluxes.v[2],Riemann_out->Fluxes.p,Riemann_out->Fluxes.B[0],Riemann_out->Fluxes.B[1],Riemann_out->Fluxes.B[2]);
#endif
#ifdef SAVE_FACE_DENSITY
    if((Riemann_out->S_M==0) || ((Riemann_out->P_M<=0)&&(!isnan(Riemann_out->P_M))))
    {
        Riemann_out->Face_Density = 0.5*(Riemann_vec.L.rho+Riemann_vec.R.rho);
    } else {
        Riemann_out->Face_Density = 0.5 * (Riemann_vec.L.rho*(S_L-v_line_L)/(S_L-Riemann_out->S_M) +
                                           Riemann_vec.R.rho*(S_R-v_line_R)/(S_R-Riemann_out->S_M));
        if(Riemann_out->S_M==S_L) {Riemann_out->Face_Density=Riemann_vec.L.rho;}
        if(Riemann_out->S_M==S_R) {Riemann_out->Face_Density=Riemann_vec.R.rho;}
    }
#endif

}


void GR_MHD_Riemann_solver_GForce(struct Input_vec_Riemann Riemann_vec,struct Riemann_outputs *Riemann_out, int gr_limiter)
{
	if(stamp==1){printf("RIEMANN SOLVER\n\n");}
	int k;
	double res_con2prim=0;
	//First rotate beta along n and the metric tensot
        double beta_n=Riemann_vec.gr.beta[0];

	//We are already aligned with x, hence we take the xx component
	double gnn=Riemann_vec.gr.ugamma[0];

	//Compute the frame line velocity
	double v_frame_line = Riemann_vec.gr.alpha*Riemann_vec.v_frame[0] - beta_n;

	//Get conservative
	struct gr_conservatives cons_L,cons_R;
	struct gr_primitives prims;
	double W_l,W_r,vL[3],vR[3],tmpvec[3];
	double vAL,vAR;
	double bL[4],bR[4];
	double lam_l[2];
	double lam_r[2];
	double S_L,S_R;
	double beta_cov[3];
	double omega_g = 1/(1+C_G);
	double BxL, BxR;
	raise_lower_indices(Riemann_vec.gr.gamma,Riemann_vec.gr.beta,beta_cov);

	//Bx continuity
	BxL = Riemann_vec.L.B[0];
	BxR = Riemann_vec.R.B[0];
	double B_normal_corrected = 0.5*(BxL + BxR);

#ifdef DIVBCLEANING_GR
	/* (B_perp, psi) subsystem, Route B (gr/ROUTE_B_cleaning_notes.pdf, eqs. 13-14): Rusanov states with characteristic speed alpha*c_h*sqrt(gamma^nn).
	   L.phi/R.phi carry Dedner's psi = c_h*psihat of each side (multiplied in hydro_core_meshless.h); c_h is the face value (max over the pair).
	   Powers of c_h: 1/c_h on the psi jump in B*, c_h on the B jump in psi*. For c_h = 1 this is Fedrigo & Lupi (2025) eq. (32). */
	double c_h = Riemann_vec.c_h;
	double ch2 = c_h*c_h;
	double phiL = Riemann_vec.L.phi, phiR = Riemann_vec.R.phi;
	double corr_norm = 1.0, corr_p = 0.;
	if(c_h > 0) {corr_p = 0.5*(phiL - phiR)*sqrt(gnn)/c_h;}
	double corr_p_abs = fabs(corr_p);
	double corr_b_abs = DEDNER_IMPLICIT_LIMITER_GR * fabs(B_normal_corrected);
	if(corr_p_abs > corr_b_abs) {corr_norm *= corr_b_abs/corr_p_abs; if(stamp==1){ printf("Correzione nel Riemann Solver AVVENUTA, fattore corr_norm %g \n", corr_norm);}}
	B_normal_corrected += corr_norm * corr_p;
	Riemann_out->phi_normal_mean = 0.5 * corr_norm * (phiL + phiR);
	Riemann_out->phi_normal_db = corr_norm * 0.5*(Riemann_vec.L.B[0] - Riemann_vec.R.B[0]) * c_h / sqrt(gnn);
#endif
	BxL = B_normal_corrected;
	BxR = B_normal_corrected;
	Riemann_out->B_normal_corrected = Riemann_vec.gr.gdet*B_normal_corrected;

	//Left state
	prims.rho0=Riemann_vec.L.rho;
	for(k=0;k<3;k++) prims.v[k]=vL[k]=Riemann_vec.L.v[k];
	prims.u=Riemann_vec.L.u;
	prims.p=Riemann_vec.L.p;
	W_l=lorentz_factor(Riemann_vec.gr.gamma,prims.v, gr_limiter);

	//Now get Rest Frame b field
	prims.bsq=0;
	raise_lower_indices(Riemann_vec.gr.gamma,prims.v,vL); //Get covariant
	for(k=1;k<3;k++) {prims.bsq += Riemann_vec.L.B[k]*vL[k]; prims.B[k]=Riemann_vec.L.B[k];}
	prims.bsq += BxL*vL[0];
	prims.B[0] = BxL;

	prims.b[0] = W_l*prims.bsq;
	for(k=0;k<3;k++) tmpvec[k] = prims.B[k]/W_l + W_l*prims.bsq*(prims.v[k]-(Riemann_vec.gr.beta[k]/Riemann_vec.gr.alpha));
	//Now lower indices for j components
	raise_lower_indices(Riemann_vec.gr.gamma,tmpvec,prims.b+1);
	for(k=1;k<=3;k++) {prims.b[k] += beta_cov[k-1]/Riemann_vec.gr.alpha*prims.b[0];}
	prims.bsq = prims.bsq*prims.bsq + vector_sqnorm(Riemann_vec.gr.gamma,prims.B)/W_l/W_l;

	vAL = sqrt(prims.bsq/(prims.rho0*(1+prims.u)+(prims.p+prims.bsq)));

	for(k=0;k<4;k++) bL[k]=prims.b[k]; //bL = (b^0, b_i)
	gr_get_conservative_variables(&prims,&cons_L,W_l,Riemann_vec.gr.gamma,Riemann_vec.gr.gdet);

	if(stamp==1) printf("Prims L: rho %g - u %g - v %g %g %g - B %g %g %g \n",prims.rho0, prims.u, prims.v[0], prims.v[1],prims.v[2], prims.B[0], prims.B[1], prims.B[2]);
	if(stamp==1) printf("Cons L: D %g - tau %g - S %g %g %g - B %g %g %g \n",cons_L.D, cons_L.tau, cons_L.S[0],cons_L.S[1],cons_L.S[2],cons_L.B[0],cons_L.B[1],cons_L.B[2]);

	//Right state
	prims.rho0=Riemann_vec.R.rho;
	for(k=0;k<3;k++) prims.v[k]=vR[k]=Riemann_vec.R.v[k];
	prims.u=Riemann_vec.R.u;
	prims.p=Riemann_vec.R.p;
	W_r=lorentz_factor(Riemann_vec.gr.gamma,prims.v, gr_limiter);

	//Now get Rest Frame b field
	prims.bsq=0;
	raise_lower_indices(Riemann_vec.gr.gamma,prims.v,vR); //Get covariant
	for(k=1;k<3;k++) {prims.bsq += Riemann_vec.R.B[k]*vR[k]; prims.B[k]=Riemann_vec.R.B[k];}
	prims.bsq += BxR*vR[0];
	prims.B[0] = BxR;

	prims.b[0] = W_r*prims.bsq;
	for(k=0;k<3;k++) tmpvec[k] = prims.B[k]/W_r + W_r*prims.bsq*(prims.v[k]-(Riemann_vec.gr.beta[k]/Riemann_vec.gr.alpha));
	//Now lower indices for j components
	raise_lower_indices(Riemann_vec.gr.gamma,tmpvec,prims.b+1);
	for(k=1;k<=3;k++) {prims.b[k] += beta_cov[k-1]/Riemann_vec.gr.alpha*prims.b[0];}
	prims.bsq = prims.bsq*prims.bsq + vector_sqnorm(Riemann_vec.gr.gamma,prims.B)/W_r/W_r;

	vAR = sqrt(prims.bsq/(prims.rho0*(1+prims.u)+(prims.p+prims.bsq)));

	for(k=0;k<4;k++) bR[k]=prims.b[k]; //bR = (b^0, b_i)
	gr_get_conservative_variables(&prims,&cons_R,W_r,Riemann_vec.gr.gamma,Riemann_vec.gr.gdet);

	if(stamp==1) printf("Prims R: rho %g - u %g - v %g %g %g - B %g %g %g \n",prims.rho0, prims.u, prims.v[0], prims.v[1],prims.v[2], prims.B[0], prims.B[1], prims.B[2]);
	if(stamp==1) printf("Cons R: D %g - tau %g - S %g %g %g - B %g %g %g \n",cons_R.D, cons_R.tau, cons_R.S[0],cons_R.S[1],cons_R.S[2],cons_R.B[0],cons_R.B[1],cons_R.B[2]);

	//Recompute coordinate velocities with upper indices
	for(k=0;k<3;k++) {vL[k] = Riemann_vec.gr.alpha*Riemann_vec.L.v[k]-Riemann_vec.gr.beta[k]; vR[k]=Riemann_vec.gr.alpha*Riemann_vec.R.v[k]-Riemann_vec.gr.beta[k];}

	//Compute lambda_hll with the corrected Bx
	lam_l[0]=lam_r[0]=lam_l[1]=lam_r[1]=S_L=S_R=0;

	compute_lambda_hll(Riemann_vec.L.v,Riemann_vec.L.cs,vAL,gnn,beta_n,Riemann_vec.gr.alpha,W_l,lam_l, gr_limiter);
	compute_lambda_hll(Riemann_vec.R.v,Riemann_vec.R.cs,vAR,gnn,beta_n,Riemann_vec.gr.alpha,W_r,lam_r, gr_limiter);

	S_L = DMIN(lam_r[1],lam_l[1]); //Lambda minus
	S_R = DMAX(lam_r[0],lam_l[0]); //Lambda plus

	double S_MAX = 1/DMAX(fabs(S_L),fabs(S_R));

	if(stamp==1) printf("Waves speed SL %g - SR %g - SMAX %g - tau %g \n",S_L,S_R,DMAX(fabs(S_L),fabs(S_R)), S_MAX);

	Riemann_out->S_M = 0.5*(S_L+S_R);
	Riemann_out->P_M = 0.5*(cons_L.p+cons_R.p);

	struct gr_conservatives fluxL,fluxR;

        //Left flux
	fluxL.D    = cons_L.D  *vL[0];
	fluxL.tau  = cons_L.tau*vL[0] + Riemann_vec.gr.gdet*cons_L.p*(vL[0]+beta_n) - Riemann_vec.gr.alpha*bL[0]*cons_L.B[0]/W_l;

	for(k=0;k<3;k++)
	{
		fluxL.S[k] = cons_L.S[k]*vL[0] - Riemann_vec.gr.alpha*bL[k+1]*cons_L.B[0]/W_l;
		fluxL.B[k] = cons_L.B[k]*vL[0] - cons_L.B[0]*vL[k];
	}
	//Null Bx Flux
	fluxL.B[0] = 0;
	fluxL.S[0] += Riemann_vec.gr.alpha*Riemann_vec.gr.gdet*cons_L.p;

#ifdef DIVBCLEANING_GR
	fluxL.phi =  cons_L.D * vL[0] * phiL;
#endif

	if(stamp==1) printf("Flussi L: D %g - tau %g - S %g %g %g - B %g %g %g \n",fluxL.D,fluxL.tau,fluxL.S[0],fluxL.S[1],fluxL.S[2],fluxL.B[0],fluxL.B[1],fluxL.B[2]);

	//Right flux
	fluxR.D    = cons_R.D  *vR[0];
	fluxR.tau  = cons_R.tau*vR[0] + Riemann_vec.gr.gdet*cons_R.p*(vR[0]+beta_n) - Riemann_vec.gr.alpha*bR[0]*cons_R.B[0]/W_r;

	for(k=0;k<3;k++)
	{
		fluxR.S[k] = cons_R.S[k]*vR[0] - Riemann_vec.gr.alpha*bR[k+1]*cons_R.B[0]/W_r;
		fluxR.B[k] = cons_R.B[k]*vR[0] - cons_R.B[0]*vR[k];
	}
	//Null Bx Flux
	fluxR.B[0]=0;
	fluxR.S[0] += Riemann_vec.gr.alpha*Riemann_vec.gr.gdet*cons_R.p;

#ifdef DIVBCLEANING_GR
	fluxR.phi =  cons_R.D * vR[0] * phiR;
#endif

	if(stamp==1) printf("Flussi R: D %g - tau %g - S %g %g %g - B %g %g %g \n",fluxR.D,fluxR.tau,fluxR.S[0],fluxR.S[1],fluxR.S[2],fluxR.B[0],fluxR.B[1],fluxR.B[2]);

	//Lax-Friederichs (LF) fluxes
	struct gr_conservatives fluxLF;

	fluxLF.D = 0.5*(fluxR.D + fluxL.D) + 0.5*(cons_L.D - cons_R.D)/S_MAX;
	fluxLF.tau = 0.5*(fluxR.tau + fluxL.tau) + 0.5*(cons_L.tau - cons_R.tau)/S_MAX;
	for(k=0;k<3;k++)
	{
		fluxLF.S[k] = 0.5*(fluxR.S[k] + fluxL.S[k]) + 0.5*(cons_L.S[k] - cons_R.S[k])/S_MAX;
		fluxLF.B[k] = 0.5*(fluxR.B[k] + fluxL.B[k]) + 0.5*(cons_L.B[k] - cons_R.B[k])/S_MAX;
	}
	fluxLF.B[0] = 0;
#ifdef DIVBCLEANING_GR
	fluxLF.phi = 0.5*(fluxR.phi + fluxL.phi) + 0.5*sqrt(ch2)*(cons_L.D*phiL - cons_R.D*phiR);
#endif

	if(stamp==1) printf("Flussi Lax-Friederichs: D %g - tau %g - S %g %g %g - B %g %g %g \n",fluxLF.D,fluxLF.tau, fluxLF.S[0],fluxLF.S[1],fluxLF.S[2], fluxLF.B[0],fluxLF.B[1],fluxLF.B[2]);

	//Lax-Wendroff (LW) conservatives
	struct gr_conservatives cons_LW;

	cons_LW.D = 0.5*(cons_L.D + cons_R.D) + 0.5*S_MAX*(fluxL.D - fluxR.D);
	cons_LW.tau = 0.5*(cons_L.tau + cons_R.tau) + 0.5*S_MAX*(fluxL.tau - fluxR.tau);
	for(k=0;k<3;k++)
	{
		cons_LW.S[k] = 0.5*(cons_L.S[k] + cons_R.S[k]) + 0.5*S_MAX*(fluxL.S[k] - fluxR.S[k]);
		cons_LW.B[k] = 0.5*(cons_L.B[k] + cons_R.B[k]) + 0.5*S_MAX*(fluxL.B[k] - fluxR.B[k]);
	}
	cons_LW.B[0] = Riemann_vec.gr.gdet*B_normal_corrected;
#ifdef DIVBCLEANING_GR
	cons_LW.phi = 0.5*(cons_L.D*phiL + cons_R.D*phiR) + 0.5*(fluxL.phi - fluxR.phi)/sqrt(ch2);
#endif

	if(stamp==1) printf("Cons Lax-Wendroff: D %g - tau %g - S %g %g %g - B %g %g %g \n",cons_LW.D, cons_LW.tau, cons_LW.S[0],cons_LW.S[1],cons_LW.S[2], cons_LW.B[0],cons_LW.B[1],cons_LW.B[2]);

	//Lax-Wendroff (LW) primitives
//	struct gr_primitives prims;
	double W_lw;
	double v_coord_LW[3];
	double v_LW_cov[3];
	res_con2prim=gr_get_primitive_variables(&cons_LW,&prims,&(Riemann_vec.gr), gr_limiter);
	if(res_con2prim<0) {printf("Warning: failed inversion with code %g \n",res_con2prim); endrun(47373);}
	for(k=0;k<3;k++){v_coord_LW[k] = (Riemann_vec.gr.alpha*prims.v[k]-Riemann_vec.gr.beta[k]);}
	W_lw = lorentz_factor(Riemann_vec.gr.gamma,prims.v, gr_limiter);

	if(stamp==1) printf("Prims LW: rho %g - u %g - v %g %g %g - B %g %g %g - bsq %g \n",prims.rho0, prims.u, prims.v[0], prims.v[1],prims.v[2], prims.B[0], prims.B[1], prims.B[2],prims.bsq);

	prims.bsq=0;
	raise_lower_indices(Riemann_vec.gr.gamma,prims.v,v_LW_cov); //Get covariant
	for(k=0;k<3;k++) {prims.bsq += prims.B[k]*v_LW_cov[k];}
	prims.b[0] = W_lw*prims.bsq;
	for(k=0;k<3;k++) tmpvec[k] = prims.B[k]/W_lw + W_lw*prims.bsq*(prims.v[k]-(Riemann_vec.gr.beta[k]/Riemann_vec.gr.alpha));
        //Now lower indices for j components
	raise_lower_indices(Riemann_vec.gr.gamma,tmpvec,prims.b+1);
	for(k=1;k<=3;k++) {prims.b[k] += beta_cov[k-1]/Riemann_vec.gr.alpha*prims.b[0];}
	prims.bsq = prims.bsq*prims.bsq + vector_sqnorm(Riemann_vec.gr.gamma,prims.B)/W_lw/W_lw;

	//Adding the magnetic pressure to the Lax-Wendroff P*
	cons_LW.p = prims.p + 0.5*prims.bsq;

	//Lax-Wendroff (LW) fluxes
	struct gr_conservatives fluxLW;

	fluxLW.D    = cons_LW.D  * v_coord_LW[0];
	fluxLW.tau  = cons_LW.tau * v_coord_LW[0] + Riemann_vec.gr.gdet*cons_LW.p*(v_coord_LW[0]+beta_n) - Riemann_vec.gr.alpha*prims.b[0]*cons_LW.B[0]/W_lw;

	for(k=0;k<3;k++)
	{
		fluxLW.S[k] = cons_LW.S[k]*v_coord_LW[0] - Riemann_vec.gr.alpha*prims.b[k+1]*cons_LW.B[0]/W_lw;
		fluxLW.B[k] = cons_LW.B[k]*v_coord_LW[0] - cons_LW.B[0]*v_coord_LW[k];
	}
	//Null Bx Flux
	fluxLW.B[0]=0;
	fluxLW.S[0] += Riemann_vec.gr.alpha*Riemann_vec.gr.gdet*cons_LW.p;
#ifdef DIVBCLEANING_GR
	fluxLW.phi = cons_LW.D * v_coord_LW[0] * cons_LW.phi;
#endif

	//OUTPUT

	if(stamp==1) printf("v frame before %g \n", v_frame_line);
#ifdef HYDRO_MESHLESS_FINITE_MASS
	Riemann_out->Fluxes.rho = 0;
	v_frame_line = (omega_g*fluxLW.D + (1-omega_g)*fluxLF.D)/cons_LW.D;
	if(stamp==1) printf("v frame after %g \n", v_frame_line);
	Riemann_out->S_M=v_frame_line;
#else
	Riemann_out->Fluxes.rho = omega_g*fluxLW.D + (1-omega_g)*fluxLF.D - v_frame_line*cons_LW.D;
#ifdef DIVBCLEANING_GR
	Riemann_out->Fluxes.phi = omega_g*fluxLW.phi + (1-omega_g)*fluxLF.phi - v_frame_line*cons_LW.phi;
#endif
#endif
	Riemann_out->Fluxes.p = omega_g*fluxLW.tau + (1-omega_g)*fluxLF.tau - v_frame_line*cons_LW.tau;
	for(k=0;k<3;k++){
		Riemann_out->Fluxes.v[k] = omega_g*fluxLW.S[k] + (1-omega_g)*fluxLF.S[k] - v_frame_line*cons_LW.S[k];
		Riemann_out->Fluxes.B[k] = omega_g*fluxLW.B[k] + (1-omega_g)*fluxLF.B[k] - v_frame_line*cons_LW.B[k];
	}
	if(stamp==1)printf("Output Fluxes D %g - p %g - v %g %g %g - B %g %g %g \n", Riemann_out->Fluxes.rho,Riemann_out->Fluxes.p,Riemann_out->Fluxes.v[0],Riemann_out->Fluxes.v[1],Riemann_out->Fluxes.v[2],Riemann_out->Fluxes.B[0],Riemann_out->Fluxes.B[1],Riemann_out->Fluxes.B[2]);
}



void GR_MHD_Riemann_solver_HLLC(struct Input_vec_Riemann Riemann_vec,struct Riemann_outputs *Riemann_out, int gr_limiter)
{
	int k;
	//First rotate beta along n and the metric tensor
	double beta_n=Riemann_vec.gr.beta[0];
	double gnn=Riemann_vec.gr.ugamma[0];

	double v_frame_line = Riemann_vec.gr.alpha*Riemann_vec.v_frame[0] - beta_n;

	//Get conservative

	struct gr_conservatives cons_L,cons_R;
	struct gr_primitives prims;
	double W_l,W_r,vL[3],vR[3],tmpvec[3];
	double vAL,vAR;
	double bL[4],bR[4];
	double lam_l[2];
	double lam_r[2];
	double S_L,S_R;
	double beta_cov[3];
	raise_lower_indices(Riemann_vec.gr.gamma,Riemann_vec.gr.beta,beta_cov);


	//Bx continuity
	double BxL,BxR;
	BxL = Riemann_vec.L.B[0];
	BxR = Riemann_vec.R.B[0];
	double B_normal_corrected = 0.5*(BxL + BxR);
#ifdef DIVBCLEANING_GR
	/* (B_perp, psi) subsystem, Route B (gr/ROUTE_B_cleaning_notes.pdf, eqs. 13-14): Rusanov states with characteristic speed alpha*c_h*sqrt(gamma^nn).
	   L.phi/R.phi carry Dedner's psi = c_h*psihat of each side (multiplied in hydro_core_meshless.h); c_h is the face value (max over the pair).
	   Powers of c_h: 1/c_h on the psi jump in B*, c_h on the B jump in psi*. For c_h = 1 this is Fedrigo & Lupi (2025) eq. (32). */
	double c_h = Riemann_vec.c_h; /* (c_h^2 is not needed in the HLLC variant) */
	double phiL = Riemann_vec.L.phi, phiR = Riemann_vec.R.phi;
	double corr_norm = 1.0, corr_p = 0.;
	if(c_h > 0) {corr_p = 0.5*(phiL - phiR)*sqrt(gnn)/c_h;}
	double corr_p_abs = fabs(corr_p);
	double corr_b_abs = DEDNER_IMPLICIT_LIMITER_GR * fabs(B_normal_corrected);
	if(corr_p_abs > corr_b_abs) {corr_norm *= corr_b_abs/corr_p_abs; if(stamp==1){ printf("Correzione nel Riemann Solver AVVENUTA, fattore corr_norm %g \n", corr_norm);}}
	B_normal_corrected += corr_norm * corr_p;
	Riemann_out->phi_normal_mean = 0.5 * corr_norm * (phiL + phiR);
	Riemann_out->phi_normal_db = corr_norm * 0.5*(Riemann_vec.L.B[0] - Riemann_vec.R.B[0]) * c_h / sqrt(gnn);
#endif
	BxL = B_normal_corrected;
	BxR = B_normal_corrected;
	Riemann_out->B_normal_corrected = Riemann_vec.gr.gdet*B_normal_corrected;

   	//Left state
    	prims.rho0=Riemann_vec.L.rho;
    	for(k=0;k<3;k++) prims.v[k]=vL[k]=Riemann_vec.L.v[k];
    	prims.u=Riemann_vec.L.u;
    	prims.p=Riemann_vec.L.p;
    	W_l=lorentz_factor(Riemann_vec.gr.gamma,prims.v, gr_limiter);

	//Now get Rest Frame b field
    	prims.bsq=0;
    	raise_lower_indices(Riemann_vec.gr.gamma,prims.v,vL); //Get covariant
    	for(k=1;k<3;k++) {prims.bsq += Riemann_vec.L.B[k]*vL[k]; prims.B[k]=Riemann_vec.L.B[k];}
    	prims.bsq += BxL*vL[0];
    	prims.B[0] = BxL;

    	prims.b[0] = W_l*prims.bsq;
    	for(k=0;k<3;k++) tmpvec[k] = prims.B[k]/W_l + W_l*prims.bsq*(prims.v[k]-(Riemann_vec.gr.beta[k]/Riemann_vec.gr.alpha));
    	//Now lower indices for j components
    	raise_lower_indices(Riemann_vec.gr.gamma,tmpvec,prims.b+1);
    	for(k=1;k<=3; k++) {prims.b[k] += beta_cov[k-1]/Riemann_vec.gr.alpha*prims.b[0];}
    	prims.bsq = prims.bsq*prims.bsq + vector_sqnorm(Riemann_vec.gr.gamma,prims.B)/W_l/W_l;

    	vAL = sqrt(prims.bsq/(prims.rho0*(1+prims.u)+(prims.p+prims.bsq)));

    	for(k=0;k<4;k++) bL[k]=prims.b[k];
    	gr_get_conservative_variables(&prims,&cons_L,W_l,Riemann_vec.gr.gamma,Riemann_vec.gr.gdet);
	//Let's compute everything with the totol energy E. We'll then subtract the D directly in the fluxes.
	cons_L.tau += cons_L.D;

    	//Right state
    	prims.rho0=Riemann_vec.R.rho;
    	for(k=0;k<3;k++) prims.v[k]=vR[k]=Riemann_vec.R.v[k];
    	prims.u=Riemann_vec.R.u;
    	prims.p=Riemann_vec.R.p;
    	W_r=lorentz_factor(Riemann_vec.gr.gamma,prims.v, gr_limiter);

    	//Now get Rest Frame b field
    	prims.bsq=0;
    	raise_lower_indices(Riemann_vec.gr.gamma,prims.v,vR); //Get covariant
    	for(k=1;k<3;k++) {prims.bsq += Riemann_vec.R.B[k]*vR[k]; prims.B[k]=Riemann_vec.R.B[k];}
    	prims.bsq += BxR*vR[0];
    	prims.B[0] = BxR;

    	prims.b[0] = W_r*prims.bsq;
    	for(k=0;k<3;k++) tmpvec[k] = prims.B[k]/W_r + W_r*prims.bsq*(prims.v[k]-(Riemann_vec.gr.beta[k]/Riemann_vec.gr.alpha));
    	//Now lower indices for j components
    	raise_lower_indices(Riemann_vec.gr.gamma,tmpvec,prims.b+1);
    	for(k=1;k<=3; k++) {prims.b[k] += beta_cov[k-1]/Riemann_vec.gr.alpha*prims.b[0];}
    	prims.bsq = prims.bsq*prims.bsq + vector_sqnorm(Riemann_vec.gr.gamma,prims.B)/W_r/W_r;

    	vAR = sqrt(prims.bsq/(prims.rho0*(1+prims.u)+(prims.p+prims.bsq)));

    	for(k=0;k<4;k++) bR[k]=prims.b[k];
    	gr_get_conservative_variables(&prims,&cons_R,W_r,Riemann_vec.gr.gamma,Riemann_vec.gr.gdet);
	//Let's compute everything with the totol energy E. We'll then subtract the D directly in the fluxes.
	cons_R.tau += cons_R.D;

    	for(k=0;k<3;k++) {vL[k] = Riemann_vec.gr.alpha*Riemann_vec.L.v[k]-Riemann_vec.gr.beta[k]; vR[k]=Riemann_vec.gr.alpha*Riemann_vec.R.v[k]-Riemann_vec.gr.beta[k];}

    	//Compute lambda_hll with the corrected Bx
    	lam_l[0]=lam_r[0]=lam_l[1]=lam_r[1]=S_L=S_R=0;

    	compute_lambda_hll(Riemann_vec.L.v,Riemann_vec.L.cs,vAL,gnn,beta_n,Riemann_vec.gr.alpha,W_l,lam_l, gr_limiter);
    	compute_lambda_hll(Riemann_vec.R.v,Riemann_vec.R.cs,vAR,gnn,beta_n,Riemann_vec.gr.alpha,W_r,lam_r, gr_limiter);

    	S_L = DMIN(lam_r[1],lam_l[1]); //Lambda minus
    	S_R = DMAX(lam_r[0],lam_l[0]); //Lambda plus

	Riemann_out->S_M = 0.5*(S_L+S_R);
	Riemann_out->P_M = 0.5*(cons_L.p+cons_R.p);

    	struct gr_conservatives fluxL,fluxR;

	//Left flux
	fluxL.D    = cons_L.D  *vL[0];
    	fluxL.tau  = cons_L.tau*vL[0] + Riemann_vec.gr.gdet*cons_L.p*(vL[0]+beta_n) - Riemann_vec.gr.alpha*bL[0]*cons_L.B[0]/W_l;

    	for(k=0;k<3;k++){
        	fluxL.S[k] = cons_L.S[k]*vL[0] - Riemann_vec.gr.alpha*bL[k+1]*cons_L.B[0]/W_l;
		fluxL.B[k] = cons_L.B[k]*vL[0] - cons_L.B[0]*vL[k];
	}
	//Null Bx Flux
	fluxL.B[0] = 0;
	fluxL.S[0] += Riemann_vec.gr.alpha*Riemann_vec.gr.gdet*cons_L.p;

	//Right flux
	fluxR.D    = cons_R.D  *vR[0];
 	fluxR.tau  = cons_R.tau*vR[0] + Riemann_vec.gr.gdet*cons_R.p*(vR[0]+beta_n) - Riemann_vec.gr.alpha*bR[0]*cons_R.B[0]/W_r;

    	for(k=0;k<3;k++){
	        fluxR.S[k] = cons_R.S[k]*vR[0] - Riemann_vec.gr.alpha*bR[k+1]*cons_R.B[0]/W_r;
		fluxR.B[k] = cons_R.B[k]*vR[0] - cons_R.B[0]*vR[k];
	}
	//Null Bx Flux
	fluxR.B[0]=0;
	fluxR.S[0] += Riemann_vec.gr.alpha*Riemann_vec.gr.gdet*cons_R.p;

	//OUTPUT

#ifdef HYDRO_MESHLESS_FINITE_VOLUME
	if(S_L>v_frame_line)
	{
	  Riemann_out->Fluxes.rho = fluxL.D - cons_L.D*v_frame_line;
	  Riemann_out->Fluxes.p = fluxL.tau - cons_L.tau*v_frame_line;
	  //Fluxes_tau = Fluxes_E - Fluxes_D
	  Riemann_out->Fluxes.p -= Riemann_out->Fluxes.rho;
	  for(k=0;k<3;k++)
	  {
		Riemann_out->Fluxes.v[k] = fluxL.S[k] - cons_L.S[k]*v_frame_line;
		Riemann_out->Fluxes.B[k] = fluxL.B[k] - cons_L.B[k]*v_frame_line;
	  }
	}
	else if(S_R<v_frame_line)
	{
	  Riemann_out->Fluxes.rho = fluxR.D - cons_R.D*v_frame_line;
	  Riemann_out->Fluxes.p = fluxR.tau - cons_R.tau*v_frame_line;
	  //Fluxes_tau = Fluxes_E - Fluxes_D
	  Riemann_out->Fluxes.p -= Riemann_out->Fluxes.rho;
	  for(k=0;k<3;k++)
	  {
		Riemann_out->Fluxes.v[k] = fluxR.S[k] - cons_R.S[k]*v_frame_line;
		Riemann_out->Fluxes.B[k] = fluxR.B[k] - cons_R.B[k]*v_frame_line;
	  }
	}
    	else if(S_L<=v_frame_line && v_frame_line<=S_R)
	{
#endif
      	  //Get HLL intermediate state and fluxes
          struct gr_conservatives U_hll;
      	  double facLR=1./(S_R-S_L);
      	  U_hll.D = (S_R*cons_R.D-S_L*cons_L.D + fluxL.D - fluxR.D)*facLR;
      	  U_hll.tau = (S_R*cons_R.tau-S_L*cons_L.tau + fluxL.tau - fluxR.tau)*facLR;
      	  double Bx_hll;
      	  for(k=0;k<3;k++)
      	  {
        	U_hll.S[k] = (S_R*cons_R.S[k]-S_L*cons_L.S[k] + fluxL.S[k] - fluxR.S[k])*facLR;
	 	U_hll.B[k] = (S_R*cons_R.B[k]-S_L*cons_L.B[k] + fluxL.B[k] - fluxR.B[k])*facLR;
	  }

	  //Now compute the star R and L states
	  struct gr_conservatives cons_star_L,cons_star_R;
	  double B_star[3], v_star[3], p_star;

	  //By and Bz are easy (L=R)
      	  B_star[0] = BxL;
	  B_star[1] = U_hll.B[1];
	  B_star[2] = U_hll.B[2];

	  //Compute vx_star (L=R)
	  double a, b, c, Flux_Sx_HLL, Flux_E_HLL, Flux_By_HLL, Flux_Bz_HLL, vB_star, denomL, denomR, gamma_sqr_inv;
	  Flux_Sx_HLL = (S_R*fluxL.S[0]-S_L*fluxR.S[0] + S_L*S_R*(cons_R.S[0] - cons_L.S[0]))*facLR;
	  Flux_E_HLL = (S_R*fluxL.tau-S_L*fluxR.tau + S_L*S_R*(cons_R.tau - cons_L.tau))*facLR;
	  if(fabs(B_star[0]) < BX_MIN){
	    Flux_By_HLL = Flux_Bz_HLL = 0;
	    a = Flux_E_HLL;
	    b = -Flux_Sx_HLL - U_hll.tau;
	    c = U_hll.S[0];
	  }else{
	    Flux_By_HLL = (S_R*fluxL.B[1]-S_L*fluxR.B[1] + S_L*S_R*(cons_R.B[1] - cons_L.B[1]))*facLR;
	    Flux_Bz_HLL = (S_R*fluxL.B[2]-S_L*fluxR.B[2] + S_L*S_R*(cons_R.B[2] - cons_L.B[2]))*facLR;

	    a = Flux_E_HLL - (U_hll.B[1]*Flux_By_HLL + U_hll.B[2]*Flux_Bz_HLL);
	    b = - Flux_Sx_HLL - U_hll.tau + (U_hll.B[1]*U_hll.B[1] + U_hll.B[2]*U_hll.B[2]) + (Flux_By_HLL*Flux_By_HLL + Flux_Bz_HLL*Flux_Bz_HLL);
	    c = U_hll.S[0] - (U_hll.B[1]*Flux_By_HLL + U_hll.B[2]*Flux_Bz_HLL);
	  } //BX_MIN check

      	  v_star[0] = (- b - sqrt(b*b - 4*a*c))/2/a;
      	  denomL = 1/(S_L - v_star[0]);
      	  denomR = 1/(S_R - v_star[0]);
      	  //compute D. It's the same whether Bx==0 or Bx=!0
      	  cons_star_L.D = cons_L.D*(S_L - Riemann_vec.L.v[0])*denomL;
      	  cons_star_R.D = cons_R.D*(S_R - Riemann_vec.R.v[0])*denomR;

	  if(fabs(B_star[0]) < BX_MIN){
	    //compute p, Sy, Sz and tau. Finally Sx
	    p_star = Flux_Sx_HLL - v_star[0]*Flux_E_HLL;
     	    cons_star_L.S[1] = cons_L.S[1]*(S_L - Riemann_vec.L.v[0])*denomL;
      	    cons_star_L.S[2] = cons_L.S[2]*(S_L - Riemann_vec.L.v[0])*denomL;
     	    cons_star_R.S[1] = cons_R.S[1]*(S_R - Riemann_vec.R.v[0])*denomR;
      	    cons_star_R.S[2] = cons_R.S[2]*(S_R - Riemann_vec.R.v[0])*denomR;
	    cons_star_L.tau = (S_L*cons_L.tau - cons_L.S[0] + p_star*v_star[0])*denomL;
	    cons_star_R.tau = (S_R*cons_R.tau - cons_R.S[0] + p_star*v_star[0])*denomR;
	    cons_star_L.B[1] = cons_L.B[1]*(S_L - Riemann_vec.L.v[0])*denomL;
	    cons_star_L.B[2] = cons_L.B[2]*(S_L - Riemann_vec.L.v[0])*denomL;
	    cons_star_R.B[1] = cons_R.B[1]*(S_R - Riemann_vec.R.v[0])*denomR;
	    cons_star_R.B[2] = cons_R.B[2]*(S_R - Riemann_vec.R.v[0])*denomR;
	    cons_star_L.S[0] = (cons_star_L.tau + p_star)*v_star[0];
	    cons_star_R.S[0] = (cons_star_R.tau + p_star)*v_star[0];
	    cons_star_L.B[0] = cons_star_R.B[0] = 0;
	    v_star[1] = v_star[2] = 0;
	    vB_star = gamma_sqr_inv = 0;
	  }else{
	    //By and Bz
	    for(k=0;k<3;k++){cons_star_L.B[k] = B_star[k]; cons_star_R.B[k] = B_star[k];}
      	    //compute vy and vz (L=R)
      	    v_star[1] = (B_star[1]*v_star[0] - Flux_By_HLL)/B_star[0];
      	    v_star[2] = (B_star[2]*v_star[0] - Flux_Bz_HLL)/B_star[0];
      	    //compute P (L=R)
      	    gamma_sqr_inv = 1 - v_star[0]*v_star[0] - v_star[1]*v_star[1] - v_star[2]*v_star[2];
      	    vB_star = v_star[0]*B_star[0] + v_star[1]*B_star[1] + v_star[2]*B_star[2];
      	    p_star = Flux_Sx_HLL + (B_star[0]*B_star[0])*gamma_sqr_inv - v_star[0]*(Flux_E_HLL - B_star[0]*vB_star);

      	    //compute tau
      	    cons_star_L.tau = (S_L*cons_L.tau - cons_L.S[0] + p_star*v_star[0] - vB_star*B_star[0])*denomL;
      	    cons_star_R.tau = (S_R*cons_R.tau - cons_R.S[0] + p_star*v_star[0] - vB_star*B_star[0])*denomR;

      	    //compute sx, sy and sz
      	    cons_star_L.S[0] = (cons_star_L.tau + p_star)*v_star[0] - vB_star*B_star[0];
     	    cons_star_L.S[1] = (-B_star[0]*(B_star[1]*gamma_sqr_inv + vB_star*v_star[1]) + S_L*cons_L.S[1] - fluxL.S[1])*denomL;
      	    cons_star_L.S[2] = (-B_star[0]*(B_star[2]*gamma_sqr_inv + vB_star*v_star[2]) + S_L*cons_L.S[2] - fluxL.S[2])*denomL;
      	    cons_star_R.S[0] = (cons_star_R.tau + p_star)*v_star[0] - vB_star*B_star[0];
      	    cons_star_R.S[1] = (-B_star[0]*(B_star[1]*gamma_sqr_inv + vB_star*v_star[1]) + S_R*cons_R.S[1] - fluxR.S[1])*denomR;
      	    cons_star_R.S[2] = (-B_star[0]*(B_star[2]*gamma_sqr_inv + vB_star*v_star[2]) + S_R*cons_R.S[2] - fluxR.S[2])*denomR;
	  } //BX_MIN check

      	  //compute the star fluxes
      	  struct gr_conservatives flux_star_L,flux_star_R;

	  //Left star flux
	  flux_star_L.D    = cons_star_L.D  *v_star[0];
      	  flux_star_L.tau  = cons_star_L.S[0];
      	  for(k=0;k<3;k++){
          	flux_star_L.S[k] = cons_star_L.S[k]*v_star[0] - cons_star_L.B[0]*cons_star_L.B[k]*gamma_sqr_inv - cons_star_L.B[0]*v_star[k]*vB_star;
		flux_star_L.B[k] = cons_star_L.B[k]*v_star[0] - cons_star_L.B[0]*v_star[k];
	  }
	  //Null Bx Flux
	  flux_star_L.B[0] = 0;
	  flux_star_L.S[0] += p_star;

	  //Right star flux
	  flux_star_R.D    = cons_star_R.D  *v_star[0];
	  flux_star_R.tau  = cons_star_R.S[0];
      	  for(k=0;k<3;k++){
          	flux_star_R.S[k] = cons_star_R.S[k]*v_star[0] - cons_star_R.B[0]*cons_star_R.B[k]*gamma_sqr_inv - cons_star_R.B[0]*v_star[k]*vB_star;
		flux_star_R.B[k] = cons_star_R.B[k]*v_star[0] - cons_star_R.B[0]*v_star[k];
	  }
	  //Null Bx Flux
	  flux_star_R.B[0] = 0;
	  flux_star_R.S[0] += p_star;

      	  double S_star;
      	  S_star = v_star[0];

      	  //middle states outputs
      	  if(S_L <= v_frame_line && v_frame_line<S_star)
	  {
		Riemann_out->Fluxes.rho = flux_star_L.D - v_frame_line*cons_star_L.D;
        	Riemann_out->Fluxes.p = flux_star_L.tau - v_frame_line*cons_star_L.tau;
	  	//Fluxes_tau = Fluxes_E - Fluxes_D
	  	Riemann_out->Fluxes.p -= Riemann_out->Fluxes.rho;
		for(k=0;k<3;k++)
		{
	 		Riemann_out->Fluxes.v[k] = flux_star_L.S[k] - v_frame_line*cons_star_L.S[k];
	 		Riemann_out->Fluxes.B[k] = flux_star_L.B[k] - v_frame_line*cons_star_L.B[k];
		}
      	  }
      	  else if(S_star <= v_frame_line && v_frame_line<=S_R)
	  {
        	Riemann_out->Fluxes.rho = flux_star_R.D - v_frame_line*cons_star_R.D;
        	Riemann_out->Fluxes.p = flux_star_R.tau - v_frame_line*cons_star_R.tau;
		//Fluxes_tau = Fluxes_E - Fluxes_D
	  	Riemann_out->Fluxes.p -= Riemann_out->Fluxes.rho;
		for(k=0;k<3;k++)
		{
	 		Riemann_out->Fluxes.v[k] = flux_star_R.S[k] - v_frame_line*cons_star_R.S[k];
	 		Riemann_out->Fluxes.B[k] = flux_star_R.B[k] - v_frame_line*cons_star_R.B[k];
		}
      	  }
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    	}
#endif
#ifdef HYDRO_MESHLESS_FINITE_MASS
	Riemann_out->Fluxes.rho = 0;
	Riemann_out->Fluxes.p = p_star*S_star - B_star[0]*vB_star;
	Riemann_out->Fluxes.v[0] = p_star - B_star[0]*B_star[0]*gamma_sqr_inv - B_star[0]*v_star[0]*vB_star;
	Riemann_out->Fluxes.v[1] = - B_star[0]*B_star[1]*gamma_sqr_inv - B_star[0]*v_star[1]*vB_star;
	Riemann_out->Fluxes.v[2] = - B_star[0]*B_star[2]*gamma_sqr_inv - B_star[0]*v_star[2]*vB_star;
	Riemann_out->Fluxes.B[0] = 0;
	Riemann_out->Fluxes.B[1] = -cons_star_L.B[0]*v_star[1];
	Riemann_out->Fluxes.B[2] = -cons_star_L.B[0]*v_star[2];
#endif

}


#endif
