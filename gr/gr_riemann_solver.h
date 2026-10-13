#ifndef MAGNETIC
void compute_lambda_hll(double v[],double cs,double g_nn,double beta_n, double alpha,double n_unit[],double Wfac,double lam[2],int gr_limiter)
{
	//Get u_n (Euclidean projection)
	double vtilde_n = v[0]*n_unit[0]+v[1]*n_unit[1]+v[2]*n_unit[2];
	double cs2=cs*cs;
	double v2max = 1.-MIN_REAL_NUMBER;
	if(gr_limiter) v2max = 1.-1./GR_GAMMA_LIMITER/GR_GAMMA_LIMITER;
	long double v2=1.-1./Wfac/Wfac;
	if(v2<0) v2=0.;


	long double discr=sqrtl((1-v2)*(g_nn*(1-v2*cs2)-vtilde_n*vtilde_n*(1-cs2)));
	lam[0] = alpha*DMIN((vtilde_n*(1-cs2)+cs*discr)/(1-v2*cs2),v2max) -beta_n;
	lam[1] = alpha*DMAX((vtilde_n*(1-cs2)-cs*discr)/(1-v2*cs2),-v2max) -beta_n;

}

void compute_starstate(struct gr_conservatives cons_L,struct gr_conservatives cons_R,double *n_unit,double vline_L,double vline_R,
	struct gr_metric *gr,double beta_n,double g_nn,double S_L,double S_R,double *S_M,double *P_M,int gr_limiter)
{

    double sXL, sXR;
    double vtmp[3];
    // the coefficients of the quadratic equaion determining lambdastar
    long double AL,AR,BL,BR,CL,CR,EL,ER;

    // the combination D*h*W which appear more than once
    //double DhWL,DhWR;
    long double lamstar_tilde; // eq (22)

    //Get conserved variables
    //DhWR = cons_R.D+cons_R.tau+ gdet*cons_R.p;
    //DhWL = cons_L.D+cons_L.tau+ gdet*cons_L.p;
    raise_lower_indices(gr->ugamma,cons_R.S,vtmp);
    sXR = (vtmp[0]*n_unit[0]+vtmp[1]*n_unit[1]+vtmp[2]*n_unit[2]);// -beta_n/gr->alpha*DhWR;
    raise_lower_indices(gr->ugamma,cons_L.S,vtmp);
    sXL = (vtmp[0]*n_unit[0]+vtmp[1]*n_unit[1]+vtmp[2]*n_unit[2]);// -beta_n/gr->alpha*DhWL;
   
    // compute the coefficients A-D of Duez (23)-(26)
    AR = (cons_R.tau+cons_R.D)*(S_R+beta_n) - (gr->alpha*sXR);//+DhWR*beta_n);
    AL = (cons_L.tau+cons_L.D)*(S_L+beta_n) - (gr->alpha*sXL);//+DhWL*beta_n);
    CR = gr->alpha*g_nn;
    CL = CR;
    //BR = CR*gdet*cons_R.p - (sXR+DhWR*beta_n/alpha)*(S_R-vline_R);
    //BL = CL*gdet*cons_L.p - (sXL+DhWL*beta_n/alpha)*(S_L-vline_L);
    BR = CR*gr->gdet*cons_R.p + (sXR)*(vline_R-S_R);
    BL = CL*gr->gdet*cons_L.p + (sXL)*(vline_L-S_L);
    ER = -(beta_n + S_R);
    EL = -(beta_n + S_L);

    // we require the negative root only (see Mignone et al), this is (27) of
    long double a=(EL*AR-AL*ER);
    long double b=-(AL*CR+BL*ER-CL*AR-EL*BR);
    long double c=BR*CL-BL*CR;
    if(a==0)
    {
	if(b==0){ printf("Error in star state; Unable to find a solution: b is 0 %2.12Lg %2.12Lg %2.12Lg\n",a,b,c); endrun(9991010);}
	lamstar_tilde = -c/b;
    }
    else{
	//we only want the negative solution
	long double discr = b * b - 4 * a * c;
    	if (discr < 0) { printf("Error in star state; Unable to find a solution: Delta<0 %2.12Lg %2.12Lg %2.12Lg - discr %2.12Lg - S %2.12g %2.12g\n",a,b,c,discr,S_L,S_R);
			 printf("Cons L D %g - S %g %g %g - tau %g - Cons R D %g - S %g %g %g - tau %g \n", cons_L.D, cons_L.S[0], cons_L.S[1], cons_L.S[2], cons_L.tau,cons_R.D, cons_R.S[0], cons_R.S[1], cons_R.S[2], cons_R.tau); endrun(9991011);}
    	else if(discr==0) {lamstar_tilde=-0.5*b/a; printf("Warning: discr=0; b=%Lg a=%Lg SM %Lg\n",b,a,lamstar_tilde);}
	else
    	{
		long double tmp;
		if(b>0) tmp= -0.5*(b + sqrtl(discr)); else {tmp = -0.5*(b - sqrtl(discr));}
		lamstar_tilde = c/tmp;
		if(fabs(lamstar_tilde)>1){
		  printf("My other speed is %2.12Lg\n",tmp/a);
		  printf("My states were %Lg %Lg %Lg\n",a,b,c);
		  printf("L: %g %g %g - R %g %g %g\n",cons_L.D,cons_L.S[0],cons_L.tau,cons_R.D,cons_R.S[0],cons_R.tau);

		  if(gr_limiter)  lamstar_tilde=sqrt(1.-1./GR_GAMMA_LIMITER/GR_GAMMA_LIMITER);
		  else endrun(9991012);
		}
    	}
    }
    // pressstar is equal left and right of contact discontinuity
    *P_M = (AL*(lamstar_tilde) + BL)/((CL+EL*(lamstar_tilde))*gr->gdet);


    //Update S_M
    *S_M = gr->alpha * (lamstar_tilde) - beta_n; // eq (22)
    if(*S_M>S_R || *S_M<S_L)
    {
#if (GENERAL_RELATIVITY==1)
	   printf("Warning! SM outside range. %g %g %g - PM %g lam %Lg\n",S_L,*S_M,S_R,*P_M,lamstar_tilde);
	   printf("My values are AL %2.12Lg AR %2.12Lg, BR+AR %2.12Lg BL-AL %2.12Lg, SL*AR-AR %2.12Lg SR*AL-BL %2.12Lg\n",AL,AR,-BR+AR,-BL-AL,S_L*AR-AR,S_R*AL+BL);
	   printf("Extra: SL*AL-BL: %2.12Lg SL*AR-BR: %2.12Lg SR*AL-BL: %2.12Lg SR*AR-BR: %2.12Lg\n",
		S_L*AL+BL,S_L*AR+BR,S_R*AL+BL,S_R*AR+BR);
	   printf("limits from MB05 %g %g\n",fmax(-BR/AR,-BL/AL),fmin(-BR/AR,-BL/AL));
	   printf("My vels are %g %g\n",vline_L,vline_R);
	   printf("My coeffs %Lg %Lg %Lg\n",a,b,c);
	   printf("Unit: %g %g %g\n",n_unit[0],n_unit[1],n_unit[2]);
#endif
    }
}

void get_HLL_flux(struct gr_conservatives cons_L,struct gr_conservatives cons_R,double vline_L,double vline_R,double S_L,double S_R,
		double v_frame_line,double alpha,double beta_n, double gdet,double *n_unit,struct Riemann_outputs *Riemann_out)
{
	int k;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
	 if(S_L>v_frame_line)
	 {
	   Riemann_out->Fluxes.rho = cons_L.D  *(vline_L - v_frame_line);
           Riemann_out->Fluxes.p   = cons_L.tau*(vline_L - v_frame_line) + gdet*cons_L.p*(vline_L+beta_n);
           for(k=0;k<3;k++)
                Riemann_out->Fluxes.v[k] = cons_L.S[k]*(vline_L - v_frame_line) + alpha*gdet*cons_L.p*n_unit[k];

	 }
	 else if(S_R<v_frame_line)
	 {
	   Riemann_out->Fluxes.rho = cons_R.D  *(vline_R - v_frame_line);
           Riemann_out->Fluxes.p   = cons_R.tau*(vline_R - v_frame_line) + gdet*cons_R.p*(vline_R+beta_n);
           for(k=0;k<3;k++)
                Riemann_out->Fluxes.v[k] = cons_R.S[k]*(vline_R - v_frame_line) + alpha*gdet*cons_R.p*n_unit[k];

	 }
	 else if(S_L<=v_frame_line && v_frame_line<=S_R)
#endif
	 {
	   //Get HLL intermediate state and fluxes
	   struct gr_conservatives U_hll;
           double facLR=1./(S_R-S_L);
	   double flux_hll;
	   U_hll.D = (S_R*cons_R.D-S_L*cons_L.D + cons_L.D*vline_L - cons_R.D*vline_R)*facLR;
	   U_hll.tau = (S_R*cons_R.tau-S_L*cons_L.tau + cons_L.tau*vline_L - cons_R.tau*vline_R +
 		gdet*(cons_L.p*(vline_L+beta_n)-cons_R.p*(vline_R+beta_n)))*facLR;
	   for(k=0;k<3;k++)
		U_hll.S[k] = (S_R*cons_R.S[k]-S_L*cons_L.S[k] + cons_L.S[k]*vline_L - cons_R.S[k]*vline_R +
				alpha*gdet*n_unit[k]*(cons_L.p-cons_R.p))*facLR;

#ifdef HYDRO_MESHLESS_FINITE_MASS
	   //This is the HLL equivalent of the zero-mass-flux face velocity
	   //Fluxrho = (SR*FL - SL*FR + SL*SR*(UR-UL))/(SR-SL)
	   //Uhll    = (SR*UR - SL*UL + FL - FR)/(SR-SL)
	   if(stamp==1) printf("v frame line %g\n ", v_frame_line);
	   v_frame_line = (S_R*cons_L.D*vline_L-S_L*cons_R.D*vline_R + S_L*S_R*(cons_R.D - cons_L.D))*facLR/U_hll.D;
	   if(stamp==1) printf("S_R %g - S_L %g - vLx %g - vRx %g - D_L %g - D_R %g \n V frame line %g\n", S_R, S_L, vline_L,vline_R , cons_L.D, cons_R.D, v_frame_line);
	   Riemann_out->Fluxes.rho = 0;
#else
	   Riemann_out->Fluxes.rho = (S_R*cons_L.D*vline_L-S_L*cons_R.D*vline_R + S_L*S_R*(cons_R.D - cons_L.D))*facLR - v_frame_line*U_hll.D;
#endif

	   Riemann_out->Fluxes.p = facLR *(S_R*(cons_L.tau*(vline_L-S_L)+gdet*cons_L.p*(vline_L+beta_n)) -
                                          S_L*(cons_R.tau*(vline_R-S_R)+gdet*cons_R.p*(vline_R+beta_n)) )-v_frame_line*U_hll.tau;

	   for(k=0;k<3;k++)
	   {
		Riemann_out->Fluxes.v[k] = facLR*((S_R)*(cons_L.S[k]*(vline_L-S_L)+alpha*gdet*cons_L.p*n_unit[k]) -
                                                  (S_L)*(cons_R.S[k]*(vline_R-S_R)+alpha*gdet*cons_R.p*n_unit[k]))
				-v_frame_line*U_hll.S[k];
	   }
	 }

}

void get_HLLC_flux(struct gr_conservatives cons_L,struct gr_conservatives cons_R,double vline_L,double vline_R,double S_L,double S_R,double v_frame_line,
		double alpha,double beta_n, double gdet,double *n_unit,struct Riemann_outputs *Riemann_out)
{
	int k;
	double S_M=Riemann_out->S_M;
	double P_M=Riemann_out->P_M;
	struct gr_conservatives U_state;

	if((S_M==S_L)||(S_M==S_R))
    	{
     	 /* trap for this case, which gives NAN below but is actually very simple */
         Riemann_out->Fluxes.rho = 0;
         Riemann_out->Fluxes.p = gdet * P_M * S_M;
         for(k=0;k<3;k++)
            Riemann_out->Fluxes.v[k] = alpha * gdet * P_M * n_unit[k];
         return;
    	}



	if(S_L>v_frame_line)
	{
	 Riemann_out->Fluxes.rho = cons_L.D  *(vline_L - v_frame_line);
	 Riemann_out->Fluxes.p   = cons_L.tau*(vline_L - v_frame_line) + gdet*cons_L.p*(vline_L+beta_n);
	 for(k=0;k<3;k++)
	 	Riemann_out->Fluxes.v[k] = cons_L.S[k]*(vline_L - v_frame_line) + alpha*gdet*cons_L.p*n_unit[k];
	}
	else if(S_L<=v_frame_line && v_frame_line<=S_M)
	{
	 double v_fac = (S_L-vline_L)/(S_L-S_M);
	 if(v_fac<0) v_fac=0;
	 U_state.D    = v_fac*cons_L.D;
	 U_state.S[0] = v_fac*cons_L.S[0] + alpha*gdet*(cons_L.p-P_M) / (S_M-S_L) * n_unit[0];
	 U_state.S[1] = v_fac*cons_L.S[1] + alpha*gdet*(cons_L.p-P_M) / (S_M-S_L) * n_unit[1];
	 U_state.S[2] = v_fac*cons_L.S[2] + alpha*gdet*(cons_L.p-P_M) / (S_M-S_L) * n_unit[2];
	 U_state.tau  = v_fac*cons_L.tau  + gdet*((vline_L+beta_n)*cons_L.p-(S_M+beta_n)*P_M)/(S_M-S_L);
	 U_state.p    = P_M;

         Riemann_out->Fluxes.rho = U_state.D  *(S_M - v_frame_line);
         Riemann_out->Fluxes.p   = U_state.tau*(S_M - v_frame_line) + gdet*P_M*(S_M+beta_n);
         for(k=0;k<3;k++)
                Riemann_out->Fluxes.v[k] = U_state.S[k]*(S_M - v_frame_line) + alpha*gdet*P_M*n_unit[k];

	}
	else if(S_M<=v_frame_line && v_frame_line<=S_R)
	{
	 //F_*R
	 double v_fac = (S_R-vline_R)/(S_R-S_M);
	 if(v_fac<0) v_fac=0;
	 U_state.D    = v_fac*cons_R.D;
	 U_state.S[0] = v_fac*cons_R.S[0] + alpha*gdet*(P_M-cons_R.p) / (S_R-S_M) * n_unit[0];
	 U_state.S[1] = v_fac*cons_R.S[1] + alpha*gdet*(P_M-cons_R.p) / (S_R-S_M) * n_unit[1];
	 U_state.S[2] = v_fac*cons_R.S[2] + alpha*gdet*(P_M-cons_R.p) / (S_R-S_M) * n_unit[2];
	 U_state.tau  = v_fac*cons_R.tau  + gdet*((S_M+beta_n)*P_M-(vline_R+beta_n)*cons_R.p)/(S_R-S_M);
	 U_state.p    = P_M;


         Riemann_out->Fluxes.rho = U_state.D  *(S_M - v_frame_line);
         Riemann_out->Fluxes.p   = U_state.tau*(S_M - v_frame_line) + gdet*P_M*(S_M+beta_n);
         for(k=0;k<3;k++)
                Riemann_out->Fluxes.v[k] = U_state.S[k]*(S_M - v_frame_line) + alpha*gdet*P_M*n_unit[k];
	}
	else if(S_R<v_frame_line)
	{
	 Riemann_out->Fluxes.rho = cons_R.D  *(vline_R - v_frame_line);
	 Riemann_out->Fluxes.p   = cons_R.tau*(vline_R - v_frame_line) + gdet*cons_R.p*(vline_R+beta_n);
	 for(k=0;k<3;k++)
	 	Riemann_out->Fluxes.v[k] = cons_R.S[k]*(vline_R - v_frame_line) + alpha*gdet*cons_R.p*n_unit[k];
	}
}

void GR_Riemann_solver(struct Input_vec_Riemann Riemann_vec,struct Riemann_outputs *Riemann_out,double* n_unit,int gr_limiter)
{

	//First rotate beta along n and the metric tensot
        double beta_n=Riemann_vec.gr.beta[0]*n_unit[0]+Riemann_vec.gr.beta[1]*n_unit[1]+Riemann_vec.gr.beta[2]*n_unit[2];

        double nxy=sqrt(n_unit[0]*n_unit[0]+n_unit[1]*n_unit[1]); //sin theta
        double grot[6];
	if(nxy>0)
        	rotate_tensor(Riemann_vec.gr.ugamma,n_unit[0]/nxy,-n_unit[1]/nxy,nxy,-n_unit[2],grot);
	else
        	rotate_tensor(Riemann_vec.gr.ugamma,1.,0.,nxy,-n_unit[2],grot);

	double gnn=grot[0]; //We align along x the face normal vector

	//Compute the line velocity
	double vline_L = Riemann_vec.gr.alpha*(Riemann_vec.L.v[0]*n_unit[0]+Riemann_vec.L.v[1]*n_unit[1]+Riemann_vec.L.v[2]*n_unit[2]) - beta_n;
	double vline_R = Riemann_vec.gr.alpha*(Riemann_vec.R.v[0]*n_unit[0]+Riemann_vec.R.v[1]*n_unit[1]+Riemann_vec.R.v[2]*n_unit[2]) - beta_n;

	double v_frame_line = Riemann_vec.gr.alpha*(Riemann_vec.v_frame[0]*n_unit[0]+Riemann_vec.v_frame[1]*n_unit[1]+
			Riemann_vec.v_frame[2]*n_unit[2]) - beta_n;

	if(stamp==1) printf("RIEMANN SOLVER: \n alpha %g - beta_n %g - Riemann_vec-v.L %g %g %g - Riemann_vec.v.R %g %g %g \n vline_L %g - vline_R %g \n", Riemann_vec.gr.alpha, beta_n, Riemann_vec.L.v[0], Riemann_vec.L.v[1], Riemann_vec.L.v[2], Riemann_vec.R.v[0], Riemann_vec.R.v[1], Riemann_vec.R.v[2], vline_L, vline_R);	

	//Get conserved
	struct gr_conservatives cons_L,cons_R;
	struct gr_primitives prims;
	double W_l,W_r;

	prims.rho0=Riemann_vec.L.rho;
	prims.v[0]=Riemann_vec.L.v[0];
	prims.v[1]=Riemann_vec.L.v[1];
	prims.v[2]=Riemann_vec.L.v[2];
	prims.u=Riemann_vec.L.u;
	prims.p=Riemann_vec.L.p;
	W_l=lorentz_factor(Riemann_vec.gr.gamma,prims.v,gr_limiter);
	gr_get_conservative_variables(&prims,&cons_L,W_l,Riemann_vec.gr.gamma,Riemann_vec.gr.gdet);

	prims.rho0=Riemann_vec.R.rho;
	prims.v[0]=Riemann_vec.R.v[0];
	prims.v[1]=Riemann_vec.R.v[1];
	prims.v[2]=Riemann_vec.R.v[2];
	prims.u=Riemann_vec.R.u;
	prims.p=Riemann_vec.R.p;
	W_r=lorentz_factor(Riemann_vec.gr.gamma,prims.v,gr_limiter);
	gr_get_conservative_variables(&prims,&cons_R,W_r,Riemann_vec.gr.gamma,Riemann_vec.gr.gdet);

	//Now get lambda_hll
	double lam_l[2];
	double lam_r[2];

	compute_lambda_hll(Riemann_vec.L.v,Riemann_vec.L.cs,gnn,beta_n,Riemann_vec.gr.alpha,n_unit,W_l,lam_l,gr_limiter);
	compute_lambda_hll(Riemann_vec.R.v,Riemann_vec.R.cs,gnn,beta_n,Riemann_vec.gr.alpha,n_unit,W_r,lam_r,gr_limiter);

        double S_L=DMIN(lam_r[1],lam_l[1]); //Lambda minus
        double S_R=DMAX(lam_r[0],lam_l[0]);//Lambda plus

	/*
	//Check if the two states are valid
	double E_L=sqrtl(vector_sqnorm(Riemann_vec.gr.ugamma,cons_L.S)+cons_L.D*cons_L.D);
	double E_R=sqrtl(vector_sqnorm(Riemann_vec.gr.ugamma,cons_R.S)+cons_R.D*cons_R.D);

	if(!(cons_L.D>0 && cons_L.tau+cons_L.D>E_L && cons_R.D>0 && cons_R.tau+cons_R.D>E_R))
	{
		printf("Invalid state entering the Riemann solver\n");
		endrun(57482);
	}
	*/

	compute_starstate(cons_L,cons_R,n_unit,vline_L,vline_R,&(Riemann_vec.gr),beta_n,gnn,S_L,S_R,&(Riemann_out->S_M),&(Riemann_out->P_M),gr_limiter);

	if(Riemann_out->P_M<0) return;

#if (GENERAL_RELATIVITY==1) //HLLC solver by MB05

#ifdef HYDRO_MESHLESS_FINITE_MASS
	Riemann_out->Fluxes.rho  = 0;
 	Riemann_out->Fluxes.v[0] = Riemann_out->P_M*n_unit[0];
 	Riemann_out->Fluxes.v[1] = Riemann_out->P_M*n_unit[1];
 	Riemann_out->Fluxes.v[2] = Riemann_out->P_M*n_unit[2];
 	Riemann_out->Fluxes.p    = Riemann_out->P_M*Riemann_out->S_M;
#else
	get_HLLC_flux(cons_L,cons_R,vline_L,vline_R,S_L,S_R,v_frame_line,Riemann_vec.gr.alpha,beta_n,Riemann_vec.gr.gdet,n_unit,Riemann_out);
#endif

#else //GENERAL_RELATIVITY==0 - HLL solver
	get_HLL_flux(cons_L,cons_R,vline_L,vline_R,S_L,S_R,v_frame_line,Riemann_vec.gr.alpha,beta_n,Riemann_vec.gr.gdet,n_unit,Riemann_out);
#endif //close GENERAL_RELATIVITY==1/0

/*-----------------------------------------*/

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

#endif
