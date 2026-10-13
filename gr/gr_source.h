#ifdef GENERAL_RELATIVITY

void get_source_term(int target,struct gr_primitives *prims,struct gr_conservatives *Src)
{
	int dir,j,k,w;
	double Tmunu[4][4];
	double Kij[3][3];
#ifdef MAGNETIC
	double v_cov[3];
	double vdotB;
	double b_con[4];	
#endif
	double d_alpha[3];
	double d_beta[3][3],d_gij[3][4][4];

	/* --- Compute the stress-energy tensor with first index up, second index down */

	//First, determine enthalpy and 4-velocity vector
	double Pstar = prims->p;
  	double rhoh=prims->rho0*(1+prims->u) + prims->p;
	int gr_limiter=0;
	double dpl[3]; for(k=0;k<3;k++) {dpl[k]=P[target].Pos[k];}
#if (GR_METRIC==2)
	for(k=0;k<3;k++) {dpl[k] = -P[target].Min_xyz_to_Sink[k];}
#endif
	if(dpl[0]*dpl[0]+dpl[1]*dpl[1]+dpl[2]*dpl[2] <= pow(GAMMA_LIMITER_RADIUS,2)) gr_limiter=1;
//	if(P[target].Pos[0]<LEFT_GAMMA_BOUNDARY||P[target].Pos[0]>RIGHT_GAMMA_BOUNDARY) gr_limiter=1;
  	double Wfac=lorentz_factor(CellP[target].gr.gamma,prims->v, gr_limiter);
  	double ucon[4]={Wfac/CellP[target].gr.alpha,0,0,0};
  	//V is the velocity seed by a Eulerian frame at rest
  	for(j=0;j<3;j++) ucon[j+1]=Wfac*(prims->v[j]-CellP[target].gr.beta[j]/CellP[target].gr.alpha);

#ifdef MAGNETIC
	rhoh += prims->bsq;
	Pstar += (prims->bsq/2);
	//Let's evaluate controvariant small b, to be used in Tmunu
	raise_lower_indices(CellP[target].gr.gamma,prims->v,v_cov);
	vdotB = v_cov[0]*prims->B[0] + v_cov[1]*prims->B[1] + v_cov[2]*prims->B[2];
	b_con[0] = Wfac*vdotB/CellP[target].gr.alpha;
	for(j=0;j<3;j++){b_con[j+1] = (prims->B[j]/Wfac) + (vdotB*ucon[j+1]);}  
#endif

  	//Now, build the stress tensor, starting from the rho*u*u terms
  	for(j=0;j<4;j++) {for(k=0;k<4;k++) {
		Tmunu[j][k] = rhoh*ucon[j]*ucon[k];
#ifdef MAGNETIC
	//Adding the -b^mu b^nu term
		Tmunu[j][k] -= b_con[j]*b_con[k];
#endif
		}
	}

  	//Now add the pressure term
  	//alpha^2 = -1/g^00 -> g^00 = -1/alpha^2
  	double alphainv2=1./(CellP[target].gr.alpha*CellP[target].gr.alpha);

  	Tmunu[0][0] += -Pstar*alphainv2;
  	//beta_i = g0i = gi0; g^0i = beta^i/alpha^2
  	for(j=0;j<3;j++) Tmunu[0][j+1]=Tmunu[j+1][0]=Tmunu[0][j+1]+Pstar*CellP[target].gr.beta[j]*alphainv2;

  	//Now matter terms (gamma^ij - beta^i beta^j/alpha^2)
  	Tmunu[1][1] += Pstar*(CellP[target].gr.ugamma[0] - CellP[target].gr.beta[0]*CellP[target].gr.beta[0]*alphainv2);
  	Tmunu[2][2] += Pstar*(CellP[target].gr.ugamma[3] - CellP[target].gr.beta[1]*CellP[target].gr.beta[1]*alphainv2);
  	Tmunu[3][3] += Pstar*(CellP[target].gr.ugamma[5] - CellP[target].gr.beta[2]*CellP[target].gr.beta[2]*alphainv2);

  	Tmunu[1][2] = Tmunu[2][1] = Tmunu[1][2]+Pstar*(CellP[target].gr.ugamma[1] - CellP[target].gr.beta[0]*CellP[target].gr.beta[1]*alphainv2);
  	Tmunu[1][3] = Tmunu[3][1] = Tmunu[1][3]+Pstar*(CellP[target].gr.ugamma[2] - CellP[target].gr.beta[0]*CellP[target].gr.beta[2]*alphainv2);
  	Tmunu[2][3] = Tmunu[3][2] = Tmunu[2][3]+Pstar*(CellP[target].gr.ugamma[4] - CellP[target].gr.beta[1]*CellP[target].gr.beta[2]*alphainv2);

	/* --- Compute the metric derivatives --- */
	get_metric_derivatives(target,Kij,d_gij,d_alpha,d_beta);

	//No mass sources (D term is 0)
	Src->D = 0.;

	//Momentum source term
	for(dir=0;dir<3;dir++)
	{
	   Src->S[dir]=0.;
	   for(j=0;j<4;j++){
              for(k=0;k<4;k++)
		{
		Src->S[dir]+=Tmunu[j][k]*d_gij[dir][j][k];
		}
	   }
#ifndef SELFGRAVITY_OFF
           Src->S[dir] += 2*Tmunu[0][0]*P[target].GravAccel[dir]; //This is the source term for the Newtonian self-gravity
#endif
	   Src->S[dir]*=0.5*CellP[target].gr.alpha*CellP[target].gr.gdet;
	}

	//Energy source term: note the extrinsic curvature Kij = -alpha*Gamma^0_ij
	Src->tau=0;
	for(j=0;j<3;j++)
	 for(k=0;k<3;k++)
	   Src->tau+= Kij[j][k] * (    Tmunu[0][0]*CellP[target].gr.beta[j]*CellP[target].gr.beta[k] +
				   2.0*Tmunu[0][j+1]*CellP[target].gr.beta[k] + Tmunu[j+1][k+1]);
	for(j=0;j<3;j++) Src->tau += (Tmunu[0][0]*CellP[target].gr.beta[j]+Tmunu[0][j+1])*(-d_alpha[j]);
#ifndef SELFGRAVITY_OFF
        //This is the source term for the Newtonian self-gravity
        for(j=0;j<3;j++) Src->tau += (Tmunu[0][0]*CellP[target].gr.beta[j]+Tmunu[0][j+1])*P[target].GravAccel[j]/CellP[target].gr.alpha;
#endif
	Src->tau *= CellP[target].gr.alpha*CellP[target].gr.gdet;


#if defined(MAGNETIC)
	for(k=0;k<3;k++) Src->B[k] = 0.;
#endif

}

void gr_add_source_term(int target,struct gr_primitives* prims,struct gr_conservatives* cons,double dtime, int do_implicit)
{
	int k;
	double res_con2prim;
	struct gr_conservatives Src;

#ifdef BOX_BND_PARTICLES
	if(P[target].ID==0){return;}//We don't want GR gravity for frozen boundary particles
#endif
	//Compute source terms
	get_source_term(target,prims,&Src);

	//Add the source terms via a predictor step

	for(k=0;k<3;k++) cons->S[k] += Src.S[k]*dtime;
	cons->tau  += Src.tau*dtime;

	//Now recover the updated primitives (prediction)
	int gr_limiter=0;
	double dpl[3]; for(k=0;k<3;k++) {dpl[k]=P[target].Pos[k];}
#if (GR_METRIC==2)
	for(k=0;k<3;k++) {dpl[k] = -P[target].Min_xyz_to_Sink[k];}
#endif
	if(dpl[0]*dpl[0]+dpl[1]*dpl[1]+dpl[2]*dpl[2] <= pow(GAMMA_LIMITER_RADIUS,2)) gr_limiter=1;
//	if(P[target].Pos[0]<LEFT_GAMMA_BOUNDARY||P[target].Pos[0]>RIGHT_GAMMA_BOUNDARY) gr_limiter=1;
	res_con2prim=gr_get_primitive_variables(cons,prims,&(CellP[target].gr), gr_limiter);
//	if(res_con2prim<0) {printf("Warning: failed inversion with code %g - particle %d at %g %g %g - kernelsize %g\n",res_con2prim, P[target].ID,P[target].Pos[0],P[target].Pos[1],P[target].Pos[2],P[target].KernelRadius);endrun(9993000);}
        if(res_con2prim<0) {printf("Warning: failed inversion with code in sources %g - ID %d - pos %g %g %g\n", res_con2prim,P[target].ID, P[target].Pos[0],P[target].Pos[1],P[target].Pos[2]);
#if (GR_METRIC==2)
                double r2_0=sqrt((P[target].Pos[0]-grP[0].Pos[0])*(P[target].Pos[0]-grP[0].Pos[0])+(P[target].Pos[1]-grP[0].Pos[1])*(P[target].Pos[1]-grP[0].Pos[1])+
                        (P[target].Pos[2]-grP[0].Pos[2])*(P[target].Pos[2]-grP[0].Pos[2]));
                double r2_1=sqrt((P[target].Pos[0]-grP[1].Pos[0])*(P[target].Pos[0]-grP[1].Pos[0])+(P[target].Pos[1]-grP[1].Pos[1])*(P[target].Pos[1]-grP[1].Pos[1])+
                        (P[target].Pos[2]-grP[1].Pos[2])*(P[target].Pos[2]-grP[1].Pos[2]));
                if(r2_0<2*EHradius(grP[0].Mass,grP[0].Spin) || r2_1<2*EHradius(grP[1].Mass,grP[1].Spin)){P[target].Mass=0; return;}
                else{endrun(9993000);}
#else
                endrun(9993000);
#endif
       }


	if(do_implicit)
	{
	 //Now correct the source term with a corrector step
	 struct gr_conservatives Src_corr=Src;

#ifdef IMPLICIT_SOLVER
	 double error=MAX_REAL_NUMBER;
	 double p_err[3];
	 for(k=0;k<3;k++) p_err[k] = cons->S[k]/cons->D;
	 //Copy Src into Src_corr to initialise it
	 Src_corr=Src;
	 int iter=0;

	 //Star implicit loop
	 do{
	  //Copy latest velocity value for error evaluation
	  for(k=0;k<3;k++) p_err[k] = cons->S[k]/cons->D;
	  //Copy last source evaluation for correction
	  Src = Src_corr;
#endif
	  //Get corrected source term
          get_source_term(target,prims,&Src_corr);

	  //Correct prediction
	  for(k=0;k<3;k++) cons->S[k] += (Src_corr.S[k]-Src.S[k])*dtime;
          cons->tau += (Src_corr.tau-Src.tau)*dtime;

	  //Get primitive quantities (corrector)
	  res_con2prim=gr_get_primitive_variables(cons,prims,&(CellP[target].gr), gr_limiter);
//	  if(res_con2prim<0) {printf("Warning: failed inversion with code %g - particle %d \n",res_con2prim, P[target].ID);endrun(9993001);}
          if(res_con2prim<0) {printf("Warning: failed inversion with code in sources %g - ID %d - pos %g %g %g\n", res_con2prim,P[target].ID, P[target].Pos[0],P[target].Pos[1],P[target].Pos[2]);
#if (GR_METRIC==2)
                double r2_0=sqrt((P[target].Pos[0]-grP[0].Pos[0])*(P[target].Pos[0]-grP[0].Pos[0])+(P[target].Pos[1]-grP[0].Pos[1])*(P[target].Pos[1]-grP[0].Pos[1])+
                        (P[target].Pos[2]-grP[0].Pos[2])*(P[target].Pos[2]-grP[0].Pos[2]));
                double r2_1=sqrt((P[target].Pos[0]-grP[1].Pos[0])*(P[target].Pos[0]-grP[1].Pos[0])+(P[target].Pos[1]-grP[1].Pos[1])*(P[target].Pos[1]-grP[1].Pos[1])+
                        (P[target].Pos[2]-grP[1].Pos[2])*(P[target].Pos[2]-grP[1].Pos[2]));
                if(r2_0<2*EHradius(grP[0].Mass,grP[0].Spin) || r2_1<2*EHradius(grP[1].Mass,grP[1].Spin)){P[target].Mass=0; return;}
                else{endrun(9993001);}
#else
                endrun(9993001);
#endif
         }



#ifdef IMPLICIT_SOLVER
	  //Compute error
	  error=0;
	  for(k=0;k<3;k++) {p_err[k]=cons->S[k]/cons->D-p_err[k]; if(p_err[k]*p_err[k]>error) error=p_err[k]*p_err[k];}
	  error = sqrt(cons->D*cons->D*error/(cons->S[0]*cons->S[0]+cons->S[1]*cons->S[1]+cons->S[2]*cons->S[2]));
	  iter++;
	  if(iter>IMPLICIT_MAX_NITER)
	  {
	   printf("WARNING: Iter >100 adding sources %2.15g %2.15g %2.15g %2.15g - %g %g %g %g\n",cons->S[0]/cons->D,
		cons->S[2]/cons->D,prims->rho0,prims->u,p_err[0],p_err[1],p_err[2],error);
	   break;
	  }
	 }while(error>IMPLICIT_SOLVER);
#endif
	}
}

#endif
