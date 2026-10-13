// function that actually computes the instantaneous fictitious acceleration due to stellar hardening (see details in Bonetti et al. 2016);
// this force is inteded to mimick the secular evolution of semi-major axis and eccentricity, i.e.:
// $\dot{a} = - a^2 G\rho H/\sigma$; $\dot{e} = a G\rho H K/\sigma$ (in the implementation below G is assumed to be 1)
// input parameters are positions (x,y,z), velocities (vx,vy,vz) and masses (m1,m2) of the two body involved;
// the function fills the two arrays containing the accelerations of body 1 (acc1_H) and body 2 (acc2_H), they must have length=3 (i.e. x,y,z components)
#ifdef GR_INCLUDE_STELLAR_HARDENING

#define N_Q 4
#define N_E 6
#define beta1_hard (0.38)
#define beta2_hard (0.055)
#define beta3_hard (8.03627856)
#define beta4_hard (0.148479509)
#define beta5_hard (1.89548651)
#define beta6_hard (0.220456858)

#define constant1_hard (4.*All.rho_hard*All.G/All.sigma_hard*M_PI)
#define constant2_hard (All.rho_hard*All.G/All.sigma_hard*M_PI)  ///sqrt(m1+m2);

//#define PREDCORR
double *q_HK;
double *e_K;

double *A_H;
double *a0_H;
double *g_H;

double **A_K;
double **a0_K;
double **g_K;
double **B_K;


double bilinear_interp(double* val_to_int, double* X_val, double* Y_val, double** values);
double linear_interp(double val_to_int, double* X_val, double* values);
void get_coefficients(double q, double e, double* A_H_interp, double* a0_H_interp, double* g_H_interp, double* A_K_interp, double* a0_K_interp, double* g_K_interp, double* B_K_interp);
void hardening_params();
double f1_ecc(double e);
double f2_ecc(double e);
double H_rate(const double aa_,double A_H_interp,double a0_H_interp,double g_H_interp, double hardening_radius);
double K_rate(const double aa_, double ee_,double A_K_interp,double a0_K_interp,double g_K_interp,double B_K_interp, double hardening_radius);
void update_hard_rate(const double t, double sma, double ecc, double hardening_radius, double* factor1_hard, double* factor2_hard);

void hardening_2B(const double r1[], const double r2[], const double v1[], const double v2[], const double m1, const double m2, double acc1_H[], double acc2_H[])
{
	// relative coordinates
	double X = r1[0]-r2[0];
	double Y = r1[1]-r2[1];
	double Z = r1[2]-r2[2];
	double Vx = v1[0]-v2[0];
	double Vy = v1[1]-v2[1];
	double Vz = v1[2]-v2[2];
	double R = sqrt( X * X + Y * Y + Z * Z );
	double V = sqrt( Vx * Vx + Vy * Vy + Vz * Vz );
	double V2 = V*V;

	double v_dot_r = X*Vx+Y*Vy+Z*Vz;
	double NORM = sqrt( R*R - v_dot_r*v_dot_r/V2 );
	double sma = -All.G*(m1+m2) / ( V2 - 2*All.G*(m1+m2)/R );
	double ecc = sqrt( (V2/(All.G*(m1+m2)) - 1./R)*(V2/(All.G*(m1+m2)) - 1./R)*R*R
	                    - V2*v_dot_r*v_dot_r/(All.G*(m1+m2)*All.G*(m1+m2))
	                    + 2*v_dot_r*v_dot_r/(All.G*R*(m1+m2)) );
        double factor1_hard, factor2_hard;
	
        // hardening radius (scale at which hardening generally proceeds at a constant rate)
	double hardening_radius = DMIN(m1,m2)/(4.*All.sigma_hard*All.sigma_hard);

        update_hard_rate(All.Time,sma,ecc,hardening_radius,&factor1_hard,&factor2_hard);

	double vx_dot_H_e = factor1_hard*v_dot_r/R*(X - v_dot_r*Vx/V2)/NORM;
	double vy_dot_H_e = factor1_hard*v_dot_r/R*(Y - v_dot_r*Vy/V2)/NORM;
	double vz_dot_H_e = factor1_hard*v_dot_r/R*(Z - v_dot_r*Vz/V2)/NORM;

	double vx_dot_H_a = -factor2_hard*V*R*Vx;
	double vy_dot_H_a = -factor2_hard*V*R*Vy;
	double vz_dot_H_a = -factor2_hard*V*R*Vz;

	acc1_H[0] = m2/(m1+m2)*(vx_dot_H_e + vx_dot_H_a);
	acc1_H[1] = m2/(m1+m2)*(vy_dot_H_e + vy_dot_H_a);
	acc1_H[2] = m2/(m1+m2)*(vz_dot_H_e + vz_dot_H_a);

	acc2_H[0] = -m1/(m1+m2)*(vx_dot_H_e + vx_dot_H_a);
	acc2_H[1] = -m1/(m1+m2)*(vy_dot_H_e + vy_dot_H_a);
	acc2_H[2] = -m1/(m1+m2)*(vz_dot_H_e + vz_dot_H_a);

}


// this function updates the variable used in the hardening implementation; best if called at each time step
void update_hard_rate(double t, double sma, double ecc, double hardening_radius, double* factor1_hard, double* factor2_hard)
{

        
        double q = DMIN(grP[0].Mass,grP[1].Mass)/DMAX(grP[0].Mass,grP[1].Mass);
        double A_H_interp, a0_H_interp, g_H_interp, A_K_interp, a0_K_interp, g_K_interp, B_K_interp;
        get_coefficients(q,ecc,&A_H_interp,&a0_H_interp,&g_H_interp,&A_K_interp,&a0_K_interp,&g_K_interp,&B_K_interp);

	double HH_hard = H_rate(sma,A_H_interp,a0_H_interp,g_H_interp,hardening_radius);
	double KK_hard = K_rate(sma, ecc,A_K_interp,a0_K_interp,g_K_interp,B_K_interp,hardening_radius);
#ifdef PREDCORR
	//All.aa_old_hard = sma, All.ee_old_hard = ecc;
	//sma = sma - sma*sma*All.rho_hard*All.HH_hard/All.sigma_hard*(t-All.t_old_hard);
	//ecc = ecc + sma*All.rho_hard*All.HH_hard*All.KK_hard/All.sigma_hard*(t-All.t_old_hard);
	//All.t_old_hard = t;
	//cout << "t = " << t << " told = " << t_old_hard << " aa = " << aa_hard << ", ee = " << ee_hard << endl;
#endif
	if(ecc <= 1e-5)
	{
		ecc = 0;
	}
	if(ecc >= 1)
	{
		ecc = 0.999;
	}

	*factor1_hard = 1./(2.*M_PI)*constant1_hard*HH_hard*KK_hard*sma*f1_ecc(ecc);
	*factor2_hard = 1./(2.*M_PI)*constant2_hard/sqrt(grP[0].Mass+grP[1].Mass)*HH_hard*sqrt(sma)*f2_ecc(ecc);

        //printf("factor1 %g, factor2 %g, H %g, K %g, q %g \n",*factor1_hard,*factor2_hard,HH_hard,KK_hard,q);
}

// auxilliary functions to get the orbital average form of the fictitious force
double f1_ecc(double e)
{
	if (e<=1e-2)
	{
		return 0;
	}
	else if(e >= 1)
	{
		printf("e_out > 1");
		exit(-1);
	}
	else
	{
		if(e >= 1)
		{
			printf("e_out > 1");
			e = 1-1e-5;
		}

		return 1./( e*pow(1.-e*e,beta1_hard)*pow(1.-pow(e,5.),beta2_hard) );
	}
}
// auxilliary functions to get the orbital average form of the fictitious force
double f2_ecc(double e)
{
	if(e >= 1)
	{
		printf("e_out > 1");
		exit(-1);
	}
	else if(e >= 1)
	{
		printf("e_out > 1");
		e = 1.-1e-5;
	}
	return pow(1.-pow(e,beta3_hard),beta4_hard)/pow(1.-pow(e,beta5_hard),beta6_hard);
}

// hardening rate
double H_rate(double aa_,double A_H_interp,double a0_H_interp,double g_H_interp, double hardening_radius)
{
	//printf("A_H_interp %g, a0_H_interp %g, g_H_interp %g, hardeningR %g \n",A_H_interp,a0_H_interp,g_H_interp,hardening_radius);
        return A_H_interp*pow(1.+aa_/(a0_H_interp*hardening_radius),g_H_interp);
}

// eccentricity growth rate
double K_rate(double aa_, double ee_,double A_K_interp,double a0_K_interp,double g_K_interp,double B_K_interp, double hardening_radius)
{
	if (ee_ == 0)
	{
		return 0.;
	}

	if(ee_ < e_K[0]) ee_ = e_K[0];
	if(ee_ > e_K[N_E-1]) ee_ = e_K[N_E-1];


	return A_K_interp*pow(1.+aa_/(a0_K_interp*hardening_radius),g_K_interp) + B_K_interp;
}

// this function initialises all the parameters used by the implementation of the hardening
void hardening_params()
{

	///////////////////////////////////////////////////////////////////////////////
	// the quantities below come from tables in Sesana,Haardt,Madau 2006 and are used to evaluate the dimensionless number H and K entering
	// a_dot and e_dot due to stellar hardening

        A_H = (double*) malloc(N_Q*sizeof(double));
        a0_H = (double*) malloc(N_Q*sizeof(double));
        g_H = (double*) malloc(N_Q*sizeof(double));
        q_HK = (double*) malloc(N_Q*sizeof(double));
        e_K = (double*) malloc(N_E*sizeof(double));
        A_K = (double**) malloc(N_E*sizeof(double*));
        a0_K = (double**) malloc(N_E*sizeof(double*));
        B_K = (double**) malloc(N_E*sizeof(double*));
        g_K = (double**) malloc(N_E*sizeof(double*));

	// values of mass-ratio
	int i, j;
	for (i = 0; i < N_E; ++i){
          A_K[i] = (double*) malloc(N_Q*sizeof(double));
          a0_K[i] = (double*) malloc(N_Q*sizeof(double));
          B_K[i] = (double*) malloc(N_Q*sizeof(double));
          g_K[i] = (double*) malloc(N_Q*sizeof(double));

        }
        
	for (i = 0; i < N_Q; ++i)
	{
		q_HK[i] = pow(3.,-i);
	}

	// values of eccentricity
	for (i = 0; i < N_E; ++i)
	{
		e_K[i] = 0.15*(i+1);
	}

	// tabulated values of parameters in functional form for H (see function H_rate)
	A_H[0]=14.55, A_H[1]=15.82, A_H[2]=17.17, A_H[3]=18.15;
	a0_H[0]=3.48, a0_H[1]=4.18, a0_H[2]=3.59, a0_H[3]=3.32;
	g_H[0]=-0.95, g_H[1]=-0.90, g_H[2]=-0.79, g_H[3]=-0.77;

	// tabulated values of parameters in functional form for K (see function K_rate)
	//e=0.15
	A_K[0][0]=0.037 ,  A_K[0][1]=0.082 ,  A_K[0][2]=0.051 ,  A_K[0][3]=0.064;
	a0_K[0][0]=0.339 , a0_K[0][1]=0.042 , a0_K[0][2]=0.385 , a0_K[0][3]=0.284;
	g_K[0][0]=-3.335 , g_K[0][1]=-0.168 , g_K[0][2]=-0.891 , g_K[0][3]=-1.206;
	B_K[0][0]=-0.012 , B_K[0][1]=-0.048 , B_K[0][2]=-0.011 , B_K[0][3]=0.021;

	//e=0.3
	A_K[1][0]=0.075 ,  A_K[1][1]=0.095 ,  A_K[1][2]=0.111 ,  A_K[1][3]=0.143;
	a0_K[1][0]=0.151 , a0_K[1][1]=0.213 , a0_K[1][2]=0.307 , a0_K[1][3]=1.033;
	g_K[1][0]=-1.548 , g_K[1][1]=-1.152 , g_K[1][2]=-1.107 , g_K[1][3]=-1.537;
	B_K[1][0]=-0.008 , B_K[1][1]=-0.012 , B_K[1][2]=-0.007 , B_K[1][3]=-0.021;

	//e=0.45
	A_K[2][0]=0.105 ,  A_K[2][1]=0.129 ,  A_K[2][2]=0.172 ,  A_K[2][3]=0.212;
	a0_K[2][0]=0.088 , a0_K[2][1]=0.137 , a0_K[2][2]=0.526 , a0_K[2][3]=0.722;
	g_K[2][0]=-0.893 , g_K[2][1]=-0.655 , g_K[2][2]=-1.174 , g_K[2][3]=-1.257;
	B_K[2][0]=-0.005 , B_K[2][1]=-0.006 , B_K[2][2]=-0.016 , B_K[2][3]=-0.022;

	//e=0.6
	A_K[3][0]=0.121 ,  A_K[3][1]=0.166 ,  A_K[3][2]=0.181 ,  A_K[3][3]=0.216;
	a0_K[3][0]=0.090 , a0_K[3][1]=0.081 , a0_K[3][2]=0.251 , a0_K[3][3]=0.430;
	g_K[3][0]=-0.895 , g_K[3][1]=-0.546 , g_K[3][2]=-1.169 , g_K[3][3]=-1.163;
	B_K[3][0]=-0.008 , B_K[3][1]=-0.006 , B_K[3][2]=-0.007 , B_K[3][3]=-0.014;

	//e=0.75
	A_K[4][0]=0.134 ,  A_K[4][1]=0.159 ,  A_K[4][2]=0.179 ,  A_K[4][3]=0.173;
	a0_K[4][0]=0.064 , a0_K[4][1]=0.079 , a0_K[4][2]=0.195 , a0_K[4][3]=0.771;
	g_K[4][0]=-0.544 , g_K[4][1]=-0.497 , g_K[4][2]=-0.846 , g_K[4][3]=-1.934;
	B_K[4][0]=-0.006 , B_K[4][1]=-0.010 , B_K[4][2]=-0.004 , B_K[4][3]=-0.014;

	//e=0.9
	A_K[5][0]=0.082 ,  A_K[5][1]=0.095 ,  A_K[5][2]=0.117 ,  A_K[5][3]=0.129;
	a0_K[5][0]=0.085 , a0_K[5][1]=0.122 , a0_K[5][2]=0.400 , a0_K[5][3]=0.329;
	g_K[5][0]=-0.663 , g_K[5][1]=-0.716 , g_K[5][2]=-1.170 , g_K[5][3]=-1.125;
	B_K[5][0]=-0.004 , B_K[5][1]=-0.008 , B_K[5][2]=-0.001 , B_K[5][3]=-0.020;

}

void get_coefficients(double q, double e, double* A_H_interp, double* a0_H_interp, double* g_H_interp, double* A_K_interp, double* a0_K_interp, double* g_K_interp, double* B_K_interp){
	// interpolation through tabulated values for H and K


        //interpolate on q
        *A_H_interp = linear_interp(q,q_HK,A_H);
        *a0_H_interp = linear_interp(q,q_HK,a0_H);
        *g_H_interp = linear_interp(q,q_HK,g_H);
       
        double pars[] = {e,q};
 
        //interpolate on q,e
	*A_K_interp = bilinear_interp(pars,e_K,q_HK,A_K);
	*a0_K_interp = bilinear_interp(pars,e_K,q_HK,a0_K);
	*g_K_interp = bilinear_interp(pars,e_K,q_HK,g_K);
	*B_K_interp = bilinear_interp(pars,e_K,q_HK,B_K);

}

double linear_interp(double val_to_int, double* X_val, double* values)
{
        int len1 = N_Q, I = len1-1;
        double xd, cc;
        int i=0;
        

        ////////////////////////////////
        for (i = 0; i < len1; ++i)
        {
                if(val_to_int > X_val[i]) // X_val (mass ratio) is encoded in reverse order in the tables
                {
                        I = i-1;
                        break;
                }
        }

        if(I == -1)
        {
                I = 0;
                xd = 0.;
        }
        else if(I == len1-1)
        {
                I = I-1;
                xd = 1.;
        }
        else
        {
                xd = (val_to_int-X_val[I])/(X_val[I+1]-X_val[I]);
        }
        
        cc = values[I]*(1.-xd) + values[I+1]*xd;
        return cc;

}

double bilinear_interp(double* val_to_int, double* X_val, double* Y_val, double** values)
{

        int len1 = N_E, I = len1-1;
        int len2 = N_Q, J = len2-1;
        double xd, yd, cc0, cc1, cc;
        int i=0,j=0;

        ////////////////////////////////
        for (i = 0; i < len1; ++i)
        {
                if(val_to_int[0] < X_val[i])
                {
                        I = i-1;
                        break;
                }
        }

        if(I == -1)
        {
                I = 0;
                xd = 0.;
        }
        else if(I == len1-1)
        {
                I = I-1;
                xd = 1.;
        }
        else
        {
                xd = (val_to_int[0]-X_val[I])/(X_val[I+1]-X_val[I]);
        }

        /////////////////////////////////
        for (j = 0; j < len2; ++j)
        {
                if(val_to_int[1] > Y_val[j]) // Y_val (mass ratio) is encoded in reverse order in the tables
                {
                        J = j-1;
                        break;
                }
        }

        if(J == -1)
        {
                J = 0;
                yd = 0.;
        }
        else if(J == len2-1)
        {
                J = J-1;
                yd = 1.;
        }
        else
        {
                yd = (val_to_int[1]-Y_val[J])/(Y_val[J+1]-Y_val[J]);
        }

        cc0 = values[I][J]  *(1.-xd) + values[I+1][J]  *xd;
        cc1 = values[I][J+1]*(1.-xd) + values[I+1][J+1]*xd;

        cc = cc0*(1.-yd) + cc1*yd;

        return cc;
}

#endif
