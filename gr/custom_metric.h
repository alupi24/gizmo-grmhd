#if defined(GENERAL_RELATIVITY) && (GR_METRIC<0)
//#define USE_ISOTROPIC_RADIUS
static double *rbins,*mass,*gtt,*riso;
static int NLTOV,log_scale;



void load_metric()
{
   FILE *fmetric;
   int i;
   double dummy;
   char fline[1000];
   fmetric=fopen("/cluster/gfedrigo/TOV/metric_table_iso.dat","r");

   if(fmetric==NULL)
   {
    	printf("Warning! Failed opening\n");
        endrun(9993002);
   }
   int check;
   fgets(fline,999,fmetric);
   sscanf(fline,"# %d %d",&NLTOV,&log_scale);
   rbins=(double*)malloc(sizeof(double)*NLTOV);
   mass=(double*)malloc(sizeof(double)*NLTOV);
   riso=(double*)malloc(sizeof(double)*NLTOV);
   gtt=(double*)malloc(sizeof(double)*NLTOV);
   for(i=0;i<NLTOV;i++)
   {
    	fgets(fline,999,fmetric);
        sscanf(fline,"%le %le %le %le",rbins+i,mass+i,gtt+i,riso+i);
   }

   fclose(fmetric);
   if(ThisTask==0){
	char spstr[4]="lin";
	if(log_scale) strcpy(spstr,"log");
	printf("Metric loaded (%d radii %s-spaced)\n",NLTOV,spstr); 
   }
}


void get_custom_metric_coefficients(double R2,double dp[3],double dp2[3],double gcov[])
{
  double r=sqrt(R2);
  int ir;
  double Mr;
#ifdef USE_ISOTROPIC_RADIUS
  double e4phi;
  if(r<=riso[0]) {r=riso[0];}

  if(log_scale) {
	if(r<riso[0]) ir=0; else ir=(log10(r/riso[0])/log10(riso[NLTOV-1]/riso[0])*(NLTOV-1));
  }
  else
  {
	int il=0;
	int ih=NLTOV-1;
	while(ih>il+1)
	{
	 ir=(il+ih)/2;
	 if(riso[ir]>r) ih=ir;
	 else il=ir;
	}
	if(riso[ir+1]<r) ir+=1;
	if(riso[ir]>r) ir=fmax(0,ir-1);
	if(!(riso[ir]<r && riso[ir+1]>=r) &&ir<NLTOV-1) {printf("Ahi! %d %g %g %g %g\n",ir,r,riso[ir-1],riso[ir],riso[ir+1]); getchar();}
  }
  //We assume linear interpolation within each bin for simplicity
  if(ir<NLTOV-1) Mr= (mass[ir+1]-mass[ir])/(riso[ir+1]-riso[ir])*(r-riso[ir])+mass[ir];
  else Mr=mass[NLTOV-1];
  if(ir<NLTOV-1) e4phi= (rbins[ir+1]-rbins[ir])/(riso[ir+1]-riso[ir])*(r-riso[ir])+rbins[ir];
  else e4phi=fmax((rbins[NLTOV-1]-rbins[NLTOV-2])/(riso[NLTOV-1]-riso[NLTOV-2])*(r-riso[NLTOV-2])+rbins[NLTOV-2],1.0);

  e4phi=e4phi*e4phi/r/r; //r^2/riso^2
  if(ir<NLTOV-1) gcov[0] = -((gtt[ir+1]-gtt[ir])/(riso[ir+1]-riso[ir])*(r-riso[ir])+gtt[ir]);
  else gcov[0] = -(1-2*Mr/r);

  gcov[1]=gcov[2]=gcov[3]=0.; //beta is 0 (no mixed terms t-j)
  //x-x
  gcov[4]=e4phi;
  //y-y
  gcov[7]=e4phi;
  //z-z
  gcov[9]=e4phi;
  //x-y
  gcov[5]=0;
  //x-z
  gcov[6]=0;
  //y-z
  gcov[8]=0;

#else //Use coordinate spherical radius

  if(r<=rbins[0]) {r=rbins[0];}

  if(log_scale) {
        if(r<rbins[0]) ir=0; else ir=(log10(r/rbins[0])/log10(rbins[NLTOV-1]/rbins[0])*(NLTOV-1));
  }
  else
  {
   	ir=r*(NLTOV-1)/rbins[NLTOV-1];
  }
  if(!(rbins[ir]<=r && r<=rbins[ir+1]) && ir<NLTOV-1) {printf("Ahi! %d %g %g %g\n",ir,rbins[ir],r,rbins[ir+1]); getchar();}

  //We assume linear interpolation within each bin for simplicity
  if(ir<NLTOV-1) Mr= (mass[ir+1]-mass[ir])/(rbins[ir+1]-rbins[ir])*(r-rbins[ir])+mass[ir];
  else Mr=mass[NLTOV-1];
  if(ir<NLTOV-1) gcov[0] = -((gtt[ir+1]-gtt[ir])/(rbins[ir+1]-rbins[ir])*(r-rbins[ir])+gtt[ir]);
  else gcov[0] = -(1-2*mass[NLTOV-1]/r);

  double grr=1./(1-2*Mr/r);
  double omega=2*Mr/(r*r*r);

  gcov[1]=gcov[2]=gcov[3]=0.; //beta is 0 (no mixed terms t-j)

  //x-x
  gcov[4]=grr*(1-omega*(dp2[1]+dp2[2]));
  //y-y
  gcov[7]=grr*(1-omega*(dp2[0]+dp2[2]));
  //z-z
  gcov[9]=grr*(1-omega*(dp2[0]+dp2[1]));
  //x-y
  gcov[5]=grr*dp[0]*dp[1]*omega;
  //x-z
  gcov[6]=grr*dp[0]*dp[2]*omega;
  //y-z
  gcov[8]=grr*dp[1]*dp[2]*omega;

#endif

}

#endif
