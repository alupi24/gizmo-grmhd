
#ifdef GR_INCLUDE_PN_CORRECTIONS

void PN_acceleration(const double r1[],const double r2[], const double v1in[], const double v2in[], double m1, double m2, double a1[], double a2[])
{
    double G = All.G, G_2 = G*G, G_3 = G_2*G;
    double c = 1., c_2 = c*c, c_4 = c_2*c_2, c_5 = c_4*c;
    int pn1 = (GR_INCLUDE_PN_CORRECTIONS)%2, pn2 = (GR_INCLUDE_PN_CORRECTIONS>>1)%2, pn25 = (GR_INCLUDE_PN_CORRECTIONS>>2)%2;
    double vxcom = (v1in[0]*m1+v2in[0]*m2)/(m1+m2), vycom = (v1in[1]*m1+v2in[1]*m2)/(m1+m2), vzcom = (v1in[2]*m1+v2in[2]*m2)/(m1+m2);
    double v1[] = {v1in[0]-vxcom, v1in[1]-vycom, v1in[2]-vzcom};
    double v2[] = {v2in[0]-vxcom, v2in[1]-vycom, v2in[2]-vzcom};
    //double v1[] = {v1in[0], v1in[1], v1in[2]};
    //double v2[] = {v2in[0], v2in[1], v2in[2]};

    double m1_2 = m1*m1, m2_2 = m2*m2;
    double m1_3 = m1*m1_2, m2_3 = m2*m2_2;

    double r12_vec[] = {r1[0]-r2[0],r1[1]-r2[1],r1[2]-r2[2]}; 
    double v12_vec[] = {v1[0]-v2[0],v1[1]-v2[1],v1[2]-v2[2]}; 

    double r12 = sqrt(r12_vec[0]*r12_vec[0] + r12_vec[1]*r12_vec[1] + r12_vec[2]*r12_vec[2]);
    //printf("r12 = %g m1 = %g m2 = %g\n",r12,m1,m2);
    double n12_vec[] = {r12_vec[0]/r12, r12_vec[1]/r12, r12_vec[2]/r12};

    double v12 = sqrt(v12_vec[0]*v12_vec[0] + v12_vec[1]*v12_vec[1] + v12_vec[2]*v12_vec[2]);
    double v12_2 = v12*v12;

    double v1_mod = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
    double v2_mod = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);

    double v1_2 = v1_mod*v1_mod, v2_2 = v2_mod*v2_mod;
    double v1_4 = v1_2*v1_2, v2_4 = v2_2*v2_2;

    double v1_v2 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    double v1_v2_2 = v1_v2*v1_v2;

    double r12_2 = r12*r12, r12_3 = r12_2*r12, r12_4 = r12_2*r12_2;

    double n12_v1 = n12_vec[0]*v1[0] + n12_vec[1]*v1[1] + n12_vec[2]*v1[2];
    double n12_v2 = n12_vec[0]*v2[0] + n12_vec[1]*v2[1] + n12_vec[2]*v2[2];
    double n12_v12 = n12_vec[0]*v12_vec[0] + n12_vec[1]*v12_vec[1] + n12_vec[2]*v12_vec[2];

    double n12_v1_2 = n12_v1*n12_v1, n12_v2_2 = n12_v2*n12_v2;
    double n12_v1_3 = n12_v1_2*n12_v1, n12_v2_3 = n12_v2_2*n12_v2;
    double n12_v1_4 = n12_v1_2*n12_v1_2, n12_v2_4 = n12_v2_2*n12_v2_2;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // body 1

    double term_N_body1_n, term_1PN_body1_n, term_1PN_body1_v, term_2PN_body1_n, term_2PN_body1_v, term_2_5PN_body1_n, term_2_5PN_body1_v;

    // Newtonian term
    term_N_body1_n = -G*m2/r12_2;


    // 1PN terms
    term_1PN_body1_n = 5*G_2*m1*m2/r12_3 + 4*G_2*m2_2/r12_3 + G*m2/r12_2*( 3./2.*n12_v2_2 - v1_2 + 4*v1_v2 - 2*v2_2 );

    term_1PN_body1_v = G*m2/r12_2*( 4*n12_v1 - 3*n12_v2 );


    // 2PN terms
    term_2PN_body1_n = 
    -57./4.*G_3*m1_2*m2/r12_4 - 69./2.*G_3*m1*m2_2/r12_4 - 9*G_3*m2_3/r12_4
    +G*m2/r12_2*( -15./8.*n12_v2_4 + 3./2.*n12_v2_2*v1_2 - 6*n12_v2_2*v1_v2 - 2*v1_v2_2 + 9./2.*n12_v2_2*v2_2 + 4*v1_v2*v2_2 - 2*v2_4 )
    +G_2*m1*m2/r12_3*( 39./2.*n12_v1_2 - 39*n12_v1*n12_v2 + 17./2.*n12_v2_2 - 15./4.*v1_2 - 5./2.*v1_v2 + 5./4.*v2_2 )
    +G_2*m2_2/r12_3*( 2*n12_v1_2 - 4*n12_v1*n12_v2 - 6*n12_v2_2 - 8*v1_v2 + 4*v2_2 );

    term_2PN_body1_v = 
    G_2*m2_2/r12_3*( -2*n12_v1 - 2*n12_v2 ) + G_2*m1*m2/r12_3*( -63./4.*n12_v1 + 55./4.*n12_v2 )
    +G*m2/r12_2*( -6*n12_v1*n12_v2_2 + 9./2.*n12_v2_3 + n12_v2*v1_2 - 4*n12_v1*v1_v2 + 4*n12_v2*v1_v2 + 4*n12_v1*v2_2 - 5*n12_v2*v2_2 );


    // 2.5PN terms
    term_2_5PN_body1_n = 208./15.*G_3*m1*m2_2/r12_4*n12_v12 - 24./5.*G_3*m1_2*m2/r12_4*n12_v12 + 12./5.*G_2*m1*m2/r12_3*n12_v12*v12_2;

    term_2_5PN_body1_v = 8./5.*G_3*m1_2*m2/r12_4 - 32./5.*G_3*m1*m2_2/r12_4 - 4./5.*G_2*m1*m2/r12_3*v12_2;

    ///////////////////////////
    for (int i = 0; i < 3; i++)
    {
        a1[i] = (0.0*term_N_body1_n*n12_vec[i] //first term already present in analytic_gravity
              + pn1/c_2*( term_1PN_body1_n*n12_vec[i]   + term_1PN_body1_v*v12_vec[i] )
              + pn2/c_4*( term_2PN_body1_n*n12_vec[i]   + term_2PN_body1_v*v12_vec[i] )
              + pn25/c_5*( term_2_5PN_body1_n*n12_vec[i] + term_2_5PN_body1_v*v12_vec[i] ))/G;
    //printf("PN_newton (a1) = %1.5g \n",term_N_body1_n*n12_vec[i]);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // body 2

    double term_N_body2_n, term_1PN_body2_n, term_1PN_body2_v, term_2PN_body2_n, term_2PN_body2_v, term_2_5PN_body2_n, term_2_5PN_body2_v;

    // Newtonian term
    term_N_body2_n = -G*m1/r12_2;


    // 1PN terms
    term_1PN_body2_n = 5*G_2*m1*m2/r12_3 + 4*G_2*m1_2/r12_3 + G*m1/r12_2*( 3./2.*n12_v1_2 - v2_2 + 4*v1_v2 - 2*v1_2 );

    term_1PN_body2_v = G*m1/r12_2*( -4*n12_v2 + 3*n12_v1 );


    // 2PN terms
    term_2PN_body2_n = 
    -57./4.*G_3*m2_2*m1/r12_4 - 69./2.*G_3*m2*m1_2/r12_4 - 9*G_3*m1_3/r12_4
    +G*m1/r12_2*( -15./8.*n12_v1_4 + 3./2.*n12_v1_2*v2_2 - 6*n12_v1_2*v1_v2 - 2*v1_v2_2 + 9./2.*n12_v1_2*v1_2 + 4*v1_v2*v1_2 - 2*v1_4 )
    +G_2*m1*m2/r12_3*( 39./2.*n12_v2_2 - 39*n12_v2*n12_v1 + 17./2.*n12_v1_2 - 15./4.*v2_2 - 5./2.*v1_v2 + 5./4.*v1_2 )
    +G_2*m1_2/r12_3*( 2*n12_v2_2 - 4*n12_v2*n12_v1 - 6*n12_v1_2 - 8*v1_v2 + 4*v1_2 );

    term_2PN_body2_v = 
    G_2*m1_2/r12_3*( 2*n12_v2 + 2*n12_v1 ) + G_2*m1*m2/r12_3*( 63./4.*n12_v2 - 55./4.*n12_v1 )
    +G*m1/r12_2*( 6*n12_v2*n12_v1_2 - 9./2.*n12_v1_3 - n12_v1*v2_2 + 4*n12_v2*v1_v2 - 4*n12_v1*v1_v2 - 4*n12_v2*v1_2 + 5*n12_v1*v1_2 );


    // 2.5PN terms
    term_2_5PN_body2_n = 208./15.*G_3*m2*m1_2/r12_4*n12_v12 - 24./5.*G_3*m2_2*m1/r12_4*n12_v12 + 12./5.*G_2*m1*m2/r12_3*n12_v12*v12_2;

    term_2_5PN_body2_v = 8./5.*G_3*m2_2*m1/r12_4 - 32./5.*G_3*m2*m1_2/r12_4 - 4./5.*G_2*m1*m2/r12_3*v12_2;

    ///////////////////////////
    for (int i = 0; i < 3; i++)
    {
        a2[i] = (-0.0*term_N_body2_n*n12_vec[i] //first term already present in analytic_gravity
              + pn1/c_2*( -term_1PN_body2_n*n12_vec[i]   - term_1PN_body2_v*v12_vec[i] )
              + pn2/c_4*( -term_2PN_body2_n*n12_vec[i]   - term_2PN_body2_v*v12_vec[i] )
              + pn25/c_5*( -term_2_5PN_body2_n*n12_vec[i] - term_2_5PN_body2_v*v12_vec[i] ))/G;
    //printf("PN_newton (a2) = %1.5g \n",term_N_body2_n*n12_vec[i]);
    }
    
}

#endif
