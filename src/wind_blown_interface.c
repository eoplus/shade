

 void
 variance
 (
   double * sigma_u_sq,
   double * sigma_c_sq,
   double const wspd
 )
 {
   *sigma_u_sq = WIND_K_U * wspd;
   *sigma_c_sq = WIND_K_C * wspd;

   double N[3]; // Outward normal of the facet

   sigma = WIND_DELTA * 
   v1 = gsl_ran_gaussian(random, sigma);
   v2 = gsl_ran_gaussian(random, sigma);
   v3 = gsl_ran_gaussian(random, sigma);

   norm = sqrt(1.0 + nu_u^2 + nu_c^2);
   N[0] = -nu_u / norm;
   N[1] = -nu_c / norm;
   N[2] = -1.0  / norm;
 }
