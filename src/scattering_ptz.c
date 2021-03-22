/*
   if (strcmp (scat_tp, "petzold") == 0) 
   {			
     fppf = fopen("anc/petzold_pf.dat", "r");					// Open phase function file
     if(fppf == NULL) 
     {
       printf ("ERROR: Failed to open Petzold data file anc/petzold_pf.dat\n");
       exit (-1);
     }
 
     for (ci = 0; ci < np_ptz; ci++) 
     {
       fscanf(fppf, "%lf %lf", &psi_ptz_pf[ci], &pfv_ptz[ci]);
       psi_ptz_pf[ci] = log(psi_ptz_pf[ci]);					// Log-transform the data
       pfv_ptz[ci] = log(pfv_ptz[ci]);
     }
     psi_ptz_c = exp(psi_ptz_pf[0]);						// Get critical angle
     coef_psi_ptz = exp(psi_ptz_pf[1]);						// Get polat angle coefficient
     coef_pf_ptz = exp(pfv_ptz[1]);						// Get phase function coefficient
     fclose(fppf);

     acc_ptz_pf = gsl_interp_accel_alloc ();					// Allocate GSL accelerator for interpolation
     itpl_ptz_pf = gsl_spline_alloc (gsl_interp_linear, np_ptz);		// Allocate GSL interpolator
     gsl_spline_init (itpl_ptz_pf, psi_ptz_pf, pfv_ptz, np_ptz);		// Initiate interpolator

     double scat_pf_ptz (const double scat_psi, const void *scat_par) 
     {										// Create scattering PF for Petzold
       double res = 0.0;
       if(scat_psi < psi_ptz_c) 
       {
         res = coef_pf_ptz * pow(coef_psi_ptz / scat_psi, 1.346);
       } else {
         res = exp(gsl_spline_eval (itpl_ptz_pf, log(scat_psi), acc_ptz_pf));
       } 
       return res;
     }

     scat_pf = scat_pf_ptz;
   
     fppf = fopen("anc/petzold_cdf_monotone.dat", "r");				// Open CDF interpolation file
     if(fppf == NULL) 
     {
       printf("ERROR: Failed to open Petzold data file "
         "anc/petzold_cdf_monotone.dat\n");
       exit(-1);
     }
  
     np_ptz = 66800;//75861							// Number of data points
     p_psi_ptz_cdf = (double*) calloc (np_ptz, sizeof(double));			// Allocate vector of polar angle points
     p_itg_ptz = (double*) calloc (np_ptz, sizeof(double));			// Allocate vector of PDF integral from 0 to psi

     for (ci = 0; ci < np_ptz; ci++)
     {
       fscanf(fppf, "%lf %lf", &p_itg_ptz[ci], &p_psi_ptz_cdf[ci]);
     }  
     p_psi_ptz_cdf[np_ptz - 1] = M_PI;						// Remove small numerical inaccuracies form readind PI from disk, that could cause GSL interpolation to fail

     acc_ptz_cdf = gsl_interp_accel_alloc ();					// Allocate GSL accelerator for interpolation
     itpl_ptz_cdf = gsl_spline_alloc (gsl_interp_linear, np_ptz);		// Allocate GSL interpolator
     gsl_spline_init (itpl_ptz_cdf, p_itg_ptz, p_psi_ptz_cdf, np_ptz);		// Initiate interpolator

     double scat_qtl_ptz (const double rud, const void *scat_par) 
     {										// Create scattering function for Petzold
       return gsl_spline_eval (itpl_ptz_cdf, rud, acc_ptz_cdf);
     }

     scat_qtl = scat_qtl_ptz;
   }

   // Truncated Petzold Scattering: 
   if (strcmp (scat_tp, "trc_petzold") == 0) 
   {
     fppf = fopen("anc/trc_petzold_pf.dat", "r");				// Open phase function file
     if(fppf == NULL)
     {
       printf("ERROR: Failed to open truncated Petzold data file "
         "anc/trunc_petzold_pf.dat\n");
       exit (-1);
     }
 
     ci = 0;
     while (ci < np_trcptz)
     {
       fscanf(fppf, "%lf %lf", &psi_trcptz_pf[ci], &pfv_trcptz[ci]);
       psi_trcptz_pf[ci] = log(psi_trcptz_pf[ci]);				// Log-transform the data
       pfv_trcptz[ci] = log(pfv_trcptz[ci]);
       ci++;
     }
     psi_trcptz_c = exp(psi_trcptz_pf[0]);					// Get critical angle
     pf_trcptz_c = exp(pfv_trcptz[0]);
     fclose(fppf);

     psi_trcptz_pf[0] = psi_trcptz_pf[0]-1E-6;					// Satisfy GSL interpolator for monotonic x
     psi_trcptz_pf[39] = log(M_PI);

     acc_ptz_pf = gsl_interp_accel_alloc ();					// Allocate GSL accelerator for interpolation
     itpl_ptz_pf = gsl_spline_alloc (gsl_interp_linear, np_trcptz); 		// Allocate GSL interpolator
     gsl_spline_init (itpl_ptz_pf, psi_trcptz_pf, pfv_trcptz, np_trcptz);	// Initiate interpolator

     double scat_pf_ptz (const double scat_psi, const void *scat_par) 	// Create scattering PF for Petzold
     {
       double res = 0.0;
       if(scat_psi < psi_trcptz_c)
       {
         res = pf_trcptz_c;
       } else {
         res = exp(gsl_spline_eval (itpl_ptz_pf, log(scat_psi), acc_ptz_pf));
       }
       return res;
     }
     scat_pf = scat_pf_ptz;

     fppf = fopen("anc/trc_petzold_cdf_monotone.dat", "r");			// Open CDF interpolation file
     if(fppf == NULL)
     {
       printf("ERROR: Failed to open truncated Petzold data file "
         "anc/trunc_petzold_cdf_monotone.dat\n");
       exit (-1);
     }
  
     np_trcptz = 33595;//32509;							// Number of data points
     p_psi_ptz_cdf = (double*) calloc (np_trcptz, sizeof(double));		// Allocate vector of polar angle points
     p_itg_ptz = (double*) calloc (np_trcptz, sizeof(double));			// Allocate vector of PDF integral from 0 to psi

     ci = 0;
     while (ci < np_trcptz)
     {
       fscanf(fppf, "%lf %lf", &p_itg_ptz[ci], &p_psi_ptz_cdf[ci]);
       ci++;
     }
     p_psi_ptz_cdf[np_trcptz - 1] = M_PI;						// Remove small numerical inaccuracies form readind PI from disk, that could cause GSL interpolation to fail

     acc_ptz_cdf = gsl_interp_accel_alloc ();					// Allocate GSL accelerator for interpolation
     itpl_ptz_cdf = gsl_spline_alloc (gsl_interp_linear, np_trcptz);		// Allocate GSL interpolator
     gsl_spline_init (itpl_ptz_cdf, p_itg_ptz, p_psi_ptz_cdf, np_trcptz);		// Initiate interpolator

     double scat_qtl_ptz (const double rud, const void *scat_par)		// Create scattering function for Petzold
     {
       return gsl_spline_eval (itpl_ptz_cdf, rud, acc_ptz_cdf);
     }

     scat_qtl = scat_qtl_ptz;
   }

*/
