
/*******************************************************************************
 scattering_ff.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 Version: 1.6
 Date: 2021-03-25
 License: GPL-3.0

 Implements the analytical phase function of Fournier & Forand (1994),
 with the modifications of Fournier & Miroslaw (1999). Note that neither
 the phase function (PF) or the cumulative distribution function (CDF) are
 used directly - they are used to build an interpolation table that is used
 as input for generic PF and CDF C functions. Also provided is a function to
 estimate the j and nh combination for a given backscattering fraction as
 presented by Mobley et al. (2002).

 REFERENCES:

 Fournier, G. R. & Forand, J. L. 1994. Analytic phase function for ocean
   water. SPIE Ocean Optics XII 2258, 8, 194-201.
   DOI: 10.1117/12.190063
 Fournier, G. R. & Jonasz, M. 1999. Computer-based underwater imaging
   analysis. Proceedings of SPIE 3761, 62-70.
   DOI: 10.1117/12.366488
 Mobley, C. D.; Sundman, L. K. & Boss, E. 2002. Phase function effects on
   oceanic light fields. Applied Optics 41, 6, 1035-1050.
   DOI: 10.1364/AO.41.001035

*******************************************************************************/

 #include <string.h>
 #include <math.h>
 #include <gsl/gsl_errno.h>
 #include <gsl/gsl_math.h>
 #include <gsl/gsl_min.h>
 #include <gsl/gsl_spline.h>

 #include "config.h"
 #include "constants.h"
 #include "aux.h"
 #include "scattering.h"
 #include "scattering_ff.h"
 #include "scattering_truncate.h"

/* scat_setup_ff

 Setup the Fournier-Forand scattering function.

 INPUT:
 scat     - Scattering parameters
            Pointer to scattering struct;
 fbb      - Pre-defined backscattering ratio for the FF function
            Constant double;
 method   - Method to retrieve scattering values for scattering functions that
            do not have analytical solution: "lut", "itp" or "min"
            Pointer to constant char;
 truncate - Flag to indicate if the diffraction peak of the phase function
            should be truncated (1) or not (0). Only applies to FF and Petzold.
            Constant int.

 OUTPUT:
 None. Allocate and initialize scattering parameters.

 */

 void
 scat_setup_ff
 (
   struct scattering * scat,
   double const scat_fbb, 
   int const truncate, 
   char const * method
 )
 {
   char stringa[STRMXLEN];
   char stringb[STRMXLEN];
   FILE *fi;

   scat->pfn = NAN;
   scat->fbb = scat_fbb;
   scat->truncate = truncate;
   strncpy(scat->method, method, STRMXLEN - 1);

//   scat->ff  = (struct scattering_ff*) malloc(sizeof(struct scattering_ff));
   scat->lut.npsi = SCAT_NPSI;

   // Get file name for ff with backscatter fraction scat_fbb:
   // Note that since only 4 digits are written, up to 4 digits only are 
   // guaranteed to generate different files.
   sprintf(stringa, "_%.4f", scat_fbb);
   strcpy(stringb, "ff");
   strcat(stringb, stringa);
   if ( truncate ) 
   {
     strcpy (stringa, stringb);
     strcpy (stringb, "trunc_");
     strcat (stringb, stringa);
   }
   strcpy (stringa, stringb);
   strcpy (stringb, SCAT_DIR);
   strcat (stringb, stringa);

   // Try to open the phase function, if not present, will be computed and 
   // saved:
   fi = fopen(stringb, "r");
   if ( !fi ) 
   {
     printf ("Creating phase function interpolation table...\n");
     fi = fopen(stringb, "w");
     if ( !fi ) 
     {
       printf ("\nERROR: Failed to create file %s\n", stringb);
       exit (-1);
     }

     // Find parameters j and iop_nh as in Mobley et al. (2002):
     scat_findpar_ff(scat_fbb, &scat->ff.j, &scat->ff.nh);

     double step_ff = M_PI / (double) (scat->lut.npsi - 1);
     double tmp_psi[2] = {0.0};
     double tmp_pf[1]  = {0.0};
     double tmp_cdf[2] = {0.0};
     for (int i = 1; i < scat->lut.npsi; i++)
     {
       scat->lut.psi[i] = scat->lut.psi[i - 1] + step_ff;
       tmp_psi[0] = scat->lut.psi[i];
       scat_pf_ff(scat, tmp_psi, tmp_pf);
       scat->lut.pf[i] = tmp_pf[0];
       scat_cdf_ff(scat, tmp_psi, tmp_cdf);
       scat->lut.cdf[i] = tmp_cdf[0]; 
     }
     scat->lut.psi[scat->lut.npsi - 1] = M_PI;
     scat->lut.cdf[scat->lut.npsi - 1] = 1.0;

     if ( truncate )
     {
       printf ("WARNING: Bias assessment of the Phase Function diffraction " 
         "peak truncation of the diffraction peak of the PF only validated "
         "for Petzold PF for Lu estimate under diffuse and direct light "
         "between 0 and 80º. It might be interesting to perform a validation "
         "run if not in this condition.\n");
       scat_trunc_ff (scat);
     }

     for (int i = 0; i < scat->lut.npsi; i++)
     {
       fprintf(fi, "%.12e\t%.12e\t%.12e\n", scat->lut.psi[i],
         scat->lut.pf[i], scat->lut.cdf[i]);
     }
   }
   else
   {
     // If scattering file exists, read from disk:
     for (int i = 0; i < scat->lut.npsi; i++) 
     {
       fscanf(fi, "%lf\t%lf\t%lf\n", &scat->lut.psi[i], &scat->lut.pf[i], 
         &scat->lut.cdf[i]);
     }
     scat->lut.psi[scat->lut.npsi - 1] = M_PI;
     scat->lut.cdf[scat->lut.npsi - 1] = 1.0;
   }
   fclose(fi);

   scat->lut.kpsi = M_PI / (scat->lut.npsi - 1.0);
   if ( strcmp(method, "lut") == 0 )
   {
     // If requested to solve based on direct access to LUT:
     scat->fun.qtl = scat_qtl_ff_lut;
     scat->fun.pf  = scat_pf_ff_lut;
   } 
   else if ( strcmp(method, "itp") == 0 )
   {

//     scat->itp       = scat_itp_alloc();
     // If requested to solve based on interpolation:
     scat->itp.accl_qtl = gsl_interp_accel_alloc();
     scat->itp.accl_cdf = gsl_interp_accel_alloc();
     scat->itp.accl_pf  = gsl_interp_accel_alloc();

     scat->itp.itpl_qtl = gsl_spline_alloc(gsl_interp_cspline, 
       scat->lut.npsi);
     scat->itp.itpl_cdf = gsl_spline_alloc(gsl_interp_cspline, 
       scat->lut.npsi);
     scat->itp.itpl_pf  = gsl_spline_alloc(gsl_interp_cspline,
       scat->lut.npsi);

     gsl_spline_init (scat->itp.itpl_qtl, scat->lut.cdf, scat->lut.psi, 
       scat->lut.npsi);
     gsl_spline_init (scat->itp.itpl_cdf, scat->lut.psi, scat->lut.cdf,
       scat->lut.npsi);
     gsl_spline_init (scat->itp.itpl_pf,  scat->lut.psi, scat->lut.pf,
       scat->lut.npsi);

     scat->fun.qtl = &scat_qtl_ff_itp;
     scat->fun.cdf = &scat_cdf_ff_itp;
     scat->fun.pdf = &scat_pdf_ff;
     scat->fun.pf  = &scat_pf_ff_itp;
   }
   scat->g = scat_av_mu(scat);
 }

/* scat_findpar_ff

 Find FF parameters (junge slope j and refractive index of hydrosols nh) for a
 given backscatter fraction as in Mobley et al. (2002). It requires the
 auxiliary function findv (a function to be minimized in the search for the
 parameter combination that give the same backscatter fraction as requested).
 
 INPUT:
 v        - Centered and scaled value of the junge slope;
 params   - GSL function compatible generic parameter;
 scat_fbb - Backscattering fraction (unitless), [0.0001, 0.15];
 j        - Pointer to power coefficient of the Junge distribution;
 iop_nh   - Pointer to refractive index of hydrosols relative to water;

 OUTPUT:
 None. Update the values pointed by j and iop_nh.

 */

 static double
 findv
 (
   double v,
   void   * params
 )
 {	
   double fbb_r = *(double *) params;		// Requested scat_fbb value, as passed in params
   double fbb_v;				// scat_fbb evaluated at the input v
   double dpi = 0.0;				// delta factor for a given iop_nh, here restricted for a fixed relation to v
   
   // dpi dependency on v, as provided by Eq. 8 in Mobley et al. 2002:
   dpi   = 2.0 / (3.0 * pow (0.01 - 0.3084 * v, 2.0));
   fbb_v = 1.0 - ((1.0 - pow(dpi, v + 1.0) - 0.5 * (1.0 - pow(dpi, v))) / 
                 ((1.0 - dpi) * pow(dpi, v)));

   return ABS(fbb_r - fbb_v);
 }

 void
 scat_findpar_ff
 (
   double const scat_fbb,
   double * j,
   double * iop_nh
 )
 {
   int    fmin_status;				// Status of the GSL minimizer
   int    iter = 0, max_iter = 1000;		// Iteration and maximum interation for the GSL minimizer
   double v = (3.0 - 3.5835) / 2.0;		// First estimate of v, with j as for Petzold as defined in Mobley et al. 2002
   double lv = -0.5, uv = 0;			// Parameter space to search for v. It include possible values for iop_nh from 1.01 to 1.15
   const  gsl_min_fminimizer_type *bfmin_t;	// Pointer to GSL minimizer type
   gsl_min_fminimizer *fmin;			// Pointer to GSL minimizer   
   gsl_function gslf;				// GSL function structure
   double fbb_r = scat_fbb;

   if ( scat_fbb < 0.0001 || scat_fbb > 0.15 ) 
   {
     printf("scat_fbb (backscattering fraction for the Fournier-Forand phase "
       "function) outside valid interval [0.0001, 0.15]");
     exit(-1);
   }
  
   gslf.function = &findv;
   gslf.params = &fbb_r;
  
   bfmin_t = gsl_min_fminimizer_brent;
   fmin    = gsl_min_fminimizer_alloc (bfmin_t);
   gsl_min_fminimizer_set (fmin, &gslf, v, lv, uv);
  
   do 
   {
     iter++;
     fmin_status = gsl_min_fminimizer_iterate (fmin);

     v = gsl_min_fminimizer_x_minimum (fmin);
     lv = gsl_min_fminimizer_x_lower (fmin);
     uv = gsl_min_fminimizer_x_upper (fmin);
     fmin_status = gsl_min_test_interval (lv, uv, 0.00001, 0.0);
   } while (fmin_status == GSL_CONTINUE && iter < max_iter);

   gsl_min_fminimizer_free (fmin);

   *j  = (-2.0 * v) + 3.0;
   *iop_nh = 1.01 + 0.1542 * (*j - 3.0);
 }

/* scat_cdf_ff, scat_pdf_ff, scat_pf_ff

 Fournier-Forand phase function.

 INPUT:
 scat_psi - Polar scattering angle in radians, [0, PI].
 j        - Power coefficient of the Junge distribution (unitless), 
            (-Inf, Inf)
 iop_nh   - Real part of the refractive index of hydrosols relative to water 
            (unitless), (1, ?]

 OUTPUT:
 Returns the normalized scattering magnitude at the requested inputs.

 */

 void
 scat_cdf_ff
 (
   struct scattering const * scat,
   double const * s_scat,
   double * prob
 )
 {
   double term1 = 0.0, term2 = 0.0;
   double hspsi = sin(s_scat[0] / 2.0);
   double u = 0.0, d = 0.0, dpi = 0.0;
   double v = 0.0, dv  = 0.0, dpiv = 0.0; 

   if ( s_scat[0] < TOLERANCE )
   {
     prob[0] = 0.0;
   }
   else
   {
     double j = scat->ff.j;
     double iop_nh = scat->ff.nh;

     v   = (3.0 - j) / 2.0;
     d   = (4.0 / (3.0 * ((iop_nh - 1) * (iop_nh - 1)))) * (hspsi * hspsi);
     dpi = (4.0 / (3.0 * ((iop_nh - 1) * (iop_nh - 1)))) * (sin(M_PI_2) * 
           sin(M_PI_2));
     dv  = pow(d, v);

     dpiv  = pow(dpi, v);
     u     = 2.0 * hspsi;
     term1 = (1.0 / ((1.0 - d) * dv)) * ((1.0 - pow(d, (v + 1.0))) - 
             (((u * u) / 4.0) * (1.0 - dv)));
     term2 = (1.0 / 8.0) * ((1.0 - dpiv) / ((dpi - 1.0) * dpiv)) *  
             cos(s_scat[0]) * (sin(s_scat[0]) * sin(s_scat[0]));

     prob[0] = term1 + term2;
   }
   prob[1] = K_1_2PI * s_scat[1];
 }

 void
 scat_pdf_ff
 (
   struct scattering const * scat,
   double const * s_scat,
   double * dens
 )
 {
   scat_pf(scat, s_scat, dens);
   dens[0] *= K_2PI * sin(s_scat[0]);
   dens[1] = K_1_2PI;
 }

 void
 scat_pf_ff
 (
   struct scattering const * scat,
   double const * s_scat,
   double * pf
 )
 {
   double term1 = 0.0, term2 = 0.0;						
   double hspsi = sin(s_scat[0] / 2.0);						
   double u = 0.0, upi = 0.0, d = 0.0, dpi = 0.0; 				
   double v = 0.0, dv  = 0.0, dpiv = 0.0; 
   double mu = cos(s_scat[0]);

   double j = scat->ff.j;
   double iop_nh = scat->ff.nh;

   v   = (3.0 - j) / 2.0;
   d   = (4.0 / (3.0 * ((iop_nh - 1) * (iop_nh - 1)))) * (hspsi * hspsi);
   dpi = (4.0 / (3.0 * ((iop_nh - 1) * (iop_nh - 1)))) * (sin(M_PI_2) * 
         sin(M_PI_2));
   dv  = pow(d, v);
   dpiv  = pow(dpi, v);
   u     = 2.0 * hspsi;
   upi   = 2.0 * sin(M_PI_2);
   term1 = (1.0 / (4.0 * M_PI)) * (1.0 / (((1.0 - d) * (1.0 - d)) * dv)) * 
           (((v * (1.0 - d)) - (1.0 - dv)) + ((4.0 / (u * u)) * 
           ((d * (1.0 - dv)) - (v * (1.0 - d)))));
   term2 = (1.0 / (4.0 * M_PI)) * 
           (1.0 / (((1.0 - dpi) * (1.0 - dpi)) * dpiv)) * 
           (((v * (1.0 - dpi)) - (1.0 - dpiv)) + ((4.0 / (upi * upi)) * 
           ((dpi * (1.0 - dpiv)) - (v * (1.0 - dpi)))));

   pf[0] = term1 + (term2 * ((mu * mu) / 4.0));
 }


/* Function: scat_qtl_ff_itp, scat_cdf_ff_itp, scat_pdf_ff_itp, scat_pf_ff_itp

 Finds the normalized scattering magnitude for the FF at a given angle by 
 interpolation.

 INPUT:
 scat_psi  - Scattering angle, [0,PI]
 scat_par - Scattering parameter structure, in this case is the str_scatitp
             structure defined in the 'scat.h' header file.

 OUTPUT:
 The normalized phase function at the requested scattering angle.

 */

 void
 scat_qtl_ff_itp
 (
   struct scattering const * scat,
   double * s_scat,
   double const * prob
 )
 {
   s_scat[0] = gsl_spline_eval(scat->itp.itpl_qtl, prob[0], 
     scat->itp.accl_qtl);
   s_scat[1] = K_2PI * prob[1];
 }

 void
 scat_cdf_ff_itp
 (
   struct scattering const *scat,
   double const *s_scat,
   double *prob
 )
 {
   prob[0] = gsl_spline_eval(scat->itp.itpl_cdf, s_scat[0], 
     scat->itp.accl_cdf);
   prob[1] = K_1_2PI * s_scat[1];
 } 

 void
 scat_pdf_ff_itp
 (
   struct scattering const *scat,
   double const *s_scat,
   double *dens
 )
 {
   dens[0] = gsl_spline_eval(scat->itp.itpl_pf, s_scat[0], scat->itp.accl_pf);
   dens[0] *= K_2PI;
   dens[1] = K_1_2PI;
 } 

 void
 scat_pf_ff_itp
 (
   struct scattering const *scat,
   double const *s_scat,
   double *pf
 )
 {
   pf[0] = gsl_spline_eval(scat->itp.itpl_pf, s_scat[0], scat->itp.accl_pf);
 } 

/* Function: scat_qtl_ff_lut, scat_cdf_ff_lut, scat_pdf_ff_lut, scat_pf_ff_lut

 Finds the normalized scattering magnitude for the FF at a given angle by 
 interval. We can use findBin here because the intervals in scat_psi are constant.

 INPUT:
 scat_psi  - Scattering angle, [0,PI]
 scat_par - Scattering parameter structure, in this case is the str_scattab
              structure defined in the 'scat.h' header file.
 
 OUTPUT:
 The normalized phase function at the requested scattering angle.

 */

 void
 scat_qtl_ff_lut
 (
   struct scattering const *scat,
   double *s_scat,
   double const *prob
 )
 {
   printf("NOT IMPLEMENTED LUT\n");
   exit(-1);
 }

 void
 scat_cdf_ff_lut
 (
   struct scattering const *scat,
   double const *s_scat,
   double *prob
 )
 {
   printf("NOT IMPLEMENTED LUT\n");
   exit(-1);
 } 

 void
 scat_pdf_ff_lut
 (
   struct scattering const *scat,
   double const *s_scat,
   double *dens
 )
 {
   printf("NOT IMPLEMENTED LUT\n");
   exit(-1);
 } 

 void
 scat_pf_ff_lut
 (
   struct scattering const *scat,
   double const *s_scat,
   double *pf
 )
 {
   printf("NOT IMPLEMENTED LUT\n");
   exit(-1);
 } 

/* Function: scat_qtl_ff_min, scat_cdf_ff_min, scat_pdf_ff_min, scat_pf_ff_min

 Finds the normalized scattering magnitude for the FF at a given angle by 
 interval. We can use findBin here because the intervals in scat_psi are constant.

 INPUT:
 scat_psi  - Scattering angle, [0,PI]
 scat_par - Scattering parameter structure, in this case is the str_scattab
              structure defined in the 'scat.h' header file.
 
 OUTPUT:
 The normalized phase function at the requested scattering angle.

 */

 void
 scat_qtl_ff_min
 (
   struct scattering const *scat,
   double *s_scat,
   double const *prob
 )
 {
   printf("NOT IMPLEMENTED LUT\n");
   exit(-1);
 }

 void
 scat_cdf_ff_min
 (
   struct scattering const *scat,
   double const *s_scat,
   double *prob
 )
 {
   printf("NOT IMPLEMENTED LUT\n");
   exit(-1);
 } 

 void
 scat_pdf_ff_min
 (
   struct scattering const *scat,
   double const *s_scat,
   double *dens
 )
 {
   printf("NOT IMPLEMENTED LUT\n");
   exit(-1);
 } 

 void
 scat_pf_ff_min
 (
   struct scattering const *scat,
   double const *s_scat,
   double *pf
 )
 {
   printf("NOT IMPLEMENTED LUT\n");
   exit(-1);
 } 


