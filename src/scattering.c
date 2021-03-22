
/*******************************************************************************
 scattering.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 Provides functions describing scattering functions and to take random samples
 of their probability distribution function.

 Two broad classes of scattering are molecular (Rayleigh) and particulate
 scattering. Rayleigh scattering occurs when interaction centers are much
 smaller than the wavelength of light and are influenced by turbulence and
 Brownian motion.

 All scattering functions are described by a set of parameters presented in
 detail in the scattering.h header file. Some parameters are unique to a
 scattering type or sampling function (e.g., for function without analytical 
 solution to scattering angle, a look-up table, interpolation, or minimization
 can be implemented). Inelastic scattering at the moment is not implemented and
 may be described by isotropic point sources, as described in the sources.c and 
 sources.h files. 

 Finally, for the PFs with a diffraction peak (e.g., FF and Petzold), the phase
 function can be "truncated" in the sense that the scattering in the diffraction
 peak is substituted by the average scattering in that diffraction peak,
 therefore without altering its integral. This is intended to reduce the
 variance of the MC simulations and is explained in details in the
 'scat_trunc.c' source file.

 In this source file, the general functions for scattering and the specific
 functions for isotropic and single term Henyey-Greestein are provided. Other
 functions are provided in specific source files named scattering_x.c, where "x"
 is the scattering type code.

 THEORY:
 As for sources, we normalize the scattering phase function to integrate to unit
 as:

 Integral{ dphi
   Integral{ dpsi
     B(psi,phi) sin(psi) } } = 1.                                           (A1)

 We use the convention that if an integral does not specifies bounds, the 
 integration is to be carried over the full range of the variable. Those ranges 
 are [0,2PI] for the azimuth angle phi and [0,PI] for the polar angle psi.
 The normalization condition results in the scattering phase function (PF),
 B(psi,phi), being a probability function, with Eq. 1 being the Cumulative
 Distribution Function (CDF).

 It is possible to separate this bi-dimensional probability function into the 
 probabilities of its components for the polar and azimuth angles:

 PDF(psi) = S(psi,phi) sin(psi) / Integrate{S(psi,phi) sin(psi) dpsi},      (A2)

 CDF(psi) = Integrate(0 to psi){PDF(psi) dpsi},                             (A3)

 PDF(phi) = Integrate{S(psi,phi) sin(psi) dpsi} / CDF(psi,phi),             (A4)

 CDF(phi) = Integrate(0 to phi){PDF(phi) dphi}.                             (A5)

 However, for the cases intended to be solved with this code (macroscopically
 isotropic medium, particle ensemble with random orientation), scattering has no
 azimuthal dependency, such that those equations simplify to:

 PDF(psi) = 2 PI B(psi,phi) sin(psi),                                       (A6)

 CDF(psi) = 2 PI Integrate(0 to psi){dpsi B(psi,phi) sin(psi)}              (A7)

 For completeness, however, the functions receive or output psi and phi. The 
 mathematical development for phi in absence of azimuthal dependency is
 presented in the sources.c file.

 In the current version, polarization is not implemented.

 IMPLEMENTATION:
 The implementation of these functions follow the paradigm of dynamic
 dispatching, as used throughout the program. This flexibility might result in
 penalty performances due to additional function pointer deferencing and unknown
 functions at compile time preventing inlining. However, limited execution time
 tests showed no loss in performance against direct call of specific functions.

 The group of functions for each scattering type are:
 scat_qtl_x_y - Scattering quantile function.
 scat_cdf_x_y - Scattering cumulative distribution function;
 scat_pdf_x_y - Scattering probability distribution function;
 scat_pf_x_y  - Scattering phase function;

 where x can be:
 iso - Isotropic;
 ray - Rayleigh (air or water);
 hg  - Henyey-Greestein (HG);
 ptz - Petzold (average);
 ff  - Fournier-Forand (FF);

 and y can be:
 min  - minimization;
 lut  - look-up table direct access;
 intp - look-up table interpolation;

 For long look-up tables, interpolation is more efficient, unless direct access
 is used with constant data intervals with findBin function.
 
 It is possible to sample a random direction by using the qtl function with
 random uniform (0,1] deviates for probabilities of psi and phi.

*******************************************************************************/

 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <math.h>

 #include <gsl/gsl_errno.h>
 #include <gsl/gsl_math.h>
 #include <gsl/gsl_min.h>
 #include <gsl/gsl_spline.h>
 #include <gsl/gsl_integration.h>

 #include "config.h"
 #include "constants.h"
 #include "scattering.h"

/* Setup functions: ***********************************************************/


/* scat_alloc, scat_free

 Allocates or free memory to/from a scattering struct.

 INPUT:
 scat - Scattering parameters
        Pointer to pointer to scattering struct;

 OUTPUT:
 Pointer to a scattering struct or none when freeing memory.

 */

 struct scattering*
 scat_alloc
 ( void )
 {
   struct scattering* scat = (struct scattering*) 
     malloc(sizeof(struct scattering));
   #ifdef CHCK_MEM_ALLOC
   if ( !scat )
   {
     printf("\nERROR: Failed to allocate memory for scat in scat_alloc "
       "(scattering.c)\n");
     exit(-1);
   }
   #endif // CHCK_MEM_ALLOC

   return scat;
 }

 void
 scat_free
 ( struct scattering **scat )
 {
   free( *scat );
   *scat = NULL;
 }

/* scat_lut_alloc, scat_lut_free

 Allocate or free memory for a scatering look-up table struct.

 INPUT:
 scat_lut - Scattering LUT
            Pointer to pointer to scattering_lut struct.

 OUTPUT:
 Pointer to scattering_lut struct or none when freeing memory.

 */

 struct scattering_lut*
 scat_lut_alloc
 ( void )
 { return (struct scattering_lut*) malloc(sizeof(struct scattering_lut)); }

 void
 scat_lut_free
 ( struct scattering_lut **scat_lut )
 {
   if ( (*scat_lut)->psi ) free((*scat_lut)->psi);
   if ( (*scat_lut)->pf )  free((*scat_lut)->pf);
   if ( (*scat_lut)->cdf ) free((*scat_lut)->cdf);
   free(*scat_lut);
   *scat_lut = NULL;
 }

/* scat_itp_alloc, scat_itp_free

 Allocate or free memory for a scatering interpolator struct.

 INPUT:
 scat_itp - Scattering interpolator
            Pointer to pointer to scattering_itp struct.

 OUTPUT:
 Pointer to scattering_itp struct or none when freeing memory.

 */

 struct scattering_itp*
 scat_itp_alloc
 ( void )
 { return (struct scattering_itp*) malloc(sizeof(struct scattering_itp)); }

 void
 scat_itp_free
 ( struct scattering_itp **scat_itp )
 {
   gsl_interp_accel_free( (*scat_itp)->accl_qtl );
   gsl_interp_accel_free( (*scat_itp)->accl_cdf );
   gsl_interp_accel_free( (*scat_itp)->accl_pf  );

   gsl_spline_free( (*scat_itp)->itpl_qtl );
   gsl_spline_free( (*scat_itp)->itpl_cdf );
   gsl_spline_free( (*scat_itp)->itpl_pf );

   free(*scat_itp);
   *scat_itp = NULL;
 }

/* scat_min_alloc, scat_min_free

 Allocate or free memory for a scatering minimization struct.

 INPUT:
 scat_min - Scattering minimization
            Pointer to pointer to scattering_min struct.

 OUTPUT:
 Pointer to scattering_min struct or none when freeing memory.

 */

 struct scattering_min*
 scat_min_alloc
 ( void )
 { return (struct scattering_min*) malloc(sizeof(struct scattering_min)); }

 void
 scat_min_free
 ( struct scattering_min **scat_min )
 {
   // NOT YET IMPLEMENTED
 }


/* scat_setup

 Setups the parameters of scattering.

 INPUT:
 scat     - Scattering parameters
            Pointer to scattering struct;
 scat_tp  - Scattering function name: "isotropic", "rayleigh", "hg", "ff", 
            "petzold"
            Pointer to constant char;
 depolr   - Depolarization ratio for Rayleigh scattering
            Constant double;
 g        - Pre-defined average cosine for the HG function
            Constant double;
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
 scat_setup
 (
   struct scattering *scat,
   char const   *scat_tp,
   double const depolr,
   double const g,
   double const fbb,
   char const   *method,
   int const    truncate
 )
 {
   strncpy(scat->tp, scat_tp, STRMXLEN);
   if ( strcmp (scat_tp, "isotropic") == 0 )
   {
     // Isotropic:
     scat_setup_iso(scat);
     return;
   }
   else if ( strcmp(scat_tp, "rayleigh") == 0 )
   {
     // Rayleigh:
     scat_setup_rlg(scat, depolr);
     return;
   }
   else if ( strcmp(scat_tp, "hg") == 0 )
   {
     // Henyey-Greestein:
     scat_setup_hg(scat, g);
     return;
   }
   else if ( strcmp(scat_tp, "ff") == 0 )
   {
     // Fournier-Forand:
     scat_setup_ff(scat, fbb, truncate, method);
     return;
   }
   else if ( strcmp(scat_tp, "ptz") == 0 )
   {
     // Petzold average:
//     scat_setup_ptz(scat, method, truncate);
     return;
   }
 }

/* scat_printf

 Prints scattering parameters.

 INPUT:
 odv    - Output device. Should be a file stream (e.g. 'stdout' to print to
          default);
 scat   - Scattering parameters
          Pointer to pointer to constant scattering struct;
 indent - Number of indentations at the header level (each indent = 2 spaces)
          Constant int;

 OUTPUT:
 None. Prints the parameters of the scattering function.

 */

 void
 scat_fprintf
 ( 
   FILE *odv,
   struct scattering const *scat,
   int const indent 
 )
 {
   // Set indentation:
   char pre_0[STRMXLEN] = "";
   char pre_1[STRMXLEN] = "";
   for (size_t i = 0; i < indent; i++)
     strcat(pre_0, "  ");
   strncpy(pre_1, pre_0, STRMXLEN);
   strcat(pre_1, "  ");

   char scat_lnm[STRMXLEN];
   if ( strcmp (scat->tp, "isotropic") == 0 )
     strncpy(scat_lnm, "Isotropic (Henyey-Greenstein with g = 0)", STRMXLEN);
   if ( strcmp (scat->tp, "rayleigh") == 0 )
     strncpy(scat_lnm, "Rayleigh", STRMXLEN);
   if ( strcmp (scat->tp, "hg") == 0 )
     strncpy(scat_lnm, "Henyey-Greestein", STRMXLEN);
   if ( strcmp (scat->tp, "ff") == 0 )
     strncpy(scat_lnm, "Fournier-Forand", STRMXLEN);
   if ( strcmp (scat->tp, "petzold") == 0 )
     strncpy(scat_lnm, "Petzold (average)", STRMXLEN);

   char mtd_lnm[STRMXLEN];
   if ( strcmp (scat->method, "analytical") == 0 )
     strncpy(mtd_lnm, "Analytical", STRMXLEN);
   if ( strcmp (scat->method, "itp") == 0 )
     strncpy(mtd_lnm, "Interpolation", STRMXLEN);
   if ( strcmp (scat->method, "min") == 0 )
     strncpy(mtd_lnm, "Minimization", STRMXLEN);
   if ( strcmp (scat->method, "lut") == 0 )
     strncpy(mtd_lnm, "LUT direct search", STRMXLEN);

   fprintf(odv, "%sScattering parameters:\n", pre_0);
   fprintf(odv, "%sType:       %s\n", pre_1, scat_lnm);
   fprintf(odv, "%sfbb:        %.4lf\n", pre_1, scat->fbb);
   fprintf(odv, "%sAverage mu: %.4lf\n", pre_1, scat->g);
   fprintf(odv, "%sTruncate:   %d\n", pre_1, scat->truncate);
   fprintf(odv, "%sMethod:     %s\n", pre_1, mtd_lnm);
 }


/* Utility functions: *********************************************************/

 static double
 scat_av_mu_f
 (
   double x, 
   void * params
 )
 {
   struct scattering *scat = (struct scattering*) params;
   double s_scat[2] = {0.0};
   double dens[2] = {0.0};

   s_scat[0] = x;
   scat->fun.pdf(scat, s_scat, dens);
   return cos(x) * dens[0];
 }

 double
 scat_av_mu
 ( struct scattering * scat )
 {
   gsl_integration_workspace * w
     = gsl_integration_workspace_alloc (1000);

   double res, error;

   gsl_function F;
   F.function = &scat_av_mu_f;
   F.params   = scat;

   gsl_integration_qag (&F, 0, M_PI, 1e-4, 0, 1000, 6, w, &res, &error);
   gsl_integration_workspace_free (w);

   return res;
 }


