
 #ifndef SCATTERING
 #define SCATTERING

 #include <gsl/gsl_errno.h>
 #include <gsl/gsl_math.h>
 #include <gsl/gsl_min.h>
 #include <gsl/gsl_spline.h>

 #include "config.h"

/*******************************************************************************
 scattering

 These structures hold information necessary to retrieve properties from
 scattering functions.

 Due to the diversity of scattering functions and their parameterizations,
 the scattering struct contain pointers to strcuts specific to each scattering
 type. The types currently implemented are:

 (1) scattering_ray - Rayleigh scattering; (scattering_ray.c source file)
 (2) scattering_hg  - Henyey-Greestein; (scattering_hg.c source file)
 (3) scattering_ff  - Fournier-Forand; (scattering_ff.c source file)

 For scattering functions that have no analytical solution for the quantile or 
 that have no parametrization (measured, tabulated values), different methods
 to retrieve values of the scattering functions are available:

 (1) scattering_lut - Direct search Look-up table;
 (2) scattering_itp - Interpolation of Look-up table;
 (3) scattering_min - Minimization;
 
 COMPONENTS (scattering):
 fun      - Scattering function pointers
            Pointer to scattering_fun;
 hg       - Henyey-Greestein function parameters
            Pointer to scattering_hg;
 ff       - Fournier-Forand function parameters
            Pointer to scattering_ff;
 lut      - Lookup table parameters
            Pointer to scattering_lut;
 itp      - Interpolation function parameters
            Pointer to scattering itp;
 min      - Minimization function parameters
            Pointer to scattering_min;
 pfn      - Normalization factor for the phase function
            Double;
 fbb      - Backscattering ratio of the phase function
            Double;
 g        - Average cosine of the phase function
            Double;
 truncate - Flag indicating if the diffraction peak of the PF was truncated
            Int;
 method   - Name of the solution method
            Pointer to char.

 COMPONENTS (scattering_fun):
 qtl      - Quantile function;
 cdf      - CDF function;
 pdf      - PDF function;
 pf       - PF function.

 COMPONENTS (scattering_ray):
 depolr   - Depolarization ratio
            Double;
 f        - (1 - depolr) / (1 + depolr)
            Double;
 f2       - f^2;
 f3       - f^3;
 c0       - Constant 0;
 c1       - Constant 1;
 c2       - Constant 2; 
 c3       - Constant 3; 
 c4       - Constant 4;
 c5       - Constant 5;
 c6       - Constant 6;

 COMPONENTS (scattering_hg):
 k_2g     - Two times the average cosine
            Double;
 k_g2     - Average cosine squared
            Double.

 COMPONENTS (scattering_ff):
 j        - Junge slope (particle size distribution)
            Double;
 nh       - Real part of the refractive index of hydrosols
            Double.

 COMPONENTS (scattering_lut):
 npsi     - Number of scattering angles
            Int;
 psi      - Scattering angles
            Pointer to double;
 pf       - Phase function values
            Pointer to double;
 cdf      - Cumulative distribution function values
            Pointer to double;
 kpsi     - Scaling for use in the findBin function for psi
            Double;
 kcdf     - Scaling for use in the findBin function for CDF
            Double.

 COMPONENTS (scattering_itp):
 accl_qtl - GSL interpolator accelerator for qtl; 
 accl_cdf - GSL interpolator accelerator for cdf;
 accl_pf  - GSL interpolator accelerator for pf;
 itpl_qtl - GSL linear interpolator for qtl;
 itpl_cdf - GSL linear interpolator for cdf;
 itpl_pf  - GSL linear interpolator for pf;

 COMPONENTS (scattering_min):

*******************************************************************************/

 struct scattering;
 struct scattering_fun;
 struct scattering_rlg;
 struct scattering_hg;
 struct scattering_ff;
 struct scattering_lut;
 struct scattering_itp;
 struct scattering_min;

 typedef void scat_qtl_fun (struct scattering const*, double*, double const*);
 typedef void scat_cdf_fun (struct scattering const*, double const*, double*);
 typedef void scat_pdf_fun (struct scattering const*, double const*, double*);
 typedef void scat_pf_fun  (struct scattering const*, double const*, double*);

 struct scattering_fun
 {
   scat_qtl_fun * qtl;
   scat_cdf_fun * cdf;
   scat_pdf_fun * pdf;
   scat_pf_fun  * pf;
 };

 struct scattering_rlg
 {
   double depolr;
   double f;
   double f2;
   double f3;
   double c0;
   double c1;
   double c2;
   double c3;
   double c4;
   double c5;
   double c6;
 };

 struct scattering_hg
 {
   double k_2g;
   double k_g2;
 };

 struct scattering_ff
 {
   double j;
   double nh;
 };

 struct scattering_lut
 {
   int    npsi;
   double psi[SCAT_NPSI];
   double pf[SCAT_NPSI];
   double cdf[SCAT_NPSI];
   double kpsi;
   double kcdf;
 };

 struct scattering_itp
 {
   gsl_interp_accel *accl_qtl;
   gsl_interp_accel *accl_cdf;
   gsl_interp_accel *accl_pf;
   gsl_spline       *itpl_qtl;
   gsl_spline       *itpl_cdf;
   gsl_spline       *itpl_pf;
 };

 struct scattering_min
 {
   void* a;
 };

 struct scattering
 {
   struct scattering_fun fun;
//   void * par;
//   void * mtd;
   struct scattering_rlg rlg;
   struct scattering_hg  hg;
   struct scattering_ff  ff;
   struct scattering_lut lut;
   struct scattering_itp itp;
   struct scattering_min min;
   double pfn;
   double fbb;
   double g;
   int    truncate;
   char   method[STRMXLEN];
   char   tp[STRMXLEN];
 };


/* Generic scattering functions: ***********************************************

 scat_qtl, scat_cdf, scat_pdf, scat_pf

 Generic functions to retrieve the quantile, CDF, PDF and phase function for
 a given scattering type.

 INPUT:
 scat   - Scattering parameters
          Pointer to constant scattering struct; 
 s_scat - Spherical direction
          [0]Theta, [1]Phi, [2]Radius
          Range: [0,PI] radians, [0,2PI] radians, 1.0 unitless
          Pointer to constant double or pointer to double;
 prob   - Cumulative probability associated with a polar and a azimuth angle
          [0]Theta, [1]Phi
          Range: [0,1] unitless
          Pointer to constant double or pointer to double;
 dens   - Probability density associated with a polar and a azimuth angle
          [0]Theta, [1]Phi
          Range: [0,Inf) unitless
          Pointer to double
 pf     - Phase function at a particular angle
          Range: [0,Inf), sr
          Pointer to double;

 OUTPUT:
 None. Updates values pointed by s_scat, prob or pf.

*******************************************************************************/

 static inline __attribute__((always_inline)) void
 scat_qtl
 (
   struct scattering const * scat,
   double * s_scat,
   double const * prob
 )
 { scat->fun.qtl(scat, s_scat, prob); }

 #define SCAT_QTL(I_scat, I_s_scat, I_prob) \
         scat->fun.qtl( (I_scat), (I_s_scat), (I_prob) );

 static inline __attribute__((always_inline)) void
 scat_cdf
 (
   struct scattering const * scat,
   double const * s_scat,
   double * prob
 )
 { scat->fun.cdf(scat, s_scat, prob); }

 #define SCAT_CDF(I_scat, I_s_scat, I_prob) \
         scat->fun.cdf( (I_scat), (I_s_scat), (I_prob) );

 static inline __attribute__((always_inline)) void
 scat_pdf
 (
   struct scattering const * scat,
   double const * s_scat,
   double * dens
 )
 { scat->fun.pdf(scat, s_scat, dens); }

 #define SCAT_PDF(I_scat, I_s_scat, I_dens) \
         scat->fun.pdf( (I_scat), (I_s_scat), (I_dens) );

 static inline __attribute__((always_inline)) void
 scat_pf
 (
   struct scattering const * scat,
   double const * s_scat,
   double * pf
 )
 { scat->fun.pf(scat, s_scat, pf); }

 #define SCAT_PF(I_scat, I_s_scat, I_pf) \
         scat->fun.pf( (I_scat), (I_s_scat), (I_pf) );

/* Function prototypes: *******************************************************/

 struct scattering *
 scat_alloc
 ( void );

 void
 scat_free
 ( struct scattering ** scat );

 struct scattering_lut *
 scat_lut_alloc
 ( void );

 void
 scat_lut_free
 ( struct scattering_lut ** scat_lut );

 struct scattering_itp *
 scat_itp_alloc
 ( void );

 void
 scat_itp_free
 ( struct scattering_itp ** scat_itp );

 struct scattering_min *
 scat_min_alloc
 ( void );

 void 
 scat_min_free
 ( struct scattering_min ** scat_min );

 void
 scat_fprintf
 ( 
   FILE * odv,
   struct scattering const * scat,
   int const indent 
 );

 void
 scat_setup
 (
   struct scattering * scat,
   char const * scat_tp,
   double const depolr,
   double const g,
   double const fbb,
   char const * method,
   int const truncate
 );

 // Utility functions:
 double
 scat_av_mu
 ( struct scattering * scat );

/* Specific scattering types: *************************************************/

 #include "scattering_iso.h"
 #include "scattering_rlg.h"
 #include "scattering_hg.h"
 #include "scattering_ff.h"
 #include "scattering_ptz.h"

 #endif // SCATTERING

