
 #ifndef SCATTERING_HG
 #define SCATTERING_HG

 #include "scattering.h"

/* Function prototypes: *******************************************************/

 void
 scat_setup_hg
 ( struct scattering * scat,
   double const g
 );

 void
 scat_qtl_hg
 (
   struct scattering const * scat,
   double * s_scat,
   double const * prob 
 );

 void
 scat_cdf_hg
 (
   struct scattering const * scat,
   double const * s_scat,
   double * prob
 );

 void
 scat_pdf_hg
 (
   struct scattering const * scat,
   double const * s_scat,
   double * dens
 );

 void
 scat_pf_hg
 (
   struct scattering const * scat,
   double const * s_scat,
   double * pf
 );

 #endif // SCATTERING_HG
