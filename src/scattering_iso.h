
 #ifndef SCATTERING_ISO
 #define SCATTERING_ISO

 #include "scattering.h"

/* Function prototypes: *******************************************************/

 void
 scat_setup_iso
 ( struct scattering * scat );

 void
 scat_qtl_iso
 (
   struct scattering const * scat,
   double * s_scat,
   double const * prob 
 );

 void
 scat_cdf_iso
 (
   struct scattering const * scat,
   double const * s_scat,
   double * prob
 );

 void
 scat_pdf_iso
 (
   struct scattering const * scat,
   double const * s_scat,
   double * dens
 );

 void
 scat_pf_iso
 (
   struct scattering const * scat,
   double const * s_scat,
   double * pf
 );

 #endif // SCATTERING_ISO
