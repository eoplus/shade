
 #ifndef SCATTERING_RLG
 #define SCATTERING_RLG

 #include "scattering.h"

/* Function prototypes: *******************************************************/

 void
 scat_setup_rlg
 (
   struct scattering * scat,
   double const depolr
 );

 // Generic Rayleigh:
 void
 scat_qtl_rlg
 (
   struct scattering const * scat,
   double * s_scat,
   double const * prob
 );

 void
 scat_cdf_rlg
 (
   struct scattering const * scat,
   double const * s_scat,
   double * prob
 );

 void
 scat_pdf_rlg
 (
   struct scattering const * scat,
   double const * s_scat,
   double * dens
 );

 void
 scat_pf_rlg
 (
   struct scattering const * scat,
   double const * s_scat,
   double * pf
 );

 // Rayleigh for standard air:
 void
 scat_qtl_rlg_air
 (
   struct scattering const * scat,
   double * s_scat,
   double const * prob
 );

 void
 scat_cdf_rlg_air
 (
   struct scattering const * scat,
   double const * s_scat,
   double * prob
 );

 void
 scat_pdf_rlg_air
 (
   struct scattering const * scat,
   double const * s_scat,
   double * dens
 );

 void
 scat_pf_rlg_air
 (
   struct scattering const * scat,
   double const * s_scat,
   double * pf
 );

 #endif // SCATTERING_RLG

