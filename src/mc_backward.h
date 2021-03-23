
 #ifndef MC
 #define MC

 #include "sources.h"
 #include "skyrad.h"
 #include "accumulators.h"
 #include "scattering.h"
 #include "bottom.h"
 #include "structures.h"

 int
 bmc
 (
   int const sim_nr, 
   int const sim_ns, 
   int const iop_nw0,
   int const btt_nbr,
   double const *sim_sza, 
   double const *sim_saa, 
   double const sim_f0,
   double CMPLX_T const iop_na, 
   double CMPLX_T const iop_nw, 
   double const iop_c, 
   double const *iop_w0,
   int const str_ncl,
   int const str_ncn,
   int const str_ncb,
   struct source const * src,
   struct scattering const * scat,
   struct bottom const * btt,
   struct skyradiance const * skr,
   str_cyln const ** cylns,
   struct str_cone const ** cones,
   str_cubd const ** cubds,
   struct accumulator_bmc * accm_df_f,
   struct accumulator_bmc * accm_df_s,
   struct accumulator_bmc * accm_dr_f,
   struct accumulator_bmc * accm_dr_s
 );

 #endif // MC

