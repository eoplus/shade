
 #ifndef AIRTIGHT_ALEX_FUN
 #define AIRTIGHT_ALEX_FUN

 #include "ray.h"
 #include "structures.h"
 #include "sources.h"
 #include "tr_interface.h"

 #include <gsl/gsl_rng.h>

 void
 check_interface_up_airtight_alex
 (
   struct light_ray * ray,
   double const l,
   str_cyln const * cyln,
   struct source const * src,
   double const muc, 
   double const *n_rat_wa,
   double const iop_na,
   double const iop_nw,
   struct Fm Fmat,
   int const btt_nbr,
   int const iop_nw0,
   double const * ref_o_1,
   double const * ref_o_2,
   double const ** M_open,
   gsl_rng * random
 );

 void
 check_interface_down_airtight_alex
 (
   struct light_ray * ray,
   double const l,
   str_cyln const * cyln,
   struct source const * src,
   double const muc, 
   double const *n_rat_wa,
   double const iop_na,
   double const iop_nw,
   struct Fm Fmat,
   int const btt_nbr,
   int const iop_nw0,
   double const * ref_o_1,
   double const * ref_o_2,
   double const ** M_open,
   gsl_rng * random
 );

 #endif // AIRTIGHT_ALEX_FUN

