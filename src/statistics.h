
 #ifndef STATISTICS
 #define STATISTICS

 #include "accumulators.h"

 void
 calc_se
 (
   struct accumulator_bmc *accm_se,
   struct accumulator_bmc *accm_mn,
   double *btt_bhr,
   char   *sufx,
   char   *ofbn
 );

 #endif // STATISTICS
