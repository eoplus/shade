
 #ifndef INTERSECT_CUBOID
 #define INTERSECT_CUBOID

 #include "structures.h"	// str_cubd

/* Function prototypes: *******************************************************/

 int
 intrs_cubd
 (
   double const * a,
   double const * b,
   double const * u,
   str_cubd const * cubd
 );

 int
 intrs_cubd_closed
 (
   double const * a,
   double const * b,
   double const * u,
   str_cubd const * cubd
 );

 #endif // INTERSECT_CUBOID
