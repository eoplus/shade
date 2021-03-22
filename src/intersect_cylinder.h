
 #ifndef INTERSECT_CYLINDER
 #define INTERSECT_CYLINDER

 #include "structures.h"	// str_cyln

/* Function prototypes: *******************************************************/

 int
 intrs_cyln
 (
   double const * a,
   double const * b,
   double const * u,
   double const l,
   str_cyln const * cyln
 );

 int
 intrs_cyln_right
 (
   double const * a,
   double const * b,
   double const * u,
   double const l,
   str_cyln const * cyln
 );

 int
 intrs_cyln_oblique
 (
   double const * a,
   double const * b,
   double const * u,
   double const l,
   str_cyln const * cyln
 );

 #endif // INTERSECT_CYLINDER
