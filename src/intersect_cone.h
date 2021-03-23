
 #ifndef INTERSECT_CONE
 #define INTERSECT_CONE

 #include "structures.h"	// str_cone

/* Function prototypes: *******************************************************/

 int
 intrs_cone
 (
   double const * a, 
   double const * b, 
   double const * u, 
   double const l, 
   struct str_cone const * cone
 );

 int
 intrs_cone_right
 (
   double const * a, 
   double const * b, 
   double const * u, 
   double const l, 
   struct str_cone const * cone
 );

 int
 intrs_cone_oblique
 (
   double const * a, 
   double const * b, 
   double const * u, 
   double const l, 
   struct str_cone const * cone
 );

 #endif // INTERSECT_CONE
