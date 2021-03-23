
 #ifndef INTERSECT
 #define INTERSECT

 #include "structures.h"
 #include "rotation.h"
 #include "intersect_cylinder.h"
 #include "intersect_cuboid.h"
 #include "intersect_cone.h"

/* Function skeletons: ********************************************************/

 int
 chk_intrs
 (
   double const * a,
   double const * b_i,
   double const * u,
   double const l_i,
   int const NCL,
   str_cyln const ** cylns,
   int const NCN,
   struct str_cone const ** cones,
   int const NCB,
   str_cubd const ** cubds
 );

 // Plane:
 void
 intrs_plane
 (
   int * intrs_f,
   double * intrs_l,
   double * intrs_p,
   double const * v,
   double const * u,
   double const l
 );

 // Rectangles:
 int
 intrs_rect
 (
   double const * v, 
   double const * u,
   double const l,
   str_rect const * rect
 );

 // Ellipses:
 int
 intrs_ellp
 (
   double const * v,
   double const * u, 
   double const l,
   str_ellp const * ellp
 );

 #endif // INTERSECT

