
 #ifndef INTERSECT
 #define INTERSECT

 #include "structures.h"
 #include "rotation.h"

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
   str_cone const ** cones,
   int const NCB,
   str_cubd const ** cubds
 );

 // Planes:
 void intrs_plane (int *intrs_f, double *intrs_l, double *intrs_p,
   double const *v, double const *u, double const l);

 // Rectangles:
 int intrs_rect (double const *v, double const *u, double const l,
   str_rect const *rect);

 // Ellipses:
 int intrs_ellp (double const *v, double const *u, double const l,
   str_ellp const *ellp);

 // Cuboids:
 int intrs_cubd (double const *a, double const *b, double const *u,
   str_cubd const *cubd);

 int intrs_cubd_closed (double const *a, double const *b, double const *u,
   str_cubd const *cubd);

 // Cylinders:
 int intrs_cyln (double const *a, double const *b, double const *u, 
  double const l, str_cyln const *cyln);

 int intrs_cyln_right (double *a, double const *b, double const *u, 
   double const l, str_cyln const *cyln);

 int intrs_cyln_oblique (double const *a, double const *b, double const *u, 
   double const l, str_cyln const *cyln);

 #endif // INTERSECT
