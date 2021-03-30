
 #ifndef BOTTOM
 #define BOTTOM

 #include "config.h"
 #include "reflectance.h"

/*******************************************************************************
 bottom

 In future versions, the idea is to expand the bottom to include a triangle mesh 
 for spatial representation, allowing for spattially resolved orientation,
 depth, and reflectance (at least a few reflectance types/values, associated to
 each triangle by a coding class - refl in this case would be allocated in an
 array). This class may be redundant. For example, the air-water interface is 
 also represented by a triangle mesh and complex surfaces could be likewise.
 however, might be useful to have specific class for an horizonthal surface that 
 might extend to infinity.
 
*******************************************************************************/

 struct bottom
 {
   double depth;
   struct reflectance refl;
 };

/* Function prototypes: *******************************************************/

 struct bottom *
 btt_alloc
 ( void );

 void
 btt_free
 ( struct bottom ** btt );

 void
 btt_setup
 (
   struct bottom * btt,
   double const depth,
   char const * type,
   int const nbhr,
   double const * bhr,
   double const k
 );

 void
 btt_fprintf
 (
   FILE * odv,
   struct bottom const * btt,
   int const indent 
 );

 #endif // BOTTOM

