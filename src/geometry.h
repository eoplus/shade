
 #ifndef GEOMETRY
 #define GEOMETRY

 #include <math.h>
 #include <gsl/gsl_rng.h>

 #include "config.h"
 #include "constants.h"
 #include "aux.h"

/* Spherical to Cartesian coordinates: *****************************************

 sph_to_cos, sph_to_cos_unit

 Convert from spherical directions to Cartesian directions (three dimensions). 
 A vector with arbitrary radius in spherical coordinates can be converted with 
 the function sph_to_cos. If the vector is a unit vector, or if the radius 
 information is not relevant, the function sph_to_cos_unit should be used as it 
 is more efficient.

 The functions are also available as macros. Though the static inline
 declaration, together with the __attribute__((always_inline)) qualifier should
 result that the functions are always inlined and effectivelly behaive as macros
 my tests showed that the macros were 3 times faster.

 INPUT:
 u - Pointer to the array receiving the Cartesian coordinates.
 s - Pointer to the array containing the spherical coordinates.

 OUTPUT:
 None. Updates the values pointed by u.

*******************************************************************************/

 static inline __attribute__((always_inline)) void
 sph_to_cos
 (
   double * u, 
   double const * s
 )
 {
   double sin_t = sin(s[0]);
   u[0] = sin_t * cos(s[1]) * s[2];
   u[1] = sin_t * sin(s[1]) * s[2];
   u[2] = cos(s[0]) * s[2];
 }

 #define SPH_TO_COS(I_u, I_s, I_sin_t) \
         (I_sin_t) = sin( (I_s)[0] ); \
         (I_u)[0]  = (I_sin_t) * cos( (I_s)[1] ) * (I_s)[2]; \
         (I_u)[1]  = (I_sin_t) * sin( (I_s)[1] ) * (I_s)[2]; \
         (I_u)[2]  = cos( (I_s)[0] ) * (I_s)[2];

 static inline __attribute__((always_inline)) void 
 sph_to_cos_unit
 (
   double * u,
   double const * s
 )
 {
   double sin_t = sin(s[0]);
   u[0] = sin_t * cos(s[1]);
   u[1] = sin_t * sin(s[1]);
   u[2] = cos(s[0]);   
 }

 #define SPH_TO_COS_UNIT(I_u, I_s, I_sin_t) \
         (I_sin_t) = sin( (I_s)[0] ); \
         (I_u)[0]  = (I_sin_t) * cos( (I_s)[1] ); \
         (I_u)[1]  = (I_sin_t) * sin( (I_s)[1] ); \
         (I_u)[2]  = cos( (I_s)[0] );

/* Cartesian to spherical coordinates: *****************************************

 cos_to_sph, cos_to_sph_unit

 Calculate the spherical coordinates for a given Cartesian coordinates. If the 
 Cartesian coordinates form a unit vector, the function cos_to_sph_unit can be 
 used and is more efficient.

 If the ray is traveling close to the Z-axis (u[Z] ~ 1.0 or -1.0), than is safer 
 to treat as if it is on the Z-axis and return a random azimuth, as the azimuth
 in this case is undefined. For this purpose, the simple rand function from the
 standard library is used.

 INPUT:
 s    - Pointer to the GSL vector to receive the spherical coordinates.
 u    - Pointer to the GSL vector containing the Cartesian coordinates.
 rand - Pointer to a GSL random number generator.

 OUTPUT:
 Updates the values pointed by s.

*******************************************************************************/

 static inline __attribute__((always_inline)) void
 cos_to_sph
 (
   double * s, 
   double const * u,
   gsl_rng *rand
 )
 {
   s[2] = sqrt( u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
   s[0] = acos(u[2] / s[2]);

   if ( (s[2]- ABS(u[2])) < TOLERANCE )
   { 
     s[1] = K_2PI * gsl_rng_uniform(rand);
   } 
   else
   {
    s[1] = atan2(u[1], u[0]);
    if ( u[1] < 0.0 )
      s[1] += K_2PI;
   }
 }

 #define COS_TO_SPH(I_s, I_u, I_random) \
         (I_s)[2] = sqrt( (I_u)[0] * (I_u)[0] + (I_u)[1] * (I_u)[1] + \
           (I_u)[2] * (I_u)[2]); \
         (I_s)[0] = acos( (I_u)[2] ); \
         if( ABS( (I_s)[2] - (I_u)[2] ) < TOLERANCE ) \
         { \
           (I_s)[1] = K_2PI * gsl_rng_uniform (I_random); \
         } else { \
            (I_s)[1] = atan2( (I_u)[1], (I_u)[0] ); \
            if ( (I_u)[1] < 0.0 ) \
            { \
              (I_s)[1] += K_2PI; \
            } \
         }

 static inline __attribute__((always_inline)) void 
 cos_to_sph_unit
 (
   double * s, 
   double const * u,
   gsl_rng *rand
 )
 {
   s[2] = 1.0;
   s[0] = acos(u[2]);

   if ( (1.0 - ABS(u[2])) < TOLERANCE )
   { 
     s[1] = K_2PI * gsl_rng_uniform(rand);
   } 
   else
   {
    s[1] = atan2(u[1], u[0]);
    if ( u[1] < 0.0 )
      s[1] += K_2PI;
   }
 }

 #define COS_TO_SPH_UNIT(I_s, I_u, I_random) \
         (I_s)[2] = 1.0; \
         (I_s)[0] = acos( (I_u)[2] ); \
         if( ABS( 1.0 - (I_u)[2] ) < TOLERANCE ) \
         { \
           (I_s)[1] = K_2PI * gsl_rng_uniform (random); \
         } else { \
            (I_s)[1] = atan2( (I_u)[1], (I_u)[0] ); \
            if ( (I_u)[1] < 0.0 ) \
            { \
              (I_s)[1] += K_2PI; \
            } \
         }

/* Matrix and vector multiplication: *******************************************

 mat_vec_prod

 Performs the multiplication of a vector v by matrix M: Mv. Note that no range 
 checking is performed - sizes of the array must be correct.

 INPUT:
 w  - Pointer to array receiving the result. Must have nr elelments;
 v  - Pointer to array containing the column vector. Must have nc elelments;
 M  - Pointer to array containing the matrix;
 nr - Number of rows of M;
 nc - Number of columns of M.

 OUTPUT:
 None. Updates the values pointed by w.

*******************************************************************************/
 
 static inline __attribute__((always_inline)) void
 mat_vec_prod 
 (
   double * ww,
   double const * vv,
   double const ** M,
   int const nr,
   int const nc
 )
 {
   for (int cr = 0; cr < nr; cr++)
   {
     ww[cr] = 0.0;
     for (int cc = 0; cc < nc; cc++)
     {
       ww[cr] += M[cr][cc] * vv[cc];
     }
   }
 }

 #define MAT_VEC_PROD(I_w, I_v, I_M, I_nr, I_nc) \
         for (int cr = 0; cr < (I_nr); cr++) \
         { \
           (I_w)[cr] = 0.0; \
           for (int cc = 0; cc < (I_nc); cc++) \
           { \
             (I_w)[cr] += (I_M)[cr][cc] * (I_v)[cc]; \
           } \
         }

/* Function prototypes: *******************************************************/

 void
 mat_inv
 (
   double ** N, 
   double const ** M,
   int const nr,
   int const nc
 );

 double
 mat_det 
 (
   double const ** M,
   int const nr,
   int const nc
 );

 void 
 mat_transpose
 (
   double ** N,
   double const ** M,
   int const nr,
   int const nc
 );

 int
 mat_equal
 (
   double const ** M, 
   double const ** N,
   int const nr,
   int const nc, 
   double const TOL
 );

 void
 vec_fprintf
 (
   FILE * odv,
   double const * v,
   int const n ,
   double const scl, 
   char const * title,
   int const indent
 );

 void
 mat_fprintf
 (
   FILE * odv,
   double const ** M,
   int const nr,
   int const nc, 
   double const scl, 
   char const * title,
   int const indent
 );

 #endif // GEOMETRY

