
 #ifndef ROTATION
 #define ROTATION

 #include "geometry.h"

/* Rotation matrix ZY0: ********************************************************

 rot_mat_ZY0

 Calculates the rotation matrix based on sequential Z and Y rotations. This is 
 the same as the rot_mat_ZYX with the last rotation being 0º, however, it is 
 simpler and uses the directional cosines as inputs. This is the generating 
 function for the rotation matrix applied to scattering angles in the MC code.

 Since the polar angles of the Z-axis are expected to be known, and the 
 directional cosines are expected to be available, the function uses the 
 directional cosines directly.
 
 INPUT:
 u - Pointer to array containing the directional cosines of the Z-axis;
 M - Pointer to pointer to array (matrix) to receive the rotation matrix.

 OUTPUT:
 Updates the values of M.

*******************************************************************************/

 static inline __attribute__((always_inline)) void 
 rot_mat_ZY0 
 (
   double ** M,
   double const * u
 )
 {
   double const sin_t = sqrt(1.0 - u[2] * u[2]);

   M[0][0] =  u[0] * u[2] / sin_t;		//  cos(p)cos(t)
   M[0][1] = -u[1] / sin_t;			// -sin(p)
   M[0][2] =  u[0]; 				//  cos(p)sin(t)
   M[1][0] =  u[1] * u[2] / sin_t;		//  sin(p)cos(t)
   M[1][1] =  u[0] / sin_t;			//  cos(p)
   M[1][2] =  u[1];				//  sin(p)sin(t)
   M[2][0] = -sin_t;				// -sin(t)
   M[2][1] =  0.0;				//  0
   M[2][2] =  u[2]; 				//  cos(t)
 }

 #define ROT_MAT_ZY0(I_M, I_u, I_sin_t) \
         (I_sin_t) = sqrt(1.0 - (I_u)[2] * (I_u)[2]); \
         (I_M)[0][0] =  (I_u)[0] * (I_u)[2] / (I_sin_t); \
         (I_M)[0][1] = -(I_u)[1] / (I_sin_t); \
         (I_M)[0][2] =  (I_u)[0]; \
         (I_M)[1][0] =  (I_u)[1] * (I_u)[2] / (I_sin_t); \
         (I_M)[1][1] =  (I_u)[0] / (I_sin_t); \
         (I_M)[1][2] =  (I_u)[1]; \
         (I_M)[2][0] = -(I_sin_t); \
         (I_M)[2][1] =  0.0; \
         (I_M)[2][2] =  (I_u)[2];

/* Rotate vector: **************************************************************

 rot_vec, rot_vec_unit

 Applies centering and rotation to a bounded vector (ray, centered to origin 
 before rotation) or Cartesian directions (unit vector, no centering).

 INPUT:
 w - Pointer to array receiving the rotated coordinates.
 v - Pointer to array with the original coordinates.
 o - Pointer to array with the center of rotation.
 M - Pointer to pointer to array (matrix) with the rotation matrix.

 OUTPUT:
 None. Updates the values pointed by w.

*******************************************************************************/

 static inline __attribute__((always_inline)) void
 rot_vec
 (
   double * w,
   double const * v,
   double const * o,
   double const ** M
 )
 {
   double c[3];
   for (size_t i = 0; i < 3; i++)
   {
     c[i] = v[i] - o[i];
   }
   mat_vec_prod(w, c, M, 3, 3);
 }

 #define ROT_VEC(I_w, I_v, I_o, I_M, I_c) \
         for (size_t i = 0; i < 3; i++) \
         { \
           (I_c)[i] = (I_v)[i] - (I_o)[i]; \
         } \
         MAT_VEC_PROD( (I_w), (I_c), (I_M), 3, 3 );

 static inline __attribute__((always_inline)) void
 rot_vec_unit
 (
   double * w,
   double const * v,
   double const ** M
 )
 {
   mat_vec_prod(w, v, M, 3, 3);
 }

 #define ROT_VEC_UNIT(I_w, I_v, I_M) \
         MAT_VEC_PROD( (I_w), (I_v), (I_M), 3, 3 );

/* Function prototypes: *******************************************************/

 void rot_mat_ZYX (double **M, double const *s);

 void rot_mat_ZYZ (double **M, double const *u, double const *s, 
   double const alpha);

 #endif // ROTATION

