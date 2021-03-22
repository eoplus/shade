
/*******************************************************************************
 intersect_cuboid.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 Version: 2021-03-12

 intrs_cubd, intrs_cubd_closed

 This source file provides functions to check for intersection of a ray against
 a cuboid.

 The specific functions assume that the ray's initial and terminal points are 
 centered on the cuboid's origin and that if the cuboid's orientation axis
 does not have the same direction as the system's Z-axis, the ray points and 
 direction were rotated to the cuboid's reference frame. Those transformations
 are performed by the wrapper function intrs_cubd which also calls the
 appropriate function depending on cuboid type.

 In version 1.6, only one cuboid type is available: right angles and closed at
 top and bottom. Future versions may include parallelepipeds and open cuboids.

 The function below implements a variant of the Williams et al. (2005)
 algorithm.

 REFERENCES:
 Williams, A.; Barrus, S.; Morley, R. K.; Shirley, P. 2005. An efficient and 
   robust ray-box intersection algorithm. In ACM SIGGRAPH 10, 9.
   DOI:10.1145/1198555.1198748

 INPUTS:
 a    - Initial point of the ray
        Dimensions: [3]: [0]X, [1]Y, [2]Z
        Range: (-Inf,Inf), meters
        Pointer to constant double;
 b    - Terminal point of the ray; 
        Dimensions: [3]: [0]X, [1]Y, [2]Z
        Range: (-Inf,Inf), meters
        Pointer to constant double;
 u    - Cartesian direction of the ray; 
        Dimensions: [3]: [0]X, [1]Y, [2]Z
        Range: [-1,1], unitless
        Pointer to constant double;
 cyln - Cuboid parameters
        Pointer to str_cubd.
 
 OUTPUTS:
 Flag indicating if an intersection occurred: (0) No, (1) Yes.

*******************************************************************************/

 #include <stddef.h>		// size_t
 #include <math.h>		// sqrt

 #include "rotation.h"		// ROT_VEC, ROT_VEC_UNIT
 #include "structures.h"	// str_cubd
 #include "intersect.h"

 int
 intrs_cubd
 (
   double const * a,
   double const * b,
   double const * u,
   str_cubd const * cubd
 )
 {
   double a_cr[3];	// Ray's initial point, centered and rotated
   double b_cr[3];	// Ray's terminal point, centered and rotated 
   double u_rot[3];	// Ray's Cartesian directions, rotated

   if ( cubd->rotate_f[0] )
   {
     double c_vec[3];	// Temporary variable for the macro ROT_VEC
     ROT_VEC(a_cr, a, cubd->o, cubd->M, c_vec);
     ROT_VEC(b_cr, b, cubd->o, cubd->M, c_vec);
     ROT_VEC_UNIT(u_rot, u, cubd->M);
   }
   else
   {
     for (size_t i = 0; i < 3; i++)
     {
       a_cr[i]  = a[i] - cubd->o[i];
       b_cr[i]  = b[i] - cubd->o[i];
       u_rot[i] = u[i];
     }
   }

   // In the current version, open cuboids and more general paralelepipeds are
   // not implemented.
   return intrs_cubd_closed(a_cr, b_cr, u_rot, cubd);
 }

 int
 intrs_cubd_closed
 (
   double const * a,
   double const * b,
   double const * u,
   str_cubd const * cubd
 )
 {
   double u_inv[3];
   double dmin[3];
   double dmax[3];

   for (size_t i = 0; i < 3; i++)
   {
     u_inv[i] = 1.0 / u[i];
     if (u_inv[i] >= 0)
     {
       if (b[i] < -cubd->hl[i] || a[i] >= cubd->hl[i])
       {
         return 0;
       }
       dmin[i] = (-cubd->hl[i] - a[i]) * u_inv[i];
       dmax[i] = ( cubd->hl[i] - a[i]) * u_inv[i];
     } 
     else
     {
       if (a[i] <= -cubd->hl[i] || b[i] > cubd->hl[i])
       {
         return 0;
       }
       dmin[i] = ( cubd->hl[i] - a[i]) * u_inv[i];
       dmax[i] = (-cubd->hl[i] - a[i]) * u_inv[i];
     }
   }
   if ( (dmin[0] > dmax[1]) || (dmin[1] > dmax[0]) )
   {
     return 0;
   } else {
     if (dmin[1] > dmin[0]) dmin[0] = dmin[1];
     if (dmax[1] < dmax[0]) dmax[0] = dmax[1];
     if ( (dmin[0] > dmax[2]) || (dmin[2] > dmax[0]) )
     {
       return 0;
     } else {
       return 1;
     }
   }
 }

