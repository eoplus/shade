
 #include <stdio.h>
 #include <stddef.h>
 #include <math.h>
 #include <stdlib.h>
 #include <gsl/gsl_vector.h>
 #include <gsl/gsl_matrix.h>
 #include <gsl/gsl_rng.h>

 #include "config.h"
 #include "aux.h"
 #include "geometry.h"
 #include "structures.h"
 #include "intersect.h"
 #include "intersect_cylinder.h"	// intrs_cyln
 #include "intersect_cuboid.h"		// intrs_cubd
 #include "intersect_cone.h"		// intrs_cone

/* Check intersection: *********************************************************

 chk_intrs

 This is the main function controling ray-structure intersection. It 

 Because of refraction, the intresection check is made in each refractive index 
 layer. In version 1.6 the system has only two layers, homogenous air and 
 homogenous water. The simple test included here is if the terminal point of the
 ray is in air and the initial point is in water, consider only the in water
 portion. To check in air, the initial point of the ray must be set at 0.0.

 INPUT:
 a,		// Pointer to previous ray position
 b_i,		// Pointer to current ray position
 u,		// Pointer to ray directional cosines
 l_i, 	        // Ray pathlength.
 NCL,		// Number of cylinders
 cylns,		// Pointer to cylinders
 NCN,		// Number of cones
 cones,		// Pointer to cones
 NCB,
 cubds

*******************************************************************************/

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
 )
 {
   int i;			// Counting variable
   int res;			// Intersection result
   int cubd_id;			// Start index for loop over boxes
   double b[3] = {0.0};		// Terminal point returned to surface if in air

   double l = l_i;
   if ( b_i[2] < 0.0 && a[2] > 0.0 )
   {
     // If the initial position was inside the water, intersection is caried
     // only in the water because refraction will change the angle of the ray
     // propagation. Therefore, we propagate the ray back to the surface. The
     // Z-axis position is set directly to avoid rounding errors.
     double l_vrt = b_i[2] / u[2];
     b[0] = b_i[0] + (-u[0] * l_vrt);
     b[1] = b_i[1] + (-u[1] * l_vrt);
     b[2] = 0.0;
     l -= l_vrt;
   } else {
     b[0] = b_i[0];
     b[1] = b_i[1];
     b[2] = b_i[2];
   }

   // If there are at least three cuboids, it means the first is the bounding 
   // box of all cuboids. So it is sufficient to test intersection against the 
   // bounding box to know if it is necessary to proceed to check the individual
   // cuboids. 
   cubd_id = 0;
   if(NCB > 2)
   {
     res = intrs_cubd(a, b, u, cubds[0]);
     cubd_id = ( res ) ? 1 : NCB;
   }
   for (i = cubd_id; i < NCB; i++)
   {

     res = intrs_cubd(a, b, u, cubds[i]);
     if ( res )
       return 1;
   }

   // Loop over cylinders:
   for (i = 0; i < NCL; i++)
   {
     res = intrs_cyln(a, b, u, l, cylns[i]);
     if ( res )
       return 1;
   }

   // Loop over cones:
   for (i = 0; i < NCN; i++)
   {
//     res = intrs_cone (ppos, cpos, cdir, s, &cns[i], segment);
//     if ( res )
//       return 1;
   }

   return 0;
 }

/* Intersection with a plane: **************************************************

 intrs_plane

 Checks if a bounded Euclidean vector (ray) intersects a plane. It is assumed 
 that the vector was first centered to the origin of the plane and rotated to 
 its reference frame. The function returns the intersection point if an 
 intersection occurs. This is a general function for intersections with 
 bi-dimensional structures by first testing intersection with their supporting 
 plane. Note that if the initial position of the ray lies on the plane but its
 direction is not parallel to the plane, no intersection occurs. If a ray ends 
 in a structure it intersects it. This convention is:
 intersection_l > 0.0 && intersection_l <= ray_l.

 The equation of the line is:
 P = U * l + A,                                                              (1)
 
 where l is a pathlength, U is the ray's direction (unit vector) and A is the 
 initial point of the ray. Any point in a ray must satisfy:
 ||P - A|| + ||P - B|| = ||B - A||,                                          (2)

 where || || is the Euclidean norm, and B is the terminal point of the ray.
 ||B - A|| is the magnitude of the ray and is its total pathlength, l, traveled 
 by the ray.

 Similarly, any point on the plane must satisfy:
 N · (P - O) = 0,                                                            (3)

 where N is the unit vector representing the Z axis of the plane (normal of the 
 plane), O is the origin of the plane and · is the dot (scalar) product, giving 
 the cosine of the angle between two vectors. The intersection is the point P 
 that satisfies Eq. 2 and 3. Substituting Eq. 1 into Eq. 3 and solving for l:
 N · (P - O) = 0 ->
 N · (U * l + A - O) = 0 ->
 (N · U) * l + N · (A - O) = 0 ->
 (N · U) * l = -(N · (A - O)) ->
 l = -(N · (A - O)) / (N · U).                                               (4)

 If (N · U) = 0, the vector is parallel to the plane and will only intersect the 
 plane if it lies on the plane, that is, if (N · (A - O)) = 0. If (N · U) != 0, 
 then an intersection occurs if l is positive and lower than the total 
 pathlength of the vector. The position of the intersection is found by using l
 from Eq. 4 with Eq. 1. 

 Those equations are simplified if the ray is first centered and rotated to the 
 reference frame of the plane:
 A_p = R * (A - O),                                                          (5)
 U_p = R * U,                                                                (6)
 
 where R is a rotation matrix. Then, Eq. 4 simplifies to:
 l = -A_p[Z] / U_p[Z].                                                       (7)

 If U_p[Z] = 0, the vector is parallel to the plane and will only intersect the 
 plane if it lies on the plane, that is, if A_p[Z] = 0. If U_p[Z] != 0, then an 
 intersection occurs if l is positive and lower than the total pathlength of the
 vector. And the position of intersection on the plane's reference frame is 
 found by substituting A_p and U_p in Eq. 1:
 P_p = U_p * l + A_p                                                         (8)

 INPUT:
 intrs_f - Pointer to flag receiving the result of the intersection;
 intrs_l - Pointer to array receiving pathlength to the intersection;
 intrs_p - Pointer to array receiving the position to the intersection;
 v       - Pointer to the initial point of the ray;
 u       - Pointer to the direction of the ray;
 l       - Total pathlength traveled by the ray.

 OUTPUT:
 None. Updates values pointed by intrs_f, intrs_l, intrs_p.

*******************************************************************************/

 void 
 intrs_plane
 (
   int * intrs_f,
   double * intrs_l,
   double * intrs_p,
   const double * v,
   const double * u,
   const double l
 )
 {
   if ( ABS(v[2]) < TOLERANCE )
   {
     // Ray starts on the plane
     if ( ABS(u[2]) < TOLERANCE )
     {
       // Ray's direction is parallel to the plane: intersection at the start
       *intrs_f = 1;
       *intrs_l = 0.0;
       for (size_t i = 0; i < 3; i++)
         intrs_p[i] = v[i];
       return;
     }
     *intrs_f = 0;
     return;
   }

   *intrs_l = -v[2] / u[2];
   if ( (*intrs_l) > 0.0 && (*intrs_l) <= l )
   {
     *intrs_f = 1;
     for (size_t i = 0; i < 3; i++)
       intrs_p[i] = v[i] + (*intrs_l) * u[i];
     return;
   }

   *intrs_f = 0;  
 }

/* Intersection with a rectangle: **********************************************

 intrs_rect

 Checks for intersection of a vector against an rectangle. It is assumed that 
 the vector was first centered to the origin of the plane and rotated to its 
 reference frame. The function returns the intersection point if an intersection
 occurs.

 For a point P to be inside an axis-aligned rectangle described in its 
 supporting plane, it must satisfy:
 |P[X]| < (0.5 * L[X]) && |P[Y]| < (0.5 * L[Y]),                             (1)
 
 where L is the length of each side of the rectangle.

 The intersection point of a vector with the supporting plane of the rectangle 
 is described in the function intrs_plane. With the position of intersection on 
 the plane, if an intersection with the plane occurs, the intersection with the 
 rectangle if the point satisfy Eq. 1. 

 INPUTS:
 v    - Pointer to the initial point of the ray, centered and rotated;
 u    - Pointer to the ray's direction, rotated;
 l    - Magnitude of the ray (line length from its initial to its terminal 
        point);
 rect - Pointer to an rectangle;
 
 OUTPUTS:
 A flag indicating if a intersection occurs.

*******************************************************************************/

 int
 intrs_rect
 ( 
   double const * v,
   double const * u,
   double const l,
   str_rect const * rect
 )
 {
   int    intrs_f;
   double intrs_l;
   double intrs_p[3];
   intrs_plane(&intrs_f, &intrs_l, intrs_p, v, u, l);

   if ( intrs_f )
   {
     if ( ABS(intrs_p[0]) <= (rect->hl[0] + TOLERANCE) &&
          ABS(intrs_p[1]) <= (rect->hl[1] + TOLERANCE) ) 
     {
       return 1;
     }
   }
   return 0;
 }

/* Intersection with an ellipse: ***********************************************

 intrs_ellp

 Checks for intersection of a ray against an ellipse. It is assumed that the 
 vector was first centered to the origin of the plane and rotated to its 
 reference frame. The function returns the intersection point if an intersection
 occurs.

 For a point P to be inside an ellipse described in its supporting plane, it
 must satisfy:
 1 <= || M * P ||^2 ->
 1 <= M[0][0]^2 * P[0]^2 + M[1][1]^2 * P[1]^2 ->
 1 <= P[X]^2 / S[X]^2 + P[Y]^2 / S[Y]^2,                                     (1)
 
 where M is a diagonal matrix with element [0][0] equal to the inverse of the 
 radius aligned with the X axis (S[X]), element [1][1] is the inverse of the 
 radius aligned with the Y-axis (S[Y]) and P is a position on the supporting 
 plane of the ellipse.

 If the ellipse reduces to a circle, the condition in Eq. 1 simplifies to:
 || P || < r ->
 sqrt( P · P ) < r,                                                          (2)
 
 where r is the radius and || || is the Euclidean norm. Intersections with a 
 circle can be tested with the ellipse function as circles are represented by 
 ellipse structures.

 The intersection point of a vector with the supporting plane of the ellipse 
 is described in the function intrs_plane. With the position of intersection
 on the plane, if an intersection with the plane occurs, the intersection with
 the ellipse if the point satisfy Eq. 1. 

 INPUTS:
 v    - Pointer to the initial point of the vector, centered and rotated;
 u    - Pointer to the vector direction direction, rotated;
 l    - Magnitude of the vector (line length from its initial to its terminal 
        point);
 ellp - Pointer to an ellipse;
 
 OUTPUTS:
 A flag indicating if a intersection occurs.

*******************************************************************************/

 int
 intrs_ellp
 (
   double const * v,
   double const * u,
   double const l,
   str_ellp const * ellp
 )
 {
   int    intrs_f;
   double intrs_l;
   double intrs_p[3];
   intrs_plane(&intrs_f, &intrs_l, intrs_p, v, u, l);

   if ( intrs_f )
   {
     double const rn_sq = (intrs_p[0] * intrs_p[0]) * ellp->rsq_inv[0] + 
                          (intrs_p[1] * intrs_p[1]) * ellp->rsq_inv[1];
     if ( (rn_sq - 1.0) < TOLERANCE ) 
     {
       return 1;
     }
   }
   return 0;
 }



