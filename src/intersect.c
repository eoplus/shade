
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

 int
 chk_intrs
 (
   double const * a,		// Pointer to previous ray position
   double const * b_i,		// Pointer to current ray position
   double const * u,		// Pointer to ray directional cosines
   double const l_i, 	            // Ray pathlength.
   int const NCL,		// Number of cylinders
   str_cyln const ** cylns,  // Pointer to cylinders
   int const NCN,		// Number of cones
   str_cone const ** cones,	// Pointer to cones
   int const NCB,
   str_cubd const ** cubds
 )
 {
   int i;				// Counting variable
   int res, res2;			// Intersection result
   int bs;				// Start index for loop over boxes
//   double invcdir[3];			// The inverse of the directional cosines
//   double invcdir_s[3];			// The sign of the inverse directional cosines
   double b[3] = {0.0};

   double l = l_i;
   if ( b_i[2] < 0.0 )
   {
     double l_vrt = b_i[2] / u[2];
     b[0] = b_i[0] + (-u[0] * l_vrt);
     b[1] = b_i[1] + (-u[1] * l_vrt);
     b[2] = 0.0; // to avoid rounding errors
     l -= l_vrt;
   } else {
     b[0] = b_i[0];
     b[1] = b_i[1];
     b[2] = b_i[2];
   }

   // Loop over boxes:
//   for (i = 0; i < 3; i++)
//   {
//     invcdir[i]   = 1.0 / u[i];
//     invcdir_s[i] = SIGN(invcdir[i]);
//   }

   // If there are at least three boxes, means the first is the bounding box of
   // all boxes. So it is sufficient to test that to know if it is necessary to
   // proceed to check the individual boxes. 
   bs = 0;
   if(NCB > 2)
   {
     res = intrs_cubd(a, b, u, cubds[0]);
     bs = ( res ) ? 1 : NCB;
   }

   for (i = bs; i < NCB; i++)
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




/* intrs_plane

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

 */

 void 
 intrs_plane
 (
   int    *intrs_f,
   double *intrs_l,
   double *intrs_p,
   const double *v,
   const double *u,
   const double l
 )
 {
   if ( ABS(v[2]) < TOLERANCE )
   {
    /* Ray starts on the plane */
     if ( ABS(u[2]) < TOLERANCE )
     {
      /* Ray's direction is parallel to the plane, intersection at the start */
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

/* intrs_rect

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

 */

 int
 intrs_rect
 ( 
   double const *v,
   double const *u,
   double const l,
   str_rect const *rect
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

/* intrs_ellp, intrs_ellp_circle

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

 */

 int
 intrs_ellp
 (
   double const *v,
   double const *u,
   double const l,
   str_ellp const *ellp
 )
 {
   int    intrs_f;
   double intrs_l;
   double intrs_p[3];
   intrs_plane(&intrs_f, &intrs_l, intrs_p, v, u, l);

   if ( intrs_f )
   {
     double rn_sq = (intrs_p[0] * intrs_p[0]) * ellp->rsq_inv[0] + 
                    (intrs_p[1] * intrs_p[1]) * ellp->rsq_inv[1];
     if ( (rn_sq - 1.0) < TOLERANCE ) 
     {
       return 1;
     }
   }
   return 0;
 }

/* intrs_cubd, intrs_cubd_closed

 Checks for intersection of a ray against a cuboid. It is assumed that the ray
 was first centered to the origin of the cuboid and rotated to its reference 
 frame. The function returns the intersection point if an intersection occurs.

 For a point P to be inside an axis-aligned rectangle described in its 
 supporting plane, it must satisfy:
 |P[X]| < (0.5 * L[X]) && |P[Y]| < (0.5 * L[Y]) && |P[Z]| < (0.5 * L[Z]),    (1)
 
 where L is the length of each dimension of the cuboid.

 The intersection point of a vector with the supporting plane of the rectangle 
 is described in the function intrs_plane. With the position of intersection on 
 the plane, if an intersection with the plane occurs, the intersection with the 
 rectangle if the point satisfy Eq. 1. 

 INPUTS:
 a    - Pointer to the initial point of the ray, centered and rotated;
 b    - Pointer to the terminal point of the ray, centered and rotated;
 u    - Pointer to the Cartesian direction of the ray;
 cubd - Pointer to an cuboid;
 
 OUTPUTS:
 A flag indicating if a intersection occurs.

 */

 int
 intrs_cubd
 (
   double const *a,
   double const *b,
   double const *u,
   str_cubd const *cubd
 )
 {
   double a_cr[3];	// Array with ray's initial point, centered and rotated
   double b_cr[3];	// Array with ray's terminal point, centered and rotated
   double u_rot[3];	// Array with ray's Cartesian directions, rotated

   if ( cubd->rotate_f[0] )
   {
     rot_vec(a_cr, a, cubd->o, (double const **) cubd->M);
     rot_vec(b_cr, b, cubd->o, (double const **) cubd->M);
     rot_vec_unit(u_rot, u, (double const **) cubd->M);
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
   double const *a,
   double const *b,
   double const *u,
   str_cubd const *cubd
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

/* intrs_cyln, intrs_cyln_right, intrs_cyln_oblique

 Checks for intersection of a ray against a circular cylinder (right or 
 oblique). Note that ray must be centered on the cylinder origin and that if the
 cylinder orientation axis is not aligned with the Z axis, the ray must be 
 rotated to the cylinder reference frame. This is normally performed by the 
 wrapper function intrs_cyln which also calls the appropriate function depending
 on cylinder type. The derivation below considers the case of a right-circular
 cylinder. Changes to an oblique-circular cylinder are considered in the end.

 In the following, all quantities are vector quantities (x, y, z). "P" is 
 the vector coordinates of an arbitrary point, "A" is the starting position 
 of the ray, "U" is the vector of cosine directions of the ray, "l" is the
 pathlength of the ray. For the cylinder, "O" is the center of the base 
 of the cylinder, "r" is its radius and "N" is a unit vector giving its 
 orientation (normal) and "h" is a distance along the cylinder axis. 
 
 The line equation is defined as:
 P = U * l + A.                                                              (1)
 
 The cylinder's lateral surface equation is given by:
 P = O + r * {cos(phi), sin(phi), h/r},                                      (2)
 
 Any point on the infinite right circular cylinder lateral surface must 
 satisfy:
 ||P - O||^2 - ((P - O) · N)^2 = r^2,                                        (3)

 where O is the center of the base of the cylinder, N is the normal to the base 
 of the cylinder (describes the orientation axis of the cylinder) and · is the
 dot (scalar) product.

 Intersections with the lateral surface of the cylinder have a P satisfying
 Eq. 1 and Eq. 3. Equating those equations we have: 
 ||P - O||^2 - ((P - O) · N)^2 = r^2 ->
 ||U * l + A - O||^2 - ((U * l + A - O) · N)^2 = r^2 ->
 ||U * l + (A - O)||^2 = ((U * l + (A - O)) · N)^2 + r^2.                    (4)

 Expanding the left-hand side we have:
 ||U * l + (A - O)||^2 ->
 ((U * l + (A - O)) · (U * l + (A - O))) ->
 (l^2 * (U · U) + 2 * l * (U · (A - O)) + ((A - O) · (A - O))^2) ->
 l^2 + 2 * l * (U · (A - O)) + ((A - O) · (A - O))^2,                        (5)
 
 where since U is a unit vector, U · U = 1. Now we can also expand the right 
 hand side of Eq. 4:
 ((U * l + (A - O)) · N)^2 + r^2 ->
 (l * (U · N) + ((A - O) · N))^2 + r^2 ->
 l^2 * (U · N)^2 + 2 * l * (U · N) * ((A - O) · N) + ((A - O) · N)^2 + r^2   (6)

 We can rearrange the expansions for the right and left-hand side of Eq. 4 by 
 combining Eq. 5 and Eq. 6:
 l^2 + 2 * l * (U · (A - O)) + ((A - O) · (A - O))^2 = 
   l^2 * (U · N)^2 + 2 * l * (U · N) * ((A - O) · N) + ((A - O) · N)^2 + 
   r^2 ->
 l^2 + 2 * l * (U · (A - O)) + ((A - O) · (A - O))^2 - 
   l^2 * (U · N)^2 - 2 * l * (U · N) * ((A - O) · N) - ((A - O) · N)^2 - 
   r^2 = 0 ->
 l^2 * (1 - (U · N)^2) + 
   2 * l * ((U · (A - O)) - (U · N) * ((A - O) · N)) + 
   ((A - O) · (A - O))^2 - ((A - O) · N)^2 - r^2 = 0,                        (7)

 which is a quadratic equation of the form a * x^2 + b * x + c = 0, with:
 a = 1 - (U · N)^2,
 b = 2 * ((U · (A - O)) - (U · N) * ((A - O) · N)),
 c = ((A - O) · (A - O))^2 - ((A - O) · N)^2 - r^2,
 delta = b^2 - 4 * a * c,
 l = - b +- sqrt(delta) / (2 * a).                                           (8)

 An intersection of the vector with the side of the cylinder occurs if l is 
 positive and lower than the magnitude of the vector. If a = 0, the b = 0 and
 the vector direction is parallel to the cylinder's axis (normal of the base).
 In this case an intersection with the infinite cylinder occurs only if the
 vector is on the surface of the cylinder, i.e., if c = 0. The intersection
 point(s) (or range) is given by substituting l from Eq. 8 in Eq. 1. If 
 (P - O) · N < than h, the height of the finite cylinder, the vector intersects
 the finite cylinder.
 
 If the vector is first translated by the cylinder origin and rotated to the 
 cylinder's reference frame, then Eq. 7 simplifies to:
 l^2 * (1 - U[Z]^2) + 
   2 * l * ((U · A) - U[Z] * A[Z]) + 
   (A · A)^2 - A[Z]^2 - r^2 = 0,                                             (9)

 which again is a quadratic equation of the form a * x^2 + b * x + c = 0, with:
 a = 1 - U[Z]^2,
 b = 2 * ((U · A) - U[Z] * A[Z]),
 c = (A · A)^2 - A[Z]^2 - r^2.

 An alternative presentation of Eq. 9 appearing in Leathers et al. (2001), 
 converts to the X and Y components:
 l^2 * (U[X]^2 + U[Y]^2) + 
   2 * l * (U[X] * A[X] + U[Y] * A[Y]) + 
   A[X]^2 + A[Y]^2 - r^2 = 0                                                (10)
 
 It is also possible that the the ray intercepts the base or top circles if the
 cylinder is closed or semi-closed. For a right circular cylinder, the
 intersection with the caps is simple to test since the caps have the same axis
 as the cylinder and no further rotation is necessary.

 For an oblique-circular or oblique in only one side, the function has to be 
 more general and therefore a bite less efficient than the simple test for the 
 right-circular cylinder, and separate functions were written.

 In an arbitrary oblique cylinder, the base and top have arbitrary orientation
 axes that may be different from the cylinder axis. In this version of the code,
 the X-axis of the opening is constrained to have the same direction as the
 X-axis of the cylinder, such that the only change is in relative polar angle to 
 the cylinder's axis. This is exemplified in the figure below for the base cap:

                      Alpha
                      ----->

                      Cylinder axis
                      |
 Loss in height |   \ |
                | ___\|___        _____ Z = 0.0 
                      \      |
                       \     | Gain in height
                        Cap
                  ____
                  Radius


 The gain in height is maximum at and azimuth of 0.0º and minimum (loss) at an
 azimuth of 180.0º. The opposite is true for the top cap. The change is 0 at 
 90.0º. The height change is +- r_x * sin(theta) so:

 DeltaH = r_x * -sin(theta) * P[x] / r_x ->
 DeltaH = -sin(theta) * P[x] (1)

 where r_x is the major axis (here equal to radius in X) of the ellipse formed 
 by the opening, theta is the angle of the cap's normal relative to the cylinder
 axis. Therefore, before the check to see if intersection occurs within the 
 height of the finite cylinder, the height is adjusted by DeltaH. Further, the
 code tests if it is necessary to rotate vectors from axis to base and from
 base to top.
 
 INPUTS:
 a    - Pointer to the initial point of the ray;
 b    - Pointer to the terminal point of the ray; 
 u    - Pointer to the direction of the ray; 
 l    - Magnitude of the ray; 
 cyln - Pointer to cylinder structure.

 OUTPUT:
 Flag indicating if an intersection occurred.

 */

 int
 intrs_cyln
 (
   double const * a, 
   double const * b, 
   double const * u, 
   double const l, 
   str_cyln const * cyln
 )
 {
   double a_cr[3];	// Array with ray's initial point, centered and rotated
   double b_cr[3];	// Array with ray's terminal point, centered and rotated 
   double u_rot[3];	// Array with ray's Cartesian directions, rotated
   double c_vec[3];

   if ( cyln->rotate_f[0] )
   {
     ROT_VEC(a_cr, a, cyln->o, cyln->M, c_vec);
     ROT_VEC(b_cr, b, cyln->o, cyln->M, c_vec);
     ROT_VEC_UNIT(u_rot, u, cyln->M);
   }
   else
   {
     for (size_t i = 0; i < 3; i++)
     {
       a_cr[i]  = a[i] - cyln->o[i];
       b_cr[i]  = b[i] - cyln->o[i];
       u_rot[i] = u[i];
     }
   }

   if ( cyln->rotate_f[1] || cyln->rotate_f[2] )
   {

     return intrs_cyln_oblique(a_cr, b_cr, u_rot, l, cyln);
   } else {

     return intrs_cyln_right(a_cr, b_cr, u_rot, l, cyln);
   }

 }

 int
 intrs_cyln_right
 (
   double *a, 
   double const *b, 
   double const *u, 
   double const l, 
   str_cyln const *cyln
 )
 {
   double intrs_l;		// Pathlength to intersection
   double intrs_p[3];		// Position of intersection
   double A, B, C, delta; 	// Terms of the Bhaskara equation
   double sdelta;		// Squared root of delta

   A = 1.0 - u[2] * u[2];
   C = a[0] * a[0] + a[1] * a[1] - cyln->rsq;

   if ( ABS(A) < TOLERANCE )
   {
     // If A == 0, B will necessarily be 0 and ray is traveling parallel to the 
     // cylinder axis. If C = 0, vector is on the surface of the infinite 
     // cylinder.
     if ( ABS(C) < TOLERANCE )
     {
       // If the ray is on the surface of the infinite cylinder, an intersection
       // with the finite cylinder occurs if points of the ray and cylinder
       // overlap in a range.
       if ( u[2] >= 0.0 )
       {
         return (b[2] >= 0.0) && (cyln->h >= a[2]);
       } else {
         return (a[2] >= 0.0) && (cyln->h >= b[2]);
       }
     }
   }
   else
   {
     B     = (2 * u[0] * a[0]) + (2 * u[1] * a[1]);
     delta = (B * B) - (4 * A * C);

     if ( delta >= 0.0 )
     {
       // If delta >= 0, the equality has real roots and the line intersects the 
       // infinite cylinder surface. To know if the ray intersects the cylinder, 
       // we check for each root, if it is positive, if it smaller than the 
       // pathlength of the ray and if the intersection vertical position is 
       // within the vertical extent of the finite cylinder.
       sdelta = sqrt(delta);

       intrs_l = ((-B - sdelta) / (2.0 * A));
       if ( (intrs_l > 0.0) && (intrs_l <= l) )
       {
         intrs_p[2] = a[2] + u[2] * intrs_l;

         if ( (intrs_p[2] > 0.0) && (intrs_p[2] <= cyln->h) )
         {
           return 1;
         }
       }

       intrs_l = ((-B + sdelta) / (2.0 * A));
       if ( (intrs_l > 0.0) && (intrs_l <= l) )
       {
         intrs_p[2] = a[2] + u[2] * intrs_l;

         if ( (intrs_p[2] > 0.0) && (intrs_p[2] <= cyln->h) )
         {
           return 1;
         }
       }
     }
   }
   // Intersection with horizontal caps:
   // If no intersections occurred at the side surface, it might still occur 
   // with the caps, if cylinder is semi-closed or closed. For the base opening,
   // no manipulation is necessary. For the top opening, the vertical position,
   // centered on the base, needs first to be re-centered on the top. Since this 
   // is a right circular cylinder, X and Y positions do not change.
   if ( cyln->closed_f[0] )
   {
     if ( intrs_ellp(a, u, l, cyln->base) )
     {
       return 1;
     }
   }
   if ( cyln->closed_f[1] )
   {
     for (size_t i = 0; i < 3; i++)
       a[i] -= cyln->h;
     return intrs_ellp(a, u, l, cyln->top);
   }
   return 0;
 }

 int
 intrs_cyln_oblique
 (
   double const *a, 
   double const *b, 
   double const *u, 
   double const l, 
   str_cyln const *cyln
 )
 {
   double a_rot[3];		// Vector of position rotated for the caps
   double u_rot[3];		// Vector of directions rotated for the caps
   double dh[2];		// Delta H, change in height due to cap inclination
   double intrs_l;		// Pathlength to intersection
   double intrs_p[3];		// Position of intersection vector
   double intrs_dhf;		// Fractional Delta H based on azimuth for each opening
   double A, B, C, delta; 	// Terms of the Baskara equation
   double sdelta;		// Squared root of delta

   A = 1.0 - u[2] * u[2];
   C = a[0] * a[0] + a[1] * a[1] - cyln->rsq;

   if ( ABS(A) < TOLERANCE )
   {
     // If A == 0, B will necessarily be 0 and ray is traveling parallel to the 
     // cylinder axis. If C = 0, vector is on the surface of the infinite 
     // cylinder.
     if ( ABS(C) < TOLERANCE )
     {
       // If the ray is on the surface of the infinite cylinder, an intersection
       // with the finite cylinder occurs if points of the ray and cylinder 
       // overlap in a range.
       intrs_dhf = a[0];
       dh[0] = cyln->dhmx[0] * intrs_dhf;
       dh[1] = cyln->dhmx[1] * intrs_dhf;
       if ( u[2] >= 0.0 )
       {
         return (b[2] >= dh[0]) && ((cyln->h + dh[1]) >= a[2]);
       } else {
         return (a[2] >= dh[0]) && ((cyln->h + dh[1]) >= b[2]);
       }
     }
   }
   else
   {
     B     = (2 * u[0] * a[0]) + (2 * u[1] * a[1]);
     delta = (B * B) - (4 * A * C);
     if ( delta >= 0.0 )
     {
       // If delta >= 0, the equality has real roots and the line intersects the 
       // infinite cylinder surface. To know if the ray intersects the cylinder, 
       // we check for each root, if it is positive, if it smaller than the 
       // pathlength of the ray and if the intersection vertical position is 
       // within the vertical extent of the finite cylinder.

       sdelta = sqrt(delta);

       intrs_l = ((-B - sdelta) / (2.0 * A));
       if ( (intrs_l > 0.0) && (intrs_l <= l) )
       {
         intrs_p[0] = a[0] + u[0] * intrs_l;
         intrs_p[2] = a[2] + u[2] * intrs_l;
         intrs_dhf = intrs_p[0];
         dh[0] = cyln->dhmx[0] * intrs_dhf;
         dh[1] = cyln->dhmx[1] * intrs_dhf;
         if ( (intrs_p[2] > dh[0]) && (intrs_p[2] <= (cyln->h + dh[1])) )
         {
           return 1;
         }
       }

       intrs_l = ((-B + sdelta) / (2.0 * A));
       if ( (intrs_l > 0.0) && (intrs_l <= l) )
       {
         intrs_p[0] = a[0] + u[0] * intrs_l;
         intrs_p[2] = a[2] + u[2] * intrs_l;
         intrs_dhf = intrs_p[0];
         dh[0] = cyln->dhmx[0] * intrs_dhf;
         dh[1] = cyln->dhmx[1] * intrs_dhf;
         if ( (intrs_p[2] > dh[0]) && (intrs_p[2] <= (cyln->h + dh[1])) )
         {
           return 1;
         }
       }
     }
   }
   // Intersection with horizontal caps:
   // If no intersections occurred at the side surface, it might still occur 
   // with the caps, if cylinder is semi-closed or closed. For the base opening,
   // no manipulation is necessary. For the top opening, the vertical position,
   // centered on the base, needs first to be re-centered on the top. 
   // Additionally, because at least one opening is oblique, it might be
   // necessary to rotate the ray from the cylinder axis to the base or top
   // axis.
   double c_vec[3];
   if ( cyln->closed_f[0] )
   {
     for (size_t i = 0; i < 3; i++)
     {
       a_rot[i] = a[i];
       u_rot[i] = u[i];
     }
     if ( cyln->rotate_f[1] )
     {
       ROT_VEC_UNIT(a_rot, a, cyln->base->M);
       ROT_VEC_UNIT(u_rot, u, cyln->base->M);
     }
     if ( intrs_ellp(a_rot, u_rot, l, cyln->base) )
     {
       return 1;
     }
   }

   if ( cyln->closed_f[1] )
   {
     for (size_t i = 0; i < 3; i++)
     {
       a_rot[i] = a[i];
       u_rot[i] = u[i];
     }
     if ( cyln->rotate_f[3] )
     {
       if ( cyln->rotate_f[2] )
       {
         ROT_VEC(a_rot, a, cyln->o, cyln->top->M, c_vec);
         ROT_VEC_UNIT(u_rot, u, cyln->top->M);
       }
       else
       {
         a_rot[2] -= cyln->h;
       }
     }
     return intrs_ellp(a_rot, u_rot, l, cyln->top);
   }
   return 0;
 }



/* Line segment and cone intersection:
 *
 * In the folling, all quantities are vector quantities (x, y, z). "pos" is the 
 * vector cordinates of an arbitrary point, "ppos" is the starting position of 
 * the ray, "cdir" is the vector of cosine directions of the ray, "s" is the
 * path length of the ray. For the cone, "vertex" is the vertex of the cone, 
 * "theta" is the halfangle of the cone and "axis" is a unit vector giving the 
 * orientation of the cone. "phi" is the azimuthal angle.
 *
 * The line equation is defined as:
 * pos = ppos + cdir * s                                                     (1)
 *
 * The cone's surface equation is given by:
 * pos = vertex +- (dir · axis) * length,                                    (2)
 *
 * with "dir" is given by the polar coordinates {theta, phi} and phi being in 
 * the range [0, 2*pi). We will restrict our self to a single sided cone, taken 
 * here as the positive part with regard to its vertex. Any point on this cone's 
 * surface must satisfy:
 * (pos - vertex) · axis = ||pos - vertex|| * cos(theta)                     (3)
 *
 * Intersections with the surface of the cone are "pos" equal for the line and 
 * for the cone's surface. By squaring Eq. 3 and substituting pos for the line 
 * equation we have:
 * ((ppos + cdir * s - vertex) · axis)^2 = 
 *   ||ppos + cdir * s - vertex||^2 * cos(theta)^2 ->
 *
 * ((s * cdir + ppos - vertex) · axis)^2 = 
 *   ||s * cdir + ppos - vertex||^2 * cos(theta)^2 ->
 * 
 * ((s * cdir + ppos_v) · axis)^2 = ||s * cdir + ppos_v||^2 * mu^2,          (4)
 * 
 * Where we are using ppos_v to represent ppos centered by the vertex and mu
 * to represent cos(theta). Expanding the left hand side we have:
 * ((s * cdir + ppos_v) · axis)^2 ->
 *
 * (s * cdir · axis + ppos_v · axis)^2 ->
 *
 * s^2 * (cdir · axis)^2 + 2 * s * (cdir · axis) * (ppos_v · axis) + 
 *   (ppos_v · axis)^2 ->
 * 
 * s^2 * cdir[2]^2 + 2 * s * cdir[2] * ppos_v[2] + ppos_v[2]^2.              (5)
 *
 * In the last step to reach Eq. 5 we have assumed the less generic case in 
 * which the axis of the cone is allinged with the z axis of the reference frame 
 * (or equivalently, that cdir has been projected into the arbitrary axis 
 * direction of the cone). Additionally, "[2]" represents the z axis position 
 * (since indexes in C strat at 0, 2 is the third index. Now we can also expand 
 * the right hand side of Eq. 4:
 * ||s * cdir + ppos_v||^2 * mu^2 ->
 *
 * ((s * cdir + ppos_v) · (s * cdir + ppos_v)) * mu^2 ->
 *
 * (s^2 * (cdir · cdir) + 2 * s * (cdir · ppos_v) + (ppos_v · ppos_v)^2) * 
 *   mu^2 ->
 * 
 * (s^2 + 2 * s * (cdir · ppos_v) + (ppos_v · ppos_v)^2) * mu^2 -> 
 * 
 * s^2 * mu^2 + 2 * s * (cdir · ppos_v) * mu^2 + (ppos_v · ppos_v)^2 * mu^2, (6)
 *
 * where since cdir is a unit vector, cdir · cdir = 1. We can rearrange the 
 * expansions for the righ and left side of Eq. 4 by combining equation 5 and 6:
 * s^2 * cdir[2]^2 + 2 * s * cdir[2] * ppos_v[2] + ppos_v[2]^2 - 
 *   (s^2 * mu^2 + 2 * s * (cdir · ppos_v) * mu^2 + 
 *   (ppos_v · ppos_v)^2 * mu^2) = 0 ->
 * 
 * s^2 * cdir[2]^2 - s^2 * mu^2 + 2 * s * cdir[2] * ppos_v[2] - 
 *   2 * s * (cdir · ppos_v) * mu^2 + ppos_v[2]^2 - 
 *   (ppos_v · ppos_v)^2 * mu^2 ->
 *
 * s^2 * (cdir[2]^2 - mu^2) + 2 * s * (cdir[2] * ppos_v[2] - (cdir · ppos_v) * 
 *   mu^2) + ppos_v[2]^2 - (ppos_v · ppos_v)^2 * mu^2,                       (7)
 * 
 * which is a quadratic equation of the form A * x^2 + B * x + C = 0, with:
 * A = cdir[2]^2 - mu^2,
 * B = 2 * cdir[2] * ppos_v[2] - (cdir · ppos_v) * mu^2,
 * C = ppos_v[2]^2 - (ppos_v · ppos_v)^2 * mu^2
 * delta = B^2 - 4 * A * C
 * s = - B +- sqrt(delta) / (2 * A)
 *
 * if A == 0, the inclination of the ray is the same as of the cone's surface 
 * and the solution is linear with:
 * s = - C / (2 * B)
 *
 * and if A == 0 and B == 0, the line segment passes through the vertex.
 *
 * if A != 0, and if delta < 0, no intersection occours. If delta = 0, the line 
 * segment intersect the cone at one point and if delta > 0, the line segment 
 * intersect the cone in two points (or many points if it lies on the surface of 
 * the cone). For a finite and optionally truncated cone, it is necessary to 
 * test it the point of intersecton is in the range of the cone height relative 
 * to vertex.
 * 
 * If it is a hollow cone, those tests are sufficient. If it is a closed or
 * semi-closed cone, it is also necessary to check for intersection at the 
 * circular caps.
 *
 * It is important to remember that this function assumes that the cone's axis 
 * is alligned with the Z axis, therefore is necessary to project the ray cosine 
 * direction into the reference frame of the cone.
 * 
 */
/*
 int intrs_cone (const double *ppos, 
                 const double *cpos, 
                 const double *cdir, 
                 const double s, 
                 const struct str_cone *cone,
                 const int    segment)
 {
   double ppos_v[3];
   double cdotpv, pvdotpv, cdota, pvdota;
   double delta, sdelta;
   double s_intrs;
   double p_intrs[3];
   double A, B, C; 
   double p_cdir;

   ppos_v[0] = ppos[0] - (*cone).vertex[0];
   ppos_v[1] = ppos[1] - (*cone).vertex[1];
   ppos_v[2] = ppos[2] - (*cone).vertex[2];

   cdotpv  = cdir[0] * ppos_v[0] + 
             cdir[1] * ppos_v[1] + 
             cdir[2] * ppos_v[2];

   pvdotpv = ppos_v[0] * ppos_v[0] + 
             ppos_v[1] * ppos_v[1] + 
             ppos_v[2] * ppos_v[2];

/*   If ray cosine directions are not projected, then the equation would be:

   cdota   = cdir[0] * (*cone).caxis[0] + 
             cdir[1] * (*cone).caxis[1] + 
             cdir[2] * (*cone).caxis[2];

   pvdota  = ppos_v[0] * (*cone).caxis[0] + 
             ppos_v[1] * (*cone).caxis[1] + 
             ppos_v[2] * (*cone).caxis[2];

   A = cdota * cdota - (*cone).mu_sq;
   B = 2.0 * (cdota * pvdota - cdotpv * (*cone).mu_sq);
   C = pvdota * pvdota - pvdotpv * (*cone).mu_sq;
*/
/*
   A = cdir[2] * cdir[2] - (*cone).mu_sq;
   B = 2.0 * (cdir[2] * ppos_v[2] - cdotpv * (*cone).mu_sq);
   C = ppos_v[2] * ppos_v[2] - pvdotpv * (*cone).mu_sq;

   if ( NUM_EQU(A, 0.0, 1E-12) )
   {
     // If A == 0, quadratic term is zero and root is linear: B * x + C = 0
     // -> x = -C / B.
     if ( NUM_EQU(B, 0.0, 1E-12) )
     {
       // If A and B are zero, line would pass through the vertex. Line segment
       // migt not, but still intersection occours as path is through the 
       // surface of the cone. All is needed is that the z ranges intersect.
       if ( NUM_EQU((*cone).height[0], 0.0, 1E-12) )
       { 
         if (cdir[2] > 0)
         {
           if ( ( (ppos[2] > ((*cone).vertex[2] + (*cone).height[0])) && 
                  (ppos[2] < ((*cone).vertex[2] + (*cone).height[1])) ) || 
                ( (ppos[2] < ((*cone).vertex[2] + (*cone).height[0])) && 
                  (cpos[2] > ((*cone).vertex[2] + (*cone).height[0])) ) )
           {
             return 1;
           }
         } else {
           if ( ( (ppos[2] < ((*cone).vertex[2] + (*cone).height[1])) && 
                  (ppos[2] > ((*cone).vertex[2] + (*cone).height[0])) ) || 
                ( (ppos[2] > ((*cone).vertex[2] + (*cone).height[1])) && 
                  (cpos[2] < ((*cone).vertex[2] + (*cone).height[1])) ) )
           {
             return 1;
           }
         }
       }
     } else {
       // If B is not zero, the root (= -C / B) gives the single point of 
       // intersection. If this root is positive and lower than the panthlength 
       // s of the ray, intersection is possible. However, we whant only the 
       // positive part of the cone, so we can project the vertex centered 
       // vector of the intersection point into the cone axis and this 
       // projection has to be positive.  
       s_intrs = -C / B;
       if ( (s_intrs > 0.0) && (s_intrs <= s) )
       {
         p_intrs[0] = ppos_v[0] + cdir[0] * s_intrs; 
         p_intrs[1] = ppos_v[1] + cdir[1] * s_intrs;
         p_intrs[2] = ppos_v[2] + cdir[2] * s_intrs;
         p_cdir = p_intrs[0] * (*cone).caxis[0] + 
                  p_intrs[1] * (*cone).caxis[1] +
                  p_intrs[2] * (*cone).caxis[2];
         if ( p_cdir > 0.0 && (p_intrs[2] > (*cone).height[0]) && 
           (p_intrs[2] < (*cone).height[1]))
         {
           return 1;
         }
       }
     }    
   } else {
     delta = B * B - 4 * A * C;
     if (delta >= 0.0) 
     {
       sdelta = sqrt(delta);

       s_intrs = ((-B - sdelta) / (2 * A));
       if ( (s_intrs > 0.0) && (s_intrs <= s) )
       {
         p_intrs[0] = ppos_v[0] + cdir[0] * s_intrs; 
         p_intrs[1] = ppos_v[1] + cdir[1] * s_intrs;
         p_intrs[2] = ppos_v[2] + cdir[2] * s_intrs;
         p_cdir = p_intrs[0] * (*cone).caxis[0] + 
                  p_intrs[1] * (*cone).caxis[1] +
                  p_intrs[2] * (*cone).caxis[2];
         if ( (p_cdir > 0.0) && (p_intrs[2] > (*cone).height[0]) && 
           (p_intrs[2] < (*cone).height[1]) )
         {
           return 1;
         }
       }

       s_intrs = ((-B + sdelta) / (2 * A));
       if ( (s_intrs > 0.0) && (s_intrs <= s) )
       {
         p_intrs[0] = ppos_v[0] + cdir[0] * s_intrs; 
         p_intrs[1] = ppos_v[1] + cdir[1] * s_intrs;
         p_intrs[2] = ppos_v[2] + cdir[2] * s_intrs;
         p_cdir = p_intrs[0] * (*cone).caxis[0] + 
                  p_intrs[1] * (*cone).caxis[1] +
                  p_intrs[2] * (*cone).caxis[2];
         if ( (p_cdir > 0.0) && (p_intrs[2] > (*cone).height[0]) && 
           (p_intrs[2] < (*cone).height[1]) )
         {
           return 1;
         }
       }
     }
   }

   // If no intersections happened at the side surface, it might still be 
   // necessary to test against circular caps if semi-closed or closed:
   if ( (*cone).closed[0] )
   {
     s_intrs = ((*cone).height[0] - ppos_v[2]) / cdir[2];
     if (s_intrs > 0 && s_intrs < s)
     {
       p_intrs[0] = ppos_v[0] + cdir[0] * s_intrs; 
       p_intrs[1] = ppos_v[1] + cdir[1] * s_intrs;
       if ( sqrt(p_intrs[0] * p_intrs[0] + p_intrs[1] * p_intrs[1]) <= 
            (*cone).radius[0] )
       {
         return 1;
       }
     }
   }

   if ( (*cone).closed[1] )
   {
     s_intrs = ((*cone).height[1] - ppos_v[2]) / cdir[2];
     if (s_intrs > 0 && s_intrs < s)
     {
       p_intrs[0] = ppos_v[0] + cdir[0] * s_intrs; 
       p_intrs[1] = ppos_v[1] + cdir[1] * s_intrs;
       if ( sqrt(p_intrs[0] * p_intrs[0] + p_intrs[1] * p_intrs[1]) <= 
            (*cone).radius[1] )
       {
         return 1;
       }
     }
   }

   return 0;
 }

 void setup_cone (struct str_cone *cone,
                  const double *cn_vertex,
                  const double *cn_axis,
                  const double *cn_height,
                  const double cn_theta,
                  const int    *cn_closed)
 {
   (*cone).vertex[0] = cn_vertex[0];
   (*cone).vertex[1] = cn_vertex[1];
   (*cone).vertex[2] = cn_vertex[2];

   (*cone).paxis[0] = cn_axis[0];
   (*cone).paxis[1] = cn_axis[1];

   (*cone).caxis[0] = sin(cn_axis[0]) * cos(cn_axis[1]);
   (*cone).caxis[1] = sin(cn_axis[0]) * sin(cn_axis[1]);
   (*cone).caxis[2] = cos(cn_axis[0]);

   (*cone).psi   = cn_theta;
   (*cone).mu    = cos(cn_theta);
   (*cone).mu_sq = (*cone).mu * (*cone).mu;

   (*cone).height[0] = cn_height[0];
   (*cone).height[1] = cn_height[1];

   (*cone).radius[0] = (*cone).height[0] * tan(cn_theta);
   (*cone).radius[1] = (*cone).height[1] * tan(cn_theta);

   (*cone).closed[0] = cn_closed[0];
   (*cone).closed[1] = cn_closed[1];

   (*cone).rotate_f = ((*cone).caxis[2] < 1.0)? 1 : 0;
 }

*/


