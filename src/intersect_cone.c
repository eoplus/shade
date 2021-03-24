
/*******************************************************************************
 intersect_cone.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 Version: 2021-03-12

 intrs_cone, intrs_cone_right, intrs_cone_oblique

 This source file provides functions to check for intersection of a ray against
 a circular cone (right or oblique).

 The specific functions assume that the ray's initial and terminal points are 
 centered on the cone's origin (vertex) and that if the cone orientation axis
 does not have the same direction as the system's Z-axis, the ray points and 
 direction were rotated to the cone's reference frame. Those transformations
 are performed by the wrapper function intrs_cone which also calls the
 appropriate function depending on cone type. The derivation below considers
 the case of a right-circular positive cone. Changes to an oblique-circular 
 positive cone are considered in the end.

 In the following, all quantities in upper case are vector quantities in three 
 dimensions (x, y, z). "P" is the vector coordinates of an arbitrary point, "A"
 is the starting position of the ray, "U" is the vector of Cartesian directions
 of the ray, and "l" is the (scalar) pathlength of the ray. For the cone,
 "O" is the vertex (origin), G is the generatrix direction (given by psi, the
 cone half-angle), "N" is a unit vector giving the cone orientation (axis) and
 "h" is a height as distance along the cone axis. 

 The line equation is defined as:
 P = U * l + A.                                                              (1)
 
 The positive cone's surface equation is given by:
 P = O + (G · N) * h,	                                                     (2)

 where O is the vertex of the cone, N is the orientation axis of the cone, G is
 the generatrix direction relative to the axis, · is the dot (scalar) product
 and (G · N) * h = cos(psi) * h gives the slant height. Any point on this cone's
 surface must satisfy:
 (P - O) · N = ||P - O|| * cos(psi),			                     (3)

 where || || is the Euclidean norm and psi is the half-angle of the cone
 (= acos(G · N)). Intersections with the lateral surface of the cone have a P
 satisfying Eq. 1 and Eq. 3. Substituting Eq. 1 into Eq. 3 and squaring it, we
 have:
 ((U * l + (A - O)) · N)^2 = ||U * l + (A - O)||^2 * mu^2,                   (4)
  
 where we are using "mu" to represent cos(psi). Expanding the left-hand side
 we have:
 ((U * l + (A - O)) · N)^2 ->
 (l * U · N + (A - O) · N)^2 ->
 l^2 * (U · N)^2 + 2 * l * (U · N) * ((A - O) · N) + ((A - O) · N)^2.        (5)

 Now we can also expand the right-hand side of Eq. 4:
 ||U * l + (A - O)||^2 * mu^2 ->
 ((U * l + (A - O)) · (U * l + (A - O))) * mu^2 ->
 (l^2 * (U · U) + 2 * l * (U · (A - O)) + ((A - O) · (A - O))^2) * mu^2 ->
 (l^2 + 2 * l * (U · (A - O)) + ((A - O) · (A - O))^2) * mu^2 -> 
 l^2 * mu^2 + 2 * l * (U · (A - O)) * mu^2 + ((A - O) · (A - O))^2 * mu^2,   (6)
 
 where since U is a unit vector, U · U = 1.0 and the term is removed. We can
 rearrange the expansions for the right and left-hand side of Eq. 4 by combining
 Eq. 5 and 6:
 l^2 * (U · N)^2 + 2 * l * (U · N) * ((A - O) · N) + ((A - O) · N)^2 = 
   l^2 * mu^2 + 2 * l * (U · (A - O)) * mu^2 + ((A - O) · (A - O))^2 * mu^2 ->
 l^2 * (U · N)^2 + 2 * l * (U · N) * ((A - O) · N) + ((A - O) · N)^2 - 
   (l^2 * mu^2 + 2 * l * (U · (A - O)) * mu^2 + 
   ((A - O) · (A - O))^2 * mu^2) = 0 ->
 l^2 * (U · N)^2 - l^2 * mu^2 + 2 * l * (U · N) * ((A - O) · N) - 
   2 * l * (U · (A - O)) * mu^2 + ((A - O) · N)^2 - 
   ((A - O) · (A - O))^2 * mu^2 = 0 ->
 l^2 * ((U · N)^2 - mu^2) + 
   2 * l * ((U · N) * ((A - O) · N) - (U · (A - O)) * mu^2) +
   ((A - O) · N)^2 - ((A - O) · (A - O))^2 * mu^2 = 0,                       (7)

 which is a quadratic equation of the form a * x^2 + b * x + c = 0, with:
 a = (U · N)^2 - mu^2,
 b = 2 * ((U · N) * ((A - O) · N) - (U · (A - O)) * mu^2),
 c = ((A - O) · N)^2 - ((A - O) · (A - O))^2 * mu^2,
 delta = b^2 - 4 * a * c,
 l = - b +- sqrt(delta) / (2 * a).                                           (8)

 If A == 0, the inclination of the ray is the same as of the cone's surface and
 the solution is linear with:
 l = - c / b.                                                                (9)
 
 If A == 0 and B == 0, the ray passes through the vertex. If A != 0, and if
 delta < 0, no intersection occurs. If delta = 0, the ray intersects the cone
 at one point and if delta > 0, the ray intersects the cone in two points (or
 many points if it lies on the surface of the cone). For a finite and optionally
 truncated cone, it is necessary to test it the point of intersection is in the
 range of the cone height relative to vertex.

 If it is a hollow cone, those tests are sufficient. If it is a closed or 
 semi-closed cone, it is also necessary to check for intersection at the caps.

 If the axis of the cone is aligned with the axis of the system (conversely,
 if the ray was centered and rotated to the cone reference frame), Eq. 7
 simplifies to:
 l^2 * (U[Z]^2 - mu^2) + 
   2 * l * (U[Z] * (A - O)[Z] - (U · (A - O)) * mu^2) +
   (A - O)[Z]^2 - ((A - O) · (A - O))^2 * mu^2 = 0,                         (10)

 which again is a quadratic equation of the form a * x^2 + b * x + c = 0, with:
 a = U[Z]^2 - mu^2,
 b = 2 * (U[Z] * (A - O)[Z] - (U · (A - O)) * mu^2),
 c = (A - O)[Z]^2 - ((A - O) · (A - O))^2 * mu^2.

 For an oblique-circular or oblique in only one side, the function has to be 
 more general and therefore a bit less efficient than the simple test for the 
 right-circular cone, and separate functions were written.

 In an arbitrary oblique cone, the base and top have arbitrary orientation
 axes that may be different from the cone axis. In this version of the code,
 the X-axis of the opening is constrained to have the same direction as the
 X-axis of the cone, such that the only change is in the relative polar angle to 
 the cone's axis. This is exemplified in the figure below for the base cap:

                      Alpha
                      ----->

                      Cone axis
                      |
 Loss in height |   \ |
                | ___\|___        _____ Z = height[0] 
                      \      |
                       \     | Gain in height
                        Cap
                  ____
                  Radius

 The gain in height is maximum at an azimuth of 0.0º and minimum (loss) at an
 azimuth of 180.0º. The same is true for the top cap. The change is 0 at 90.0º.
 The height change is +- r_x * sin(theta) so:

 DeltaH = r_x * tan(theta) * P[x] / r_x ->
 DeltaH = tan(theta) * P[x]                                                (11)

 where r_x is the major axis (here equal to radius in X) of the ellipse formed 
 by the opening, theta is the angle of the cap's normal relative to the cone
 axis. Therefore, before the check to see if intersection occurs within the 
 height of the finite cone, the height is adjusted by DeltaH. Further, the
 code tests if it is necessary to rotate vectors from axis to base and from
 base to top.
 
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
 l    - Pathlength of the ray
        Range: (0,Inf), meters
        Constant double;
 cone - Cone parameters
        Pointer to constant str_cone struct.

 OUTPUT:
 Flag indicating if an intersection occurred: (0) No, (1) Yes.

*******************************************************************************/

 #include <stddef.h>		// size_t
 #include <math.h>		// sqrt

 #include "config.h"		// TOLERANCE
 #include "aux.h"		// ABS
 #include "rotation.h"		// ROT_VEC, ROT_VEC_UNIT
 #include "structures.h"	// str_cone
 #include "intersect.h"		// intrs_ellp

 int
 intrs_cone
 (
   double const * a, 
   double const * b, 
   double const * u, 
   double const l, 
   struct str_cone const * cone
 )
 {
   double a_cr[3];	// Ray's initial point, centered and rotated
   double b_cr[3];	// Ray's terminal point, centered and rotated 
   double u_rot[3];	// Ray's Cartesian directions, rotated
   double c_vec[3];	// Temporary variable for the macro ROT_VEC

   if ( cone->rotate_f[0] )
   {
     ROT_VEC(a_cr, a, cone->o, cone->M, c_vec);
     ROT_VEC(b_cr, b, cone->o, cone->M, c_vec);
     ROT_VEC_UNIT(u_rot, u, cone->M);
   }
   else
   {
     for (size_t i = 0; i < 3; i++)
     {
       a_cr[i]  = a[i] - cone->o[i];
       b_cr[i]  = b[i] - cone->o[i];
       u_rot[i] = u[i];
     }
   }

   if ( cone->rotate_f[1] || cone->rotate_f[2] )
   {
     return intrs_cone_oblique(a_cr, b_cr, u_rot, l, cone);
   } else {
     return intrs_cone_right(a_cr, b_cr, u_rot, l, cone);
   }
 }


 int
 intrs_cone_right
 (
   double const * a, 
   double const * b, 
   double const * u, 
   double const l, 
   struct str_cone const * cone
 )
 {
   double intrs_l;
   double intrs_p[3];
   double intrs_mu;

   double const cdotpv  = u[0] * a[0] + 
                          u[1] * a[1] + 
                          u[2] * a[2];

   double const pvdotpv = a[0] * a[0] + 
                          a[1] * a[1] + 
                          a[2] * a[2];

   double const A = u[2] * u[2] - cone->musq;
   double const B = 2.0 * (u[2] * a[2] - cdotpv * cone->musq);
   double const C = a[2] * a[2] - pvdotpv * cone->musq;

   if ( ABS(A) < TOLERANCE )
   {
     // If A == 0, quadratic term is zero and root is linear: B * x + C = 0
     // -> x = -C / B.
     if ( ABS(B) < TOLERANCE  )
     {
       // If A and B are zero, line would pass through the vertex. Line segment
       // migt not, but still intersection occours as path is through the 
       // surface of the cone. All is needed is that the z ranges intersect.
       if ( cone->h[0]  < TOLERANCE )
       { 
         if ( u[2] > 0.0 )
         {
           if ( ( (a[2] > cone->h[0]) && (a[2] < cone->h[1]) ) || 
                ( (a[2] < cone->h[0]) && (b[2] > cone->h[0]) ) )
           {
             return 1;
           }
         } else {
           if ( ( (a[2] < cone->h[1]) && (a[2] > cone->h[0]) ) || 
                ( (a[2] > cone->h[1]) && (b[2] < cone->h[1]) ) )
           {
             return 1;
           }
         }
       }
     } else {
       // If B is not zero, the root (= -C / B) gives the single point of 
       // intersection. If this root is positive and lower than the panthlength 
       // l of the ray, intersection is possible. However, we want only the 
       // positive part of the cone, so we can project the vertex centered 
       // vector of the intersection point into the cone axis and this 
       // projection has to be positive.  
       intrs_l = -C / B;
       if ( (intrs_l > 0.0) && (intrs_l <= l) )
       {
         for (size_t i = 0; i < 3; i++)
           intrs_p[i] = a[i] + u[i] * intrs_l;
         intrs_mu = intrs_p[0] * cone->u[0] + 
                    intrs_p[1] * cone->u[1] +
                    intrs_p[2] * cone->u[2];
         if ( (intrs_mu > 0.0) && (intrs_p[2] > cone->h[0]) && 
              (intrs_p[2] < cone->h[1]) )
           return 1;
       }
     }    
   }
   else
   {
     double const delta = B * B - 4.0 * A * C;
     if ( delta >= 0.0 ) 
     {
       double sdelta = sqrt(delta);
       double const k_1_2A = 1.0 / (2.0 * A);

       intrs_l = (-B - sdelta) * k_1_2A;
       if ( (intrs_l > 0.0) && (intrs_l <= l) )
       {
         for (size_t i = 0; i < 3; i++)
           intrs_p[i] = a[i] + u[i] * intrs_l;
         intrs_mu = intrs_p[0] * cone->u[0] + 
                    intrs_p[1] * cone->u[1] +
                    intrs_p[2] * cone->u[2];
         if ( (intrs_mu > 0.0) && (intrs_p[2] > cone->h[0]) && 
              (intrs_p[2] < cone->h[1]) )
           return 1;
       }

       intrs_l = (-B + sdelta) * k_1_2A;
       if ( (intrs_l > 0.0) && (intrs_l <= l) )
       {
         for (size_t i = 0; i < 3; i++)
           intrs_p[i] = a[i] + u[i] * intrs_l;
         intrs_mu = intrs_p[0] * cone->u[0] + 
                    intrs_p[1] * cone->u[1] +
                    intrs_p[2] * cone->u[2];
         if ( (intrs_mu > 0.0) && (intrs_p[2] > cone->h[0]) && 
           (intrs_p[2] < cone->h[1]) )
           return 1;
       }
     }
   }

   // Intersection with horizontal caps:
   // If no intersections occurred at the side surface, it might still occur 
   // with the caps, if the cone is semi-closed or closed. The initial point of
   // the ray is first re-centered on the base and on the top. Since this is a 
   // right-circular cone, X and Y positions do not change as the ray
   // is already centered and rotated to the cone axis.
   if ( cone->closed_f[0] )
   {
     double a_base[3] = {a[0], a[1], a[2] - cone->h[0]};
     if ( intrs_ellp(a_base, u, l, cone->base) )
       return 1;
   }

   if ( cone->closed_f[1] )
   {
     double a_top[3] = {a[0], a[1], a[2] - cone->h[1]};
     return intrs_ellp(a_top, u, l, cone->top);
   }

   return 0;
 }

 int
 intrs_cone_oblique
 (
   double const * a, 
   double const * b, 
   double const * u, 
   double const l, 
   struct str_cone const * cone
 )
 {
   // NOT IMPLEMENTED YET
   return 1;
 }

