
/*******************************************************************************
 intersect_cylinder.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 Version: 2021-03-12

 intrs_cyln, intrs_cyln_right, intrs_cyln_oblique

 This source file provides functions to check for intersection of a ray against
 a circular cylinder (right or oblique).

 The specific functions assume that the ray's initial and terminal points are 
 centered on the cylinder's origin and that if the cylinder orientation axis
 does not have the same direction as the system's Z-axis, the ray points and 
 direction were rotated to the cylinder's reference frame. Those transformations
 are performed by the wrapper function intrs_cyln which also calls the
 appropriate function depending on cylinder type. The derivation below considers
 the case of a right-circular cylinder. Changes to an oblique-circular cylinder
 are considered in the end.

 In the following, all quantities in upper case are vector quantities in three 
 dimensions (x, y, z). "P" is the vector coordinates of an arbitrary point, "A"
 is the starting position of the ray, "U" is the vector of Cartesian directions
 of the ray, and "l" is the (scalar) pathlength of the ray. For the cylinder,
 "O" is the center of the base of the cylinder, "r" is its radius and "N" is a
 unit vector giving its orientation (axis) and "h" is the cylinder's height
 given as a distance along the cylinder axis. 
 
 The line equation is defined as:
 P = U * l + A.                                                              (1)
 
 The (infinitely tall) cylinder's lateral surface equation is given by:
 P = O + r * {cos(phi), sin(phi), h/r},                                      (2)
 
 and any point on the infinite right-circular cylinder lateral surface must 
 satisfy:
 ||P - O||^2 - ((P - O) · N)^2 = r^2,                                        (3)

 where O is the center of the base of the cylinder, N is the orientation axis of
 the cylinder, || || is the l2 norm and · is the dot (scalar) product.

 Intersections with the lateral surface of the cylinder have a P satisfying
 Eq. 1 and Eq. 3. Substituting Eq. 1 into Eq. 3 we have: 
 ||P - O||^2 - ((P - O) · N)^2 = r^2 ->
 ||U * l + A - O||^2 - ((U * l + A - O) · N)^2 = r^2 ->
 ||U * l + (A - O)||^2 = ((U * l + (A - O)) · N)^2 + r^2.                    (4)

 Expanding the left-hand side we have:
 ||U * l + (A - O)||^2 ->
 ((U * l + (A - O)) · (U * l + (A - O))) ->
 (l^2 * (U · U) + 2 * l * (U · (A - O)) + ((A - O) · (A - O))^2) ->
 l^2 + 2 * l * (U · (A - O)) + ((A - O) · (A - O))^2,                        (5)
 
 where since U is a unit vector, U · U = 1 and is removed. Now we can also 
 expand the right hand side of Eq. 4:
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

 An intersection of the ray with the side of the cylinder occurs if l is 
 positive and lower than the total pathlength of the ray. If a = 0, then b = 0
 and the ray's direction is parallel to the cylinder's axis. In this case an
 intersection with the infinite cylinder occurs only if the ray is on the
 surface of the cylinder, i.e., if c = 0. The intersection point(s) (or range)
 is given by substituting l from Eq. 8 in Eq. 1. If (P - O) · N < than h, the
 height of the finite cylinder, the vector intersects the finite cylinder.
 
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
 azimuth of 180.0º. The same is true for the top cap. The change is 0 at 90.0º.
 The height change is +- r_x * sin(theta) so:

 DeltaH = r_x * -sin(theta) * P[x] / r_x ->
 DeltaH = -sin(theta) * P[x]                                                (11)

 where r_x is the major axis (here equal to radius in X) of the ellipse formed 
 by the opening, theta is the angle of the cap's normal relative to the cylinder
 axis. Therefore, before the check to see if intersection occurs within the 
 height of the finite cylinder, the height is adjusted by DeltaH. Further, the
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
 cyln - Cylinder parameters
        Pointer to str_cyln.

 OUTPUT:
 Flag indicating if an intersection occurred: (0) No, (1) Yes.

*******************************************************************************/

 #include <stddef.h>		// size_t
 #include <math.h>		// sqrt

 #include "config.h"		// TOLERANCE
 #include "aux.h"		// ABS
 #include "rotation.h"		// ROT_VEC, ROT_VEC_UNIT
 #include "structures.h"	// str_cyln
 #include "intersect.h"		// intrs_ellp

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
   double a_cr[3];	// Ray's initial point, centered and rotated
   double b_cr[3];	// Ray's terminal point, centered and rotated 
   double u_rot[3];	// Ray's Cartesian directions, rotated
   double c_vec[3];	// Temporary variable for the macro ROT_VEC

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
   double const * a, 
   double const * b, 
   double const * u, 
   double const l, 
   str_cyln const * cyln
 )
 {
   double intrs_p[3];		// Position of intersection

   double const A = 1.0 - u[2] * u[2];
   double const C = a[0] * a[0] + a[1] * a[1] - cyln->rsq;

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
     double const B     = 2.0 * ( (u[0] * a[0]) + (u[1] * a[1]) );
     double const delta = (B * B) - (4.0 * A * C);

     if ( delta >= 0.0 )
     {
       // If delta >= 0, the equality has real roots and the line intersects the 
       // infinite cylinder surface. To know if the ray intersects the cylinder, 
       // we check for each root, if it is positive, if it smaller than the 
       // pathlength of the ray and if the intersection vertical position is 
       // within the vertical extent of the finite cylinder.
       double const sdelta = sqrt(delta);
       double const k_1_2A = 1.0 / (2.0 * A);

       double intrs_l = (-B - sdelta) * k_1_2A;
       if ( (intrs_l > 0.0) && (intrs_l <= l) )
       {
         intrs_p[2] = a[2] + u[2] * intrs_l;

         if ( (intrs_p[2] > 0.0) && (intrs_p[2] <= cyln->h) )
         {
           return 1;
         }
       }

       intrs_l = (-B + sdelta) * k_1_2A;
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
   // is a right circular cylinder, X and Y positions do not change as the ray
   // is already rotated to the cylinder axis.
   if ( cyln->closed_f[0] )
   {
     if ( intrs_ellp(a, u, l, cyln->base) )
     {
       return 1;
     }
   }
   if ( cyln->closed_f[1] )
   {
     double a_top[3] = {a[0], a[1], a[2] - cyln->h};
     return intrs_ellp(a_top, u, l, cyln->top);
   }
   return 0;
 }


 int
 intrs_cyln_oblique
 (
   double const * a, 
   double const * b, 
   double const * u, 
   double const l, 
   str_cyln const * cyln
 )
 {
   double dh[2];		// Delta H, change in height due to cap inclination
   double intrs_dhf;		// Fractional Delta H based on azimuth for each opening

   double const A = 1.0 - u[2] * u[2];
   double const C = a[0] * a[0] + a[1] * a[1] - cyln->rsq;

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
     double const B     = 2.0 * ( (u[0] * a[0]) + (u[1] * a[1]) );
     double const delta = (B * B) - (4.0 * A * C);
     if ( delta >= 0.0 )
     {
       // If delta >= 0, the equality has real roots and the line intersects the 
       // infinite cylinder surface. To know if the ray intersects the cylinder, 
       // we check for each root, if it is positive, if it smaller than the 
       // pathlength of the ray and if the intersection vertical position is 
       // within the vertical extent of the finite cylinder.
       double const sdelta = sqrt(delta);
       double const k_1_2A = 1.0 / (2.0 * A);
       double intrs_p[3];		// Position of intersection

       double intrs_l = (-B - sdelta) * k_1_2A;
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

       intrs_l = (-B + sdelta) * k_1_2A;
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
   double a_rot[3];		// Vector of position rotated for the caps
   double u_rot[3];		// Vector of directions rotated for the caps

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
         double c_vec[3];	// Temporary variable for the ROT_VEC macro
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


