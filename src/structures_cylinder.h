
 #ifndef STRUCTURES_CYLINDER
 #define STRUCTURES_CYLINDER

 #include <stdio.h>		// FILE

 #include "structures.h"	// str_ellp

/******************************************************************************* 
 cylinder (str_cyln)

 This structure implements a general positive (in)finite cylinder (right 
 circular or oblique) with arbitrary orientation axis of the cylinder's normal 
 and rotation around the normal. The cylinder is defined by its origin, radius, 
 height along its axis, azimuth rotation in its own reference frame (alpha), and
 orientation of the Z-axis (normal) in the system's reference frame.

 Additionally, the cylinder may be open or closed in the top and bottom sides. 
 Those "caps" are represented by ellipses. The inclination with regard to the
 cuboid axis can be any value between 0 and 60º. Therefore it is possible to 
 have an oblique-circular cylinder. The other limit of the caps is that they 
 must have the same alpha as the cylinder.

 With the convention used for intersections in this code, the cylinder 
 information is retained in its own reference frame and the ray is centered and
 rotated if necessary. Therefore, all the relevant information becomes the
 origin, the radius, the height, rotation flag and the rotation matrix. 
 Additionally, the properties of the caps are stored as ellipses and 
 transformations are relative to the cylinder axis. Those are all calculated
 with a setup function.

          Radius
          ____

         Origin
      ___.___
     /      /   / Cylinder axis
    /      /   /
   /      /   /
  /      /   /
 /______/   /_Height[1]
 

 COMPONENTS:
 o        - The origin (center of the base) of the cylinder
            Dimensions: [3]: [0]X, [1]Y, [2]Z
            Range: (-Inf,Inf), meters
            Array of double;
 rsq      - Squared radius
            Range: (0,Inf), meters^2
            Double;
 h        - Height of the cylinder along its axis
            Range: (0,Inf), meters
            Double;
 dhmx     - Maximum relative change in height, relative to the center of the
            cap, due to cap's inclination
            Dimensions: [2]: [0]Base, [1]Top
            Range: (0,sin(60º)), unitless
            Array of double;
 rotate_f - Flag to indicate if the ray needs to be rotated before intersection
            check
            Dimensions: [4]: [0]Cone, [1]Base, [2]Top, [3]Base to top
            Values: 0 (FALSE), 1 (TRUE)
            Array of int;
 M        - Rotation matrix from the system's reference frame to the cylinder's
            reference frame
            Dimensions: [3][3]
            Range: [-1,1], unitless
            Pointer to pointer to double;
 closed_f - Flag to indicate if the cylinder is closed at the base and top
            openings
            Dimensions: [2]: [0]Base, [1]Top
            Values: 0 (FALSE), 1 (TRUE)
            Array of int;
 base     - Ellipse of the base opening
            Pointer to str_ellp;
 top      - Ellipse of the top opening
            Pointer to str_ellp;
 r        - Radius of the circular cylinder
            Range: (0,Inf), meters
            Double;
 u        - Cartesian directions of the axis in the system's reference frame
            Dimensions: [3]: [0]X, [1]Y, [2]Z
            Range: [-1,1], unitless
            Array of double;
 s        - Spherical directions of the axis in the system's reference frame
            Dimensions: [3]: [0]Theta, [1]Phi, [2]Radius (always 1.0)
            Range: [0,PI], [0,2PI], radians
            Array of double;
 alp      - The azimuthal angle in the cylinder's reference frame
            Range: [0,2PI], radians
            Double;
 s_base   - Spherical directions of the normal (Z-axis) of the base opening,
            relative to the cylinder axis
            Dimensions: [3]: [0]Theta, [1]Phi, [2]Radius (always 1.0)
            Range: [0,PI], [0,2PI], radians
            Array of double;
 s_top    - Spherical directions of the normal (Z-axis) of the top opening,
            relative to the cylinder axis
            Dimensions: [3]: [0]Theta, [1]Phi, [2]Radius (always 1.0)
            Range: [0,PI], [0,2PI], radians
            Array of double;

 The cylinder is setup by the str_cyln_setup function with values that can be
 read from the standard input file with the function str_cyln_read.

 STANDARD INPUT FILE SPECIFICATION:
 O[X] O[Y] O[Z] S[t] S[p] Alpha Radius Height S_base[t] S_base[p] S_top[t] 
   S_top[p] Closed[Base] Closed[Top]

 The function str_cyln_setup describes the parameters that must be passed to 
 create a cylinder structure.

*******************************************************************************/

 struct str_cyln
 {
   double o[3];
   double rsq;
   double h;
   double dhmx[2];
   int rotate_f[4];
   double ** M;
   int closed_f[2];
   str_ellp * base;
   str_ellp * top;
   double r;
   double u[3];
   double s[2];
   double alp;
   double s_base[2];
   double s_top[2];
 };

/* Function prototypes: *******************************************************/

 struct str_cyln *
 str_cyln_alloc
 ( void );

 void
 str_cyln_free
 ( struct str_cyln ** cyln );

 void
 str_cyln_read
 (
   FILE * fi,
   long int * fpos,
   double * origin,
   double * axis,
   double * alpha,
   double * radius,
   double * height,
   double * base_s,
   double * top_s,
   int * closed
 );

 void
 str_cyln_setup
 (
   struct str_cyln * cyln,
   double * origin,
   double const * axis,
   double const alpha,
   double const radius,
   double const height,
   double const * s_base,
   double const * s_top,
   int const * closed
 );

 void
 str_cyln_fprintf 
 (
   FILE * odv,
   struct str_cyln const * cyln,
   int const indent
 );

 void
 str_cylns_fprintf
 (
   FILE * odv,
   struct str_cyln const ** cylns, 
   int n,
   int const indent
 );


 #endif // STRUCTURES_CYLINDER

