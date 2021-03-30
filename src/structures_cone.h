
 #ifndef STRUCTURES_CONE
 #define STRUCTURES_CONE

 #include <stdio.h>		// FILE

 #include "structures.h"	// str_ellp

/*******************************************************************************
 cone  (str_cone)

 The structure below implements a positive (in)finite cone. The cone is 
 described by the position of its vertex, the orientation of its axis, the 
 positive heights, relative to vertex, along the axis and by the cosine of it
 half angle (i.e., the angle between the axis and the surface).

 For example, an infinite positive cone has height[0] = 0 and 
 height[1] = INFINITY. A truncated infinite cone has height[0] > 0 and a 
 truncated finite cone has height[0] > 0 and height[1] < INFINITY.

 Two additional descriptions are necessary to caracterize (and describe) if 
 the cone is open, semi-closed (one end) or closed (both ends).

 With the convention used for intersections in this code, the cone parameters 
 are retained in its own reference frame and the ray is centered and rotated if
 necessary. Therefore, all the relevant information becomes the vertex, its half
 angle, the heights, rotation flag and the rotation matrix. Additionally, the 
 properties of the caps are stored as ellipses and transformations are relative
 to the cone axis. Those are all calculated with a setup function.
         
         . Origin (vertex)   
                    | Cone axis
                    |
      _______       |_Height[0]
     /       \      |             
    /         \     |
   /           \    |
  /             \   |
 /_______________\  |_Height[1]
  
         |\                    
         | \ Psi

 The function setup_cone describes the parameters that must be passed to 
 create a cone structure.

 COMPONENTS:
 o        - The origin (vertex) of the cone
            Dimensions: [3]: [0]X, [1]Y, [2]Z
            Range: (-Inf,Inf), meters
            Array of double;
 h        - Height of the cone along its axis
            Dimensions: [2]: [0]Base, [1]Top
            Range: (0,Inf), meters
            Array of double;
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
 M        - Rotation matrix from the system's reference frame to the cone's
            reference frame
            Dimensions: [3][3]
            Range: [-1,1], unitless
            Pointer to pointer to double;
 closed_f - Flag to indicate if the cone is closed at the base and top opnenings
            Dimensions: [2]: [0]Base, [1]Top
            Values: 0 (FALSE), 1 (TRUE)
            Array of int;
 base     - Ellipse of the base opening
            Pointer to str_ellp;
 top      - Ellipse of the top opening
            Pointer to str_ellp;
 u        - Cartesian directions of the axis in the system's reference frame
            Dimensions: [3]: [0]X, [1]Y, [2]Z
            Range: [-1,1], unitless
            Array of double;
 s        - Spherical directions of the axis in the system's reference frame
            Dimensions: [3]: [0]Theta, [1]Phi, [2]Radius (always 1.0)
            Range: [0,PI], [0,2PI], radians
            Array of double;
 psi      - Cone half-angle
            Range: (0,PI), radians
            Double;
 mu       - Cosine of the cone's half-angle
            Range: (0,1), unitless
            Double;
 musq     - Square of the cosine of the cone's half-angle
            Range: (0,1), unitless
            Double;
 alp      - The azimuthal angle in the cones's reference frame
            Range: [0,2PI], radians
            Double;
 s_base   - Spherical directions of the normal (Z-axis) of the base opening,
            relative to the cone axis
            Dimensions: [3]: [0]Theta, [1]Phi, [2]Radius (always 1.0)
            Range: [0,PI], [0,2PI], radians
            Array of double;
 s_top    - Spherical directions of the normal (Z-axis) of the top opening,
            relative to the cone axis
            Dimensions: [3]: [0]Theta, [1]Phi, [2]Radius (always 1.0)
            Range: [0,PI], [0,2PI], radians
            Array of double;

 The cone is setup by the str_cone_setup function with values that can be read
 from the standard input file with the function str_cone_read.

 STANDARD INPUT FILE SPECIFICATION:
 O[X] O[Y] O[Z] S[t] S[p] Alpha Psi Height[Base] Height[top] S_base[t] S_base[p]
   S_top[t] S_top[p] Closed[Base] Closed[Top]

 The function str_cone_setup describes the parameters that must be passed to 
 create a cone structure.

*******************************************************************************/

 struct str_cone
 {
   double o[3];
   double h[2];
   double dhmx[2];
   int rotate_f[4];
   double ** M;
   int closed_f[2];
   str_ellp * base;
   str_ellp * top;
   double mu;
   double musq;
   double psi;
   double u[3];
   double s[2];
   double alp;
   double s_base[2];
   double s_top[2];
 };

/* Function prototypes: *******************************************************/

 struct str_cone *
 str_cone_alloc
 ( void );

 void
 str_cone_free
 ( struct str_cone ** cone );

 void
 str_cone_read
 (
   FILE * fi,
   long int * fpos,
   double * origin,
   double * axis,
   double * alpha,
   double * psi,
   double * height,
   double * base_s,
   double * top_s,
   int * closed
 );

 void
 str_cone_setup
 (
   struct str_cone * cone,
   double * origin,
   double const * axis,
   double const alpha,
   double const psi,
   double const * height,
   double const * s_base,
   double const * s_top,
   int const * closed
 );

 void
 str_cone_fprintf 
 (
   FILE * odv,
   struct str_cone const * cone,
   int const indent
 );

 void
 str_cones_fprintf
 (
   FILE * odv,
   struct str_cone const ** cones,
   int n,
   int const indent
 );

 #endif // STRUCTURES_CONE

