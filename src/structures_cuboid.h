
 #ifndef STRUCTURES_CUBOID
 #define STRUCTURES_CUBOID

 #include <stdio.h>		// FILE

/*******************************************************************************
 cuboid (str_cubd)

 NOTE: In the current version, parallelepipeds and open cuboids are not 
 implemented.

 This structure implements a general cuboid with arbitrary orientation axis 
 of the cuboid' normal and rotation around the normal. The cuboid is defined by 
 its origin, length of its sides, azimuth rotation in its own reference frame 
 (alpha), and orientation of the Z-axis (normal) in the system's reference 
 frame.

 Additionally, the cuboid may be open or closed in the top and bottom sides. 
 Those "caps" are represented by rectangles. The inclination with regard to the
 cuboid axis can be any value between 0 and 60º. Therefore it is possible to 
 have an oblique cuboid (parallelepiped). The other limit of the caps is that 
 they must have the same alpha as the cuboid.

 With the convention used for intersections in this code, the cuboid 
 information is retained in its own reference frame and the ray is centered and
 rotated if necessary. Therefore, all the relevant information becomes the
 origin, the half-length of each axis, the rotation flag and the rotation 
 matrix. Additionally, the properties of the caps are stored as rectangles and
 transformations are relative to the cuboid axis. Those are all calculated with
 a setup function.

                 _____________
                /     Origin /|        Alpha
               /     ·      / |        -------> 
              /____________/  |       |
 Z_length  |  |            |  |       | Cuboid axis   
           |  |            |  |       |
           |  |            |  |       |             
           |  |            |  /   /   |
           |  |            | /   /    
           |  |____________|/   / Y_length
                      
               ____________
               X_length                

 COMPONENTS:
 o        - Pointer to double (vector) containing the origin, the position of  
            the center of the rectangle, meters, (-INFINITY,INFINITY), 
            [0]X, [1]Y, [2]Z;
 hl       - Half-length in each axis, meters, (-INFINITY,INFINITY),
            [0]X, [1]Y, [2]Z;
 rotate_f - Pointer to int (vector) containing flags to indicate if ray needs to
            be rotated before intersection check.
            [0]Cuboid, [1]Base, [2]Top, [3]Base to top;
 M        - Pointer to pointer to double (matrix) with the rotation matrix from 
            the system's reference frame to the cuboid' reference frame.
            [3][3];
 closed_f - Pointer to int (vector) containing flags to indicate if cuboid base
            and top are open (0) or closed (1),
            [0]Base, [1]Top;
 base     - Pointer to rectangle representing the base; 
 top      - Pointer to rectangle representing the top;
 x        - Pointer to double (vector) with X-axis positions in it own reference 
            frame, meters, (-INFINITY,INFINITY), 
            [0] Lower, [1] Upper;
 y        - Pointer to double (vector) with Y-axis positions in it own reference 
            frame Y-axis, meters, (-INFINITY,INFINITY),
            [0] Lower, [1] Upper;
 z        - Pointer to double (vector) with Z-axis positions in it own reference 
            frame, meters, (-INFINITY,INFINITY), 
            [0] Lower, [1] Upper;
 u        - Pointer to double (vector) containing the Cartesian coordinates of 
            the normal (Z-axis), unitless, [-1, 1],
            [0] X, [1] Y, [2] Z;
 s        - Pointer to double (vector) containing the spherical coordinates of 
            the normal (Z-axis), radians, [0,PI] or [0,2PI],
            [0] Theta [1] Phi [2] Radius (always 1.0 and never used);
 alp      - The azimuthal angle of rotation in the rectangle's reference frame,
            radians, [0,2PI];
 s_base   - Pointer to double (vector) containing the spherical coordinates of 
            the normal (Z-axis) of the base opening, radians, [0,60º] and 0º or
            180º,
            [0]Theta [1]Phi;
 s_top    - Pointer to double (vector) containing the spherical coordinates of 
            the normal (Z-axis) of the base opening, radians, [0,60º] and 0º or
            180º,
            [0]Theta [1]Phi;

 The cuboid is setup by the str_cubd_setup function, requiring only o, x, y, 
 s and alp to be provided. Those values can be read from the standard input file 
 with the function str_cubd_read.

 STANDARD INPUT FILE SPECIFICATION:
 O[X] O[Y] O[Z] S[t] S[p] Alpha L[X] L[Y] L[Z] S_base[t] S_base[p] S_top[t] 
   S_top[p] Closed[Base] Closed[Top]

 The function str_cubd_setup describes the parameters that must be passed to 
 create a cuboid structure.

*******************************************************************************/

 struct str_cubd
 {
   double o[3];
   double hl[3];
   int rotate_f[4];
   double ** M;
   int closed_f[2];
   str_rect * base;
   str_rect * top;
   double x[2];
   double y[2];
   double z[2];
   double u[3];
   double s[3];
   double alp;
   double s_base[2];
   double s_top[2];
 };

/* Function prototypes: *******************************************************/

 struct str_cubd *
 str_cubd_alloc
 ( void );

 void
 str_cubd_free
 ( struct str_cubd ** cubd );

 void
 str_cubd_read
 (
   FILE * fi,
   long int * fpos,
   double * origin,
   double * axis,
   double * alpha,
   double * lengths,
   double * base_s,
   double * top_s,
   int * closed
 );

 void
 str_cubd_setup
 (
   struct str_cubd * cubd,
   double * origin,
   double const * axis,
   double const alpha,
   double const * lengths,
   double const * s_base, 
   double const * s_top,
   int const * closed
 );

 void
 str_cubd_fprintf 
 (
   FILE * odv,
   struct str_cubd const * cubd,
   int const indent 
 );

 void
 str_cubds_fprintf
 (
   FILE * odv,
   struct str_cubd const ** cubds,
   int n,
   int const indent 
 );

 #endif // STRUCTURES_CUBOID

