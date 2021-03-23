
 #ifndef STRUCTURES
 #define STRUCTURES

 #include <stdio.h>
 #include <stddef.h>

/*******************************************************************************
 Rectangle (str_rect)
 
 This structure implements a general rectangle with arbitrary orientation axis 
 of the rectangle's normal, and rotation around that normal. The rectangle is
 defined by its origin, length of its sides, azimuth rotation in its own 
 reference frame (alpha), and orientation of the Z-axis (normal) in the system's
 reference frame.

 With the convention used for intersections in this code, the rectangle 
 information is retained in its own reference frame and the ray is centered and
 rotated if necessary. Therefore, all the relevant information becomes the
 origin, the half-length of each side, the rotation flag and the rotation 
 matrix. Those are calculated with a setup function.


    ________    Z
    \    o  \   |
     \   ·   \  |------> 
      \_______\  \  alpha
                  \
                   X
 
 COMPONENTS:
 o   - Pointer to double (vector) containing the origin, the position of the 
       center of the rectangle, meters, (-INFINITY,INFINITY), 
       [0] X, [1] Y, [2] Z;
 x   - Pointer to double (vector) with X-axis positions in it own reference 
       frame, meters, (-INFINITY,INFINITY), 
       [0] Lower, [1] Upper;
 y   - Pointer to double (vector) with Y-axis positions in it own reference 
       frame Y-axis, meters, (-INFINITY,INFINITY),
       [0] Lower, [1] Upper;
 M   - Pointer to pointer to double (matrix) with the rotation matrix from the 
       system's reference frame to the rectangle's reference frame.
       [3][3];
 u   - Pointer to double (vector) containing the Cartesian coordinates of the 
       normal (Z-axis), unitless, [-1, 1],
       [0] X, [1] Y, [2] Z;
 s   - Pointer to double (vector) containing the spherical coordinates of the 
       normal (Z-axis), radians, [0,PI] or [0,2PI],
       [0] Theta [1] Phi [2] Radius (always 1.0 and never used);
 alp - The azimuthal angle of rotation in the rectangle's reference frame,
       radians, [0,2PI];
 hl  - Half-length in X and Y-axis, meters, (-INFINITY,INFINITY),
       [0]X [1]Y;
 rotate_f - Flag to indicate if ray needs to be centered and rotated before 
            intersection check;

 The rectangle is setup by the str_rect_setup function, requiring only o, x, y, 
 s and alp to be provided. Those values can be read from the standard input file 
 with the function str_rect_read.

 STANDARD INPUT FILE SPECIFICATION:
 O[X] O[Y] O[Z] S[t] S[p] Alpha L[X] L[Y]

 Examples:

 RECTANGLE WITH NORMAL AT {60.00º,30.00º} AND ALPHA AT 10.00º:
 INPUT FILE:  0.0 0.0 0.0 60.0 30.0 10.0 12.0 3.0
 Rectangle parameters:
   Origin:    0.00e+00,  0.00e+00,  0.00e+00 (m)
   X-range:  -6.00e+00,  6.00e+00  (1.20e+01 m)
   Y-range:  -1.50e+00,  1.50e+00  (3.00e+00 m)
   Axis:       60.00º,  30.00º ( 7.50e-01,  4.33e-01,  5.00e-01)
   Alpha:      10.00º
   Rotate:    1

 SQUARE ALIGNED WITH Z-AXIS AND ALPHA AT 30.00º AND ORIGIN AT -5 M (X):
 INPUT FILE: -5.0 0.0 0.0 0.0 0.0 30.0 3.0 3.0
 Rectangle parameters:
   Origin:   -5.00e+00,  0.00e+00,  0.00e+00 (m)
   X-range:  -1.50e+00,  1.50e+00  (3.00e+00 m)
   Y-range:  -1.50e+00,  1.50e+00  (3.00e+00 m)
   Axis:        0.00º,   0.00º ( 0.00e+00,  0.00e+00,  1.00e+00)
   Alpha:      30.00º
   Rotate:    1

*******************************************************************************/

 typedef struct
 {
   double o[3];
   int rotate_f;
   double **M;
   double hl[2];
   double x[2];
   double y[2];
   double s[3];
   double u[3];
   double alp;
 } str_rect;

/* Function prototypes: *******************************************************/

 str_rect* str_rect_alloc (void);

 void str_rect_free (str_rect **rect);

 void str_rect_read (FILE *fi, long int *fpos, double *origin, double *axis,
   double *alpha, double *lengths);

 void str_rect_setup (str_rect *rect, double const *origin, double const *axis,
   double const alpha, double const *lengths);

 void
 str_rect_fprintf 
 (
   FILE *odv,
   str_rect const *rect,
   int const indent 
 );

 void str_rects_printf (str_rect const **rects, size_t n);

/*******************************************************************************
 Triangle (str_trig)

 Not implemented.
         
        /\
       /  \
      /    \
     /      \
    /        \
   /          \
  /            \
 /______________\
  

*******************************************************************************/

/* Function prototypes: *******************************************************/

/*******************************************************************************
 Ellipse (str_ellp)

 This structure implements a general ellipse with arbitrary orientation axis 
 of the ellipse's normal and rotation around the normal. The ellipse is
 defined by its origin, radius of semi-major and minor axis, azimuth rotation 
 in its own reference frame (alpha), and orientation of the Z-axis (normal) in 
 the system's reference frame.

 With the convention used for intersections in this code, the ellipse 
 information is retained in its own reference frame and the ray is centered and
 rotated if necessary. Therefore, all the relevant information becomes the
 origin, the inverse of the squared minor and major radius, the rotation flag 
 and the rotation matrix. Those are calculated with a setup function.

 COMPONENTS:
 o        - Pointer to double (vector) containing the origin, the position of 
            the center of the rectangle, meters, (-INFINITY,INFINITY), 
            [0] X, [1] Y, [2] Z;
 rotate_f - Flag to indicate if ray needs to be centered and rotated before 
            intersection check;
 rsq_inv  - Pointer to double (vector) with the square of the minor and major 
            radius of the ellipse, meters^2, (-INFINITY,INFINITY), 
            [0] X, [1] Y.
 M        - Pointer to pointer to double (matrix) with the rotation matrix from
            the system's reference frame to the rectangle's reference frame.
            [3][3];
 r        - Pointer to double (vector) with minor and major radius of the 
            ellipse, meters, (-INFINITY,INFINITY), 
            [0] X, [1] Y.
 u        - Pointer to double (vector) containing the Cartesian coordinates of 
            the normal (Z-axis), unitless, [-1, 1],
            [0] X, [1] Y, [2] Z;
 s        - Pointer to double (vector) containing the spherical coordinates of 
            the normal (Z-axis), radians, [0,PI] or [0,2PI],
            [0] Theta [1] Phi [2] Radius (always 1.0 and never used);
 alp      - The azimuthal angle in the ellipse's reference frame.

 The ellipse is setup by the str_ellp_setup function, requiring only o, r, s 
 and alp to be provided. Those values can be read from the standard input file 
 with the function str_ellp_read.

 STANDARD INPUT FILE SPECIFICATION:
 O[X] O[Y] O[Z] S[t] S[p] Alpha Radius[X] Radius[Y]

 Examples:

 CIRCLE WITH NORMAL AT {60.00º,30.00º} AND ALPHA AT 10.00º:
 Ellipse parameters:
   Origin:    0.00e+00,  0.00e+00,  0.00e+00 (m)
   Radius:   1.00e+00,  1.00e+00  (m)
   Axis:       60.00º,  30.00º ( 7.50e-01,  4.33e-01,  5.00e-01)
   Alpha:      10.00º
   Rotate:    1
   Rotation matrix:
   Dimension: 3, 3
   Data:  3.396e-01  3.966e-01 -8.529e-01 
         -5.676e-01  8.095e-01  1.504e-01 
          7.500e-01  4.330e-01  5.000e-01 

 ELLIPSE ALIGNED WITH Z-AXIS AND ALPHA AT 30.00º AND ORIGIN AT -3 M (X):
 Ellipse parameters:
   Origin:   -3.00e+00,  0.00e+00,  0.00e+00 (m)
   Radius:   3.00e+00,  1.00e+00  (m)
   Axis:        0.00º,   0.00º ( 0.00e+00,  0.00e+00,  1.00e+00)
   Alpha:      30.00º
   Rotate:    1
   Rotation matrix:
   Dimension: 3, 3
   Data:  8.660e-01  5.000e-01  0.000e+00 
         -5.000e-01  8.660e-01  0.000e+00 
          0.000e+00  0.000e+00  1.000e+00 

*******************************************************************************/

 typedef struct
 {
   double o[3];
   double rsq_inv[2];
   int rotate_f;
   double **M;
   double r[2];
   double u[3];
   double s[3];
   double alp;
 } str_ellp;

/* Function prototypes: *******************************************************/

 str_ellp* str_ellp_alloc (void);

 void str_ellp_free (str_ellp **ellp);

 void str_ellp_read (FILE *fi, long int *fpos, double *origin, double *axis,
   double *alpha, double *radius);

 void str_ellp_setup (str_ellp *ellp, double const *origin, double const *axis,
   double const alpha, double const *radius);

 void
 str_ellp_fprintf 
 (
   FILE *odv,
   str_ellp const *ellp,
   int const indent 
 );

 void str_ellps_printf (str_ellp const **ellp, size_t n);


/* Header for three-dimensional types: ****************************************/

 #include "structures_cylinder.h"
 #include "structures_cone.h"
 #include "structures_cuboid.h"

 #endif // STRUCTURES

