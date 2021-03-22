
 #ifndef STRUCTURE
 #define STRUCTURE

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

/*******************************************************************************
 Cuboid (str_cubd)

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

 typedef struct
 {
   double o[3];
   double hl[3];
   int    rotate_f[4];
   double **M;
   int    closed_f[2];
   str_rect *base;
   str_rect *top;
   double x[2];
   double y[2];
   double z[2];
   double u[3];
   double s[2];
   double alp;
   double s_base[2];
   double s_top[2];
 } str_cubd;

/* Function prototypes: *******************************************************/

 str_cubd* str_cubd_alloc (void);

 void str_cubd_free (str_cubd **cubd);

 void str_cubd_read (FILE *fi, long int *fpos, double *origin, double *axis,
   double *alpha, double *lengths, double *base_s, double *top_s, int *closed);

 void str_cubd_setup (str_cubd *cubd, double *origin, double const *axis,
   double const alpha, double const *lengths, double const *s_base, 
   double const *s_top, int const *closed);

 void
 str_cubd_fprintf 
 (
   FILE *odv,
   str_cubd const *cubd,
   int const indent 
 );

 void str_cubds_printf (str_cubd const **cubds, size_t n);

/******************************************************************************* 
 Cylinder (str_cyln)

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
 o        - Pointer to double (vector) containing the origin, the position of 
            the center of the cylinder's base, meters, (-INFINITY,INFINITY), 
            [0] X, [1] Y, [2] Z;
 rsq      - Pointer to double (vector) with the square of the radius of the
            cylinder, meters^2, (INFINITY, 0).
 h        - Height of the cylinder along its axis, meters, (0,INFINITY);
 dhmx     - Maximum change in height, relative to the center of the cap, due to 
            cap's inclination, meters, (0,sin(60º));
 rotate_f - Pointer to int (vector) containing flags to indicate if ray needs to
            be rotated before intersection check.
            [0]Cylinder, [1]Base, [2]Top, [3]Base to top;
 M        - Pointer to pointer to double (matrix) with the rotation matrix from 
            the system's reference frame to the cylinder's reference frame.
            [3][3];
 closed_f - Pointer to int (vector) containing flags to indicate if cylinder's
            base and top are open (0) or closed (1),
            [0]Base, [1]Top;
 r        - Pointer to double (vector) with the radius of the cylinder, meters,
            (0, INFINITY);
 u        - Pointer to double (vector) containing the Cartesian coordinates of 
            the normal (Z-axis), unitless, [-1, 1],
            [0] X, [1] Y, [2] Z;
 s        - Pointer to double (vector) containing the spherical coordinates of 
            the normal (Z-axis), radians, [0,PI] or [0,2PI],
            [0] Theta [1] Phi [2] Radius (always 1.0 and never used);
 alp      - The azimuthal angle in the cylinder's reference frame.
 s_base   - Pointer to double (vector) containing the spherical coordinates of 
            the normal (Z-axis) of the base opening, radians, [0,60º] and 0º or
            180º,
            [0]Theta [1]Phi;
 s_top    - Pointer to double (vector) containing the spherical coordinates of 
            the normal (Z-axis) of the base opening, radians, [0,60º] and 0º or
            180º,
            [0]Theta [1]Phi;

 The cylinder is setup by the str_cyln_setup function with values that can be
 read from the standard input file with the function str_cyln_read.

 STANDARD INPUT FILE SPECIFICATION:
 O[X] O[Y] O[Z] S[t] S[p] Alpha Radius Height S_base[t] S_base[p] S_top[t] 
   S_top[p] Closed[Base] Closed[Top]

 The function str_cyln_setup describes the parameters that must be passed to 
 create a cylinder structure.

*******************************************************************************/

 typedef struct
 {
   double o[3];
   double rsq;
   double h;
   double dhmx[2];
   int    rotate_f[4];
   double **M;
   int    closed_f[2];
   str_ellp *base;
   str_ellp *top;
   double r;
   double u[3];
   double s[2];
   double alp;
   double s_base[2];
   double s_top[2];
 } str_cyln;

/* Function prototypes: *******************************************************/

 str_cyln* str_cyln_alloc (void);

 void str_cyln_free (str_cyln **cyln);

 void str_cyln_read (FILE *fi, long int *fpos, double *origin, double *axis,
   double *alpha, double *radius, double *height, double *base_s, double *top_s,
   int *closed);

 void str_cyln_setup (str_cyln *cyln, double *origin, double const *axis,
   double const alpha, double const radius, double const height,
   double const *s_base, double const *s_top, int const *closed);

 void
 str_cyln_fprintf 
 (
   FILE *odv,
   str_cyln const *cyln,
   int const indent
 );

 void str_cylns_printf (str_cyln const **cylns, size_t n);

/*******************************************************************************
 Cone  (str_cone)

         
         . Origin    
                    | Cone axis
                    |
      _______       |_Height[0]
     /       \      |             
    /         \     |
   /           \    |
  /             \   |
 /_______________\  |_Height[1]
  
         |\                    
         | \ Theta

 The function setup_cone describes the parameters that must be passed to 
 create a cone structure.

 Cone
 *
 * The structure below implements a positive (in)finite cone. The cone is 
 * described by the position of its vertex, the orientation of its axis, the 
 * positive heights, relative to vertex, along the axis and by the cosine of it
 * half angle (i.e., the angle between the axis and the surface.
 *
 * For example, an infinite positive cone has height[0] = 0 and 
 * height[1] = INFINITY. A truncated infinite cone has height[0] > 0 and a 
 * truncated finite cone has height[0] > 0 and height[1] < INFINITY.
 *
 * Two additional descriptions are necessary to caracterize (and describe) if 
 * the cone is open, semi-closed (one end) or closed (both ends).
 * 
 

*******************************************************************************/

 struct str_cone
 {
   double vertex[3];	// Position of the vertex. [0] X, [1] Y, [2] Z
   double paxis[2];	// Polar directions of the axis. [0] theta [1] phi
   double caxis[3];     // Direction cosines of the axis. [0] X, [1] Y, [2] Z
   double psi;          // Cone's half angle;
   double mu;           // Cosine of the cone's half angle (angle between axis and side)
   double mu_sq;        // Square of mu. 
   double height[2];    // Height along the axis, from the vertex. [0] min, [1] max
   double radius[2];    // Radius at the minimum height [0] and at the maximum height [1].
   int    closed[2];    // Logical flag to indicate if its closed in the minimum height [0] and maximum height [1].
   int    rotate_f;	// Logical: Should the ray direction be rotated?
 };

 typedef struct
 {
   int a;
 } str_cone;

/* Function prototypes: *******************************************************/

 #endif // STRUCTURE

