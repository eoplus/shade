
/*******************************************************************************
 rotation.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 General function to compute rotation matrices and apply then to ray direction
 and Stokes vector.

 Rotations in two dimensions are simple to describe. We take a virtual Z axis
 orthogonal to the XY plane. When looking into the Z-axis, the active rotation
 counterclockwise around Z in a right-handed coordinate system is given by the
 pre-multiplication matrix:

 R(a) = | cos(a) -sin(a) |
        | sin(a)  cos(a) |, (a stands for angle alpha)                       (1)

 where a (alpha) is the rotation angle. The following represents a rotation of 
 45º for a vector laying over the X axis:

  Y   .
  |  /        
  | /         R(45) * | 1 | = | 0.7071068 |
  |/___. X            | 0 |   | 0.7071068 |,

 where ^T represents a transpose of the row vector to a column vector. This idea
 is generalized in three dimensions by rotating around each axis at a time. The 
 matrix R is therefore expanded to a 3x3 matrix, with the reference axis, around 
 which the rotation is taking place, receiving only 0 and 1. For example, the 
 above rotation in two dimensions used a "virtual" Z axis to explain the 
 concept, and in three dimensions the matrix is the same, but expanded:

 R[Z](alpha) = | cos(a) -sin(a)      0  |
               | sin(a)  cos(a)      0  |
               |     0       0       1  |. (alpha on XY plane)               (2)

 The same is true for the Y and X axis:

 R[Y](beta) =  |  cos(b)     0   sin(b) |
               |      0      1       0  |
               | -sin(b)     0   cos(b) |, (beta on XZ plane)                (3)

 R[X](gamma) = |      1      0       0  |
               |      0  cos(g) -sin(g) |
               |      0  sin(g)  cos(g) |, (gamma on YZ plane)               (4)

 It is common to combine those matrices through multiplication, resulting in a
 single rotation matrix around Z, Y and X:

 R(alpha,beta,gamma) = R[Z](alpha) * R[Y](beta) * R[X](gamma) =              (5)
| cos(a)cos(b) cos(a)sin(b)sin(g)-sin(a)cos(g) cos(a)sin(b)cos(g)+sin(a)sin(g) |
| sin(a)cos(b) sin(a)sin(b)sin(g)+cos(a)cos(g) sin(a)sin(b)cos(g)-cos(a)sin(g) |
|      -sin(b)                    cos(b)sin(g)                    cos(b)cos(g) | 
                                                                             
 The rotation angles alpha, beta and gamma are also known as row, yaw and pitch.

 One use of the rotation matrix is to rotate a direction from the scattering 
 reference frame to the system's reference frame. Since the ray's reference 
 frame has an undefined azimuth (ray is a vector aligned with its Z-axis), the
 only rotations that need to be specified are around the Z-axis of the system 
 (azimuth of travel) and around the Y-axis (polar angle of travel). In this case
 The rotation matrix simplifies to:

 R(phi,theta) = R[Z](alpha) * R[Y](beta) = 
 | cos(p)cos(t) -sin(p) cos(p)sin(t) |
 | sin(p)cos(t)  cos(p) sin(p)sin(t) |
 |      -sin(t)      0        cos(t) |,                                      (6)
 
 where theta and phi are the polar and azimuthal angle of the ray in the 
 system's reference frame. Note that the last column are the directional cosines 
 of the ray in X, Y and Z-axis, which is the orientation of the Z-axis of the 
 scattering reference frame. This simplified rotation means that the orientation
 of the X-axis of the ray's reference frame will have a fixed relation with the
 orientation of its Z-axis. In other words, the direction of the ray's X-axis 
 varies in a predictable manner with the direction of ray travel. This however 
 is not an issue in natural waters since the medium is macroscopically isotropic 
 and has no azimuthal dependency of scattering (for an assemble of particles 
 with random orientation). The matrix in Eq. 6 is the one provided by Leathers 
 et al. (2004).

 When testing intersections of a ray with a structure that has an arbitrary 
 orientation axis and rotation around that axis, one approach is to first rotate 
 the ray to the structure's reference frame. Here the azimuthal symmetry depends 
 on the azimuthal symmetry of the structure (e.g., a circle or a right-circular 
 cylinder) but the general case requires that the orientation of the X-axis of
 the structure reference frame is also tracked. However, Eq. 5 is not very 
 convenient to express orientation of a structure. The more natural way is 
 through the polar coordinates of the axis of the structure an a rotation around 
 that axis, which implies two rotations around the Z-axis:

 R(phi,theta,alpha) = R[Z](phi) * R[Y](theta) * R[Z](alpha) =                (7) 
| cos(p)cos(t)cos(a)-sin(p)sin(a) -cos(p)cos(t)sin(a)-sin(p)cos(a) cos(p)sin(t) |
| sin(p)cos(t)cos(a)+cos(p)sin(a) -sin(p)cos(t)sin(a)+cos(p)cos(a) sin(p)sin(t) |
|                   -sin(t)cos(a)                     sin(t)sin(a)       cos(t) |.

 As expected, the first column is the direction cosines of the structure's 
 X-axis, the second column is the Y-axis and the last column the Z-axis. Note
 that Eq. 7 is not equal to Eq. 6 with R(phi + alpha, theta). To project a ray
 into the reference frame of this structure, we multiply its direction and its
 origin-cenetered position by the inverse of Eq. 7. For those rotation matrices
 the inverse is equal to the transpose.

 Lets take the example of a cuboid with dimensions 5, 3 and 2 in X, Y and Z of
 its own reference frame. Lets define that this cuboid has its Z-azis with 
 polar coordinates theta = 30, phi = 45. Additionally, we define that this 
 cuboid is rotated around its own Z-axis by an angle of alpha = 90º, then we 
 have from Eq. 7:

 R(45,30,90)   = | -0.7071068 -0.6123724 0.3535534 |
                 |  0.7071068 -0.6123724 0.3535534 |
                 |          0  0.5000000 0.8660254 |,

 and

 R(45,30,90)^T = | -0.7071068  0.7071068         0 |
                 | -0.6123724 -0.6123724 0.5000000 |
                 |  0.3535534  0.3535534 0.8660254 |,


 which is the rotation matrix to rotate a ray to the structure's reference 
 frame. Finally, we note that rotations are about a reference point of the 
 structure, such that the initial (A) and terminal (B) points of the ray need to 
 be translated to the origin (O) point of the structure: A_s = A - O, 
 B_s = B - O, before the rotation is applied. Rotation matrix for structures are
 fixed and computed only once and stored. If a very large number of structures 
 is present this method might be inefficient memory wise.

 The rotation is performed with a matrix to column vector multiplication, as in
 the function rot_vec.

 A property of those rotation matrices is that the inverse matrix is equal to 
 the transpose of the rotation matrix and its determinant is equal to 1 (proper)
 or -1 (improper).

*******************************************************************************/

 #include <math.h>

 #include "config.h"
 #include "aux.h"
 #include "geometry.h"

/* rot_mat_ZYX

 Calculates the rotation matrix based on sequential Z, Y and X rotations (roll, 
 yaw, pitch).

 INPUT:
 s - Pointer to array containing the couterclockwise rotation angles around 
     the Z, Y and X-axis, radians, [0, 2 PI];
 M - Pointer to pointer to array (matrix) to receive the rotation matrix.

 OUTPUT:
 Updates the values of array pointed by M.

*/

 void 
 rot_mat_ZYX 
 (
   double **M,
   double const *s
 )
 {
   double cos_a = cos(s[0]);
   double sin_a = sin(s[0]);
   double cos_b = cos(s[1]);
   double sin_b = sin(s[1]);
   double cos_g = cos(s[2]);
   double sin_g = sin(s[2]);

   M[0][0] =  cos_a * cos_b;				//  cos(a)cos(b)
   M[0][1] =  cos_a * sin_b * sin_g - sin_a * cos_g;	//  cos(a)sin(b)sin(g)-sin(a)cos(g)
   M[0][2] =  cos_a * sin_b * cos_g + sin_a * sin_g;	//  cos(a)sin(b)cos(g)+sin(a)sin(g)
   M[1][0] =  sin_a * cos_b;				//  sin(a)cos(b)
   M[1][1] =  sin_a * sin_b * sin_g + cos_a * cos_g;	//  sin(a)sin(b)sin(g)+cos(a)cos(g)
   M[1][2] =  sin_a * sin_b * cos_g - cos_a * sin_g;	//  sin(a)sin(b)cos(g)-cos(a)sin(g)
   M[2][0] = -sin_b;					// -sin(b)
   M[2][1] =  cos_b * sin_g;				//  cos(b)sin(g)
   M[2][2] =  cos_b * cos_g;				//  cos(b)cos(g)
 }

/* rot_mat_ZYZ

 Calculates the rotation matrix based on sequential Z, Y and Z rotations. Since
 the polar angles of the structure's Z-axis are provided by the input file and
 the directional cosines will already be available, the uses the directional
 cosines directly. The spherical coordinates are also passed however, for the 
 case that Y-axis rotation is 0º or 180º (theta = 0º or 180º) as sin(theta) = 0.
 
 INPUT:
 u     - Pointer to array containing the directional cosines of the Z-axis;
 s     - Pointer to array containing the spherical coordinates of the Z-axis, 
         as theta,phi,radius;
 M     - Pointer to the GSL matrix to receive the rotation matrix;
 alpha - Rotation angle around the Z-axis, radians, [0,2 PI].

 OUTPUT:
 Updates the values of M;
*/

 void 
 rot_mat_ZYZ 
 (
   double **M,
   double const *u,
   double const *s,
   double const alpha
 )
 {
   double cos_a = cos(alpha);
   double sin_a = sin(alpha);

   if ( 1.0 - ABS(u[2]) < TOLERANCE )
   {
     double cos_p = cos(s[1]);
     double sin_p = sin(s[1]);

     M[0][0] =  cos_p * cos_a - sin_p * sin_a;				//  cos(p)cos(a)-sin(p)sin(a)
     M[0][1] = -cos_p * sin_a - sin_p * cos_a;				// -cos(p)sin(a)-sin(p)cos(a)
     M[0][2] =  0.0; 							//  0.0
     M[1][0] =  sin_p * cos_a + cos_p * sin_a;				//  sin(p)cos(a)+cos(p)sin(a)
     M[1][1] = -sin_p * sin_a + cos_p * cos_a;				// -sin(p)sin(a)+cos(p)cos(a)
     M[1][2] =  0.0;		 					//  0.0
     M[2][0] =  0.0;							//  0.0
     M[2][1] =  0.0;							//  0.0
     M[2][2] =  1.0;							//  1.0
   }
   else
   {
     double sin_t = sqrt(1.0 - u[2] * u[2]);
     M[0][0] =  u[0] * u[2] * cos_a / sin_t - u[1] * sin_a / sin_t;	//  cos(p)cos(t)cos(a)-sin(p)sin(a)
     M[0][1] = -u[0] * u[2] * sin_a / sin_t - u[1] * cos_a / sin_t;	// -cos(p)cos(t)sin(a)-sin(p)cos(a)
     M[0][2] =  u[0]; 							//  cos(p)sin(t)
     M[1][0] =  u[1] * u[2] * cos_a / sin_t + u[0] * sin_a / sin_t;	//  sin(p)cos(t)cos(a)+cos(p)sin(a)
     M[1][1] = -u[1] * u[2] * sin_a / sin_t + u[0] * cos_a / sin_t;	// -sin(p)cos(t)sin(a)+cos(p)cos(a)
     M[1][2] =  u[1]; 							//  sin(p)sin(t)
     M[2][0] = -sin_t * cos_a;						// -sin(t)cos(a)
     M[2][1] =  sin_t * sin_a;						//  sin(t)sin(a)
     M[2][2] =  u[2];							//  cos(t)
   }
 }


