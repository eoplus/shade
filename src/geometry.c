
/*******************************************************************************
 geometry.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 Version: 1.6
 Date: 2021-03-25
 License: GPL-3.0

 Provides basic function for vectors and matrices. Small functions that are used
 many times in the code are declared as static inline and defined in the header
 file geometry.h.

 A bounded Euclidean vector is typically represented by its initial point A and 
 its terminal point B. A free vector has the initial point A at the origin of 
 the reference frame. Other properties of vectors are its magnitude (or 
 Euclidean norm, length), represented by the scalar m:

 m = ||B-A|| = sqrt( sum( (B[i]-A[i])^2 ) ),                                 (1)

 and its direction, represented by the unit vector u:
 
 u = (B-A) / ||B-A||,                                                        (2)
 
 with the unit vector u being a free vector (initial point at the origin) with
 unit Euclidean norm. The adjectives "row" and "column" before "vector" 
 represents a vector in matrix notation. A row vector is a matrix with 1 row and
 a column vector is a matrix with 1 column.

 The geometry of the simulation is defined in three-dimensions. A right-hand 
 coordinate system is used, with the Z-axis pointing downwards. This is the 
 convention used in the reference work of Mobley (1994). The axis vectors are 
 defined by the unit vectors:

      x  y  z
 X = {1, 0, 0};
 Y = {0, 1, 0};
 Z = {0, 0, 1};

 Spherical coordinates give a position in terms of:

 theta  - Angle from the Z-axis;
 phi    - Angle from the X-axis; 
 radius - Distance from the origin.

 The relation between the spherical and Cartesian coordinates is:

 x = radius * sin(theta) * cos(phi)
 y = radius * sin(theta) * cos(phi)
 z = radius * cos(theta),                                                    (3)

 with the reverse being:

 theta  = acos(z / ||{x,y,z}||),
 phi    = atan2(y, x),
 radius = ||{x,y,z}||.                                                       (4)

 Any position in the coordinate system can be reached by scaling the Z-axis and 
 performing two consecutive rotations, starting with phi (counterclockwise 
 rotation around the Z-axis) and then theta (counterclockwise rotation around 
 the Y-axis). The rotation.c source file has further information on rotations 
 and rotation matrices. Matrices are used in this code mainly for rotation of 
 directions.

 No specific data type is defined for vectors and matrices. When GSL functions 
 are needed, the calling function performs a conversion from arrays to the GSL 
 vector and matrix data types as necessary. 

*******************************************************************************/

 #include <stdlib.h>
 #include <string.h>
 #include <math.h>
 #include <gsl/gsl_vector.h>
 #include <gsl/gsl_matrix.h>
 #include <gsl/gsl_permutation.h>
 #include <gsl/gsl_linalg.h>
 #include <gsl/gsl_rng.h>

 #include "config.h"
 #include "constants.h"
 #include "aux.h"
 #include "geometry.h"

/* Matrix inversion: ***********************************************************

 mat_inv

 Calculate the inverse of a matrix M.

 INPUT:
 N  - Pointer to an array of pointers (matrix) receiving the results
 M  - Pointer to an array of pointers (matrix) for which the inverse is to be 
      calculated;
 nr - Number of rows of the matrix;
 nc - Number of columns of the matrix;

 OUTPUT:
 Updates the values of matrix pointed by N with the inverse of the matrix 
 pointed by M.

*******************************************************************************/

 void
 mat_inv
 (
   double ** N, 
   double const ** M,
   int const nr,
   int const nc
 )
 {
   /* Copy matrices to arrays and setup GSL matrix views from array */
   int n = nc * nr;
   double * Na = (double*) malloc(n * sizeof(double));
   double * Ma = (double*) malloc(n * sizeof(double));
   for (int cr = 0; cr < nr; cr++)
   {
     for (int cc = 0; cc < nc; cc++)
     { 
       Na[cr * nc + cc] = N[cr][cc];
       Ma[cr * nc + cc] = M[cr][cc];
     }
   }

   gsl_matrix_view Nv = gsl_matrix_view_array(Na, nr, nc);
   gsl_matrix_view Mv = gsl_matrix_view_array(Ma, nr, nc);

   int s;
   gsl_matrix * LU = gsl_matrix_alloc( nr, nc );
   gsl_matrix_memcpy( LU, &Mv.matrix ) ;
   gsl_permutation * p  = gsl_permutation_alloc( nc );
   gsl_linalg_LU_decomp( LU, p, &s );    
   gsl_linalg_LU_invert( LU, p, &Nv.matrix );

   for (int cr = 0; cr < nr; cr++)
   {
     for (int cc = 0; cc < nc; cc++)
     { 
       N[cr][cc] = gsl_matrix_get(&Nv.matrix, cr, cc);
     }
   }

   gsl_matrix_free( LU );
   gsl_permutation_free( p );
   free(Na);
   free(Ma);
 }

/* Matrix determinant: *********************************************************

 mat_det

 Calculates the determinant of a matrix.

 INPUT:
 M  - Pointer to an array of pointers (matrix) for which the determinant is to 
      be calculated;
 nr - Number of rows of the matrix;
 nc - Number of columns of the matrix;

 OUTPUT:
 The determinant of matrix pointed by M.

*******************************************************************************/

 double
 mat_det 
 (
   double const ** M,
   int const nr,
   int const nc
 )
 {
   /* Copy matrices to arrays and setup GSL matrix views from array */
   int n = nc * nr;
   double * Ma = (double*) malloc(n * sizeof(double));
   for (int cr = 0; cr < nr; cr++)
   {
     for (int cc = 0; cc < nc; cc++)
     { 
       Ma[cr * nc + cc] = M[cr][cc];
     }
   }
   gsl_matrix_view Mv = gsl_matrix_view_array(Ma, nr, nc);

   int s;
   gsl_matrix * LU = gsl_matrix_alloc( nr, nc );
   gsl_matrix_memcpy( LU, &Mv.matrix ) ;
   gsl_permutation * p  = gsl_permutation_alloc( nc );
   gsl_linalg_LU_decomp( LU, p, &s );
   
   double det = gsl_linalg_LU_det( LU, s );
   gsl_matrix_free( LU );
   gsl_permutation_free( p );
   free(Ma);

   return det;
 }

/* Matrix transpose: ***********************************************************

 mat_transpose

 Transpose a square matrix M.

 INPUTS:
 N - Pointer to matrix pointer (not allocated);
 M - Matrix pointer (to be transposed);

 OUTPUT:
 None. Updates the values of matrix N with the transpose of M.

*******************************************************************************/

 void 
 mat_transpose
 (
   double ** N,
   double const ** M,
   int const nr,
   int const nc
 )
 {
   for (int cr = 0; cr < nr; cr++)
   {
     for (int cc = 0; cc < nc; cc++)
     {
       N[cc][cr] = M[cr][cc];
     }
   }
 }

/* Matrix (near) equality: ***************************************************** 

 mat_equal

 Function to compare the (approximate) equality of two matrices. It is similar
 to the GSL function gsl_matrix_equal, but provides a tolerance to compensate
 for rounding errors.

 INPUT:
 N   - Pointer to an array of pointers (matrix);
 M   - Pointer to an array of pointers (matrix);
 TOL - Tolerance for approximate equality;
 nr  - Number of rows of the matrix;
 nc  - Number of columns of the matrix;

 OUTPUT:
 0 if not (approximate) equal and 1 otherwise.

*******************************************************************************/

 int
 mat_equal
 (
   double const ** M, 
   double const ** N,
   int const nr,
   int const nc, 
   double const TOL
 )
 {
   for (int cr = 0; cr < nr; cr++)
   {
     for (int cc = 0; cc < nc; cc++)
     {
       if ( ABS(M[cr][cc] - N[cr][cc]) > TOL ) 
       {
         return 0;
       }
     }
   }
   return 1;
 }


/* Print vector: ***************************************************************

 vec_fprintf

 Formatted print of vector to file.

 INPUT:
 odv    - Output device. Should be a file stream (e.g. 'stdout' to print to
          default);
 v      - Pointer to an array;
 n      - Length of the array;
 scl    - Arbitrary multiplier; For example, to print angles in degrees, use
          DEG;
 title  - Pointer to a title string or NULL.
 indent - Number of indentations at the header level (each indent = 2 spaces)
          Constant int;
 OUTPUT:
 None. Prints vector. 

*******************************************************************************/

 void
 vec_fprintf
 (
   FILE * odv,
   double const * v,
   int const n,
   double const scl, 
   char const * title,
   int const indent
 )
 {
   // Set indentation:
   char pre_0[STRMXLEN] = "";
   char pre_1[STRMXLEN] = "";
   for (int i = 0; i < indent; i++)
     strcat(pre_0, "  ");
   strncpy(pre_1, pre_0, STRMXLEN);
   strcat(pre_1, "  ");

   if ( title )
   {
     fprintf(odv, "%s%s:\n", pre_0, title);
   } else {
     fprintf(odv, "%sVector:\n", pre_0);
   }
   fprintf(odv, "%sDimension: %d\n", pre_1, n);
   fprintf(odv, "%sData: ", pre_1);
   for (int i = 0; i < n; i++)
   {
    if ( i && !(i % 5) )
      fprintf(odv, "\n%s      ", pre_1);
     fprintf(odv, "% .3e ", v[i] * scl );
   }
   fprintf(odv, "\n");
 }

/* Print matrix: ***************************************************************
 
 mat_fprintf

 Formatted print of matrix to file.

 INPUT:
 odv    - Output device. Should be a file stream (e.g. 'stdout' to print to
          default);
 M      - Pointer to an array of pointers (matrix);
 nr     - Number of rows of M;
 nc     - Number of columns of M;
 scl    - Arbitrary multiplier; For example, to print angles in degrees, use 
          DEG;
 title  - Pointer to a title string or NULL.
 indent - Number of indentations at the header level (each indent = 2 spaces)
          Constant int;
 OUTPUT:
 None. Prints matrix.

*******************************************************************************/

 void
 mat_fprintf
 (
   FILE * odv,
   double const ** M,
   int const nr,
   int const nc, 
   double const scl, 
   char const * title,
   int const indent
 )
 {
   // Set indentation:
   char pre_0[STRMXLEN] = "";
   char pre_1[STRMXLEN] = "";
   for (int i = 0; i < indent; i++)
     strcat(pre_0, "  ");
   strncpy(pre_1, pre_0, STRMXLEN);
   strcat(pre_1, "  ");

   if ( title )
   {
     fprintf(odv, "%s%s:\n", pre_0, title);
   } else {
     fprintf(odv, "%sMatrix:\n", pre_0);
   }
   fprintf(odv, "%sDimension: %d, %d\n", pre_1, nr, nc);
   fprintf(odv, "%sData: ", pre_1);

   for (int i = 0; i < nr; i++)
   {
     if (i > 0) fprintf(odv, "%s      ", pre_1);
     for (int j = 0; j < nc; j++)
     {
       fprintf(odv, "% .3e ", M[i][j] * scl );
     }
     fprintf(odv, "\n");
   }
 }

