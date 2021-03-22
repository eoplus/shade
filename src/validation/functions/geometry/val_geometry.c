
/*******************************************************************************
 val_geometry.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 This is a driver program to test if the operations in vectors and metrices are
 correctly performed.

 Direction for which results are easily verified were used.

 Compile instruction:
 gcc val_geometry.c ../../../geometry.c ../../../memory.c -lm -lgsl -O3 -o val_geometry.o

 ./val_geometry.o
*******************************************************************************/

 #include <stdio.h>
 #include <sys/time.h>

 #include "../../../config.h"
 #include "../../../constants.h"
 #include "../../../memory.h"
 #include "../../../geometry.h"

 void main ( void )
 {
 
   double *v = calloc_1d(3, "v in val_geometry.c");
   double *u = calloc_1d(3, "u in val_geometry.c");
   double *s = calloc_1d(3, "s in val_geometry.c");

   double **M = calloc_2d(3, 3, "M in val_geometry.c");
   double **N = calloc_2d(3, 3, "N in val_geometry.c");
   double **I = calloc_2d(3, 3, "I in val_geometry.c");

   // Initiate random number generator:
   struct timeval tv;
   unsigned long int seed;
   gettimeofday (&tv, 0);
   seed   = tv.tv_sec + tv.tv_usec;
   gsl_rng *random = gsl_rng_alloc (gsl_rng_ranlxs0);
   gsl_rng_set (random, seed);

  /* Validation for sph and cos ***********************************************/
   printf("\n%s", ANSI_BOLD);
   printf("sph and cos ===================================================\n%s",
     ANSI_RESET);

   printf("\n  Unit radius: \n");
   s[0] = 60.0 * RAD;
   s[1] = 45.0 * RAD;
   s[2] = 1.0;
   sph_to_cos_unit(u, s);
   cos_to_sph_unit(v, u, random);
   s[2] *= RAD;
   vec_fprintf(stdout, s, 3, DEG, "Spherical", 2);
   vec_fprintf(stdout, u, 3, 1.0, "Cartesian", 2);
   v[2] *= RAD;
   vec_fprintf(stdout, v, 3, DEG, "Back to spherical", 2);

   printf("\n  Radius = 5 \n");
   s[2] = 5.0;
   sph_to_cos(u, s);
   cos_to_sph(v, u, random);
   s[2] *= RAD;
   vec_fprintf(stdout, s, 3, DEG, "Spherical", 2);
   vec_fprintf(stdout, u, 3, 1.0, "Cartesian", 2);
   v[2] *= RAD;
   vec_fprintf(stdout, v, 3, DEG, "Back to spherical", 2);

  /* Validation for mat_inv ***************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("mat_inv =======================================================\n%s",
     ANSI_RESET);

   M[0][0] = 2.0;
   M[1][1] = 2.0;
   M[2][2] = 2.0;
   mat_fprintf(stdout, (double const **) M, 3, 3, 1.0, "Original matrix", 1);   
   mat_inv (N, (double const **) M, 3, 3);
   mat_fprintf(stdout, (double const **) N, 3, 3, 1.0, "Inverse matrix", 1);   

  /* Validation for mat_det ***************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("mat_det =======================================================\n%s",
     ANSI_RESET);

   M[0][0] = 2.0;
   M[1][1] = 2.0;
   M[2][2] = 2.0;
   mat_fprintf(stdout, (double const **) M, 3, 3, 1.0, "Original matrix", 1);
   printf("\n  Determinant: %.3e\n", mat_det((double const **) M, 3, 3));
  
  /* Validation for mat_equal *************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("mat_equal =====================================================\n%s",
     ANSI_RESET);

   N[0][0] = 2.0;
   N[1][1] = 2.0;
   N[2][2] = 2.0;
   mat_fprintf(stdout, (double const **) M, 3, 3, 1.0, "Matrix 1", 1);
   mat_fprintf(stdout, (double const **) N, 3, 3, 1.0, "Matrix 2", 1);
   printf("  Equal: %d\n", mat_equal((double const **) M, (double const **) N, 3, 
     3, TOLERANCE));

  /* Validation for mat_vec_prod **********************************************/
   printf("\n%s", ANSI_BOLD);
   printf("mat_vec_prod ==================================================\n%s",
     ANSI_RESET);

   I[0][0] = 1.0;
   I[1][1] = 1.0;
   I[2][2] = 1.0;
   mat_vec_prod(v, u, (double const **) I, 3, 3);
   mat_fprintf(stdout, (double const **) I, 3, 3, 1.0, "Identity matrix", 1);
   printf("\n");
   vec_fprintf(stdout, u, 3, 1.0, "Original vector", 1);
   printf("\n");
   vec_fprintf(stdout, v, 3, 1.0, "Product", 1);
 }

