
/*******************************************************************************
 val_rotation.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 This is a driver program to test if the calculation of the rotation matrices 
 are correctly implemented.

 Rotation matrices have the following properties:
 (1) The determinant is equal to -1 (improper) or 1 (proper);
 (2) The inverse is equal to the transpose;

 To test if the calculated rotations correspond with the desired directions,
 simple angle choices are made for which orientation is simple to determine. 
 Visual checks are performed in R.

 Compile instruction:
 gcc val_rotation.c ../../../rotation.c ../../../geometry.c ../../../memory.c \
   -lm -lgsl -O3 -o val_rotation.o

 ./val_rotation.o
*******************************************************************************/

 #include <stdio.h>
 #include <math.h>
 #include <sys/time.h>

 #include "../../../config.h"
 #include "../../../aux.h"
 #include "../../../constants.h"
 #include "../../../memory.h"
 #include "../../../rotation.h"
 #include "../../../geometry.h"

 void main ( void )
 {
   double **M_ori = calloc_2d(3, 3, "1");
   double **M_inv = calloc_2d(3, 3, "1");
   double **M_tsp = calloc_2d(3, 3, "1");
   double *u_ori  = calloc_1d(3, "1");
   double *u_rot  = calloc_1d(3, "1");
   double *v_ori  = calloc_1d(3, "1");
   double *v_rot  = calloc_1d(3, "1");
   double *s_ori  = calloc_1d(3, "1");
   double *s_rot  = calloc_1d(3, "1");
   double *o      = calloc_1d(3, "1");
 
   double det;
   int ci;

   gsl_rng *random;
   struct timeval tv;
   unsigned long int seed = tv.tv_sec + tv.tv_usec;
   gettimeofday (&tv, 0);
   random = gsl_rng_alloc (gsl_rng_ranlxs0);
   gsl_rng_set (random, seed);

  /* Validation for rot_mat_ZYX ***********************************************/
   int n_set01 = 5;
   int *res_set01 = (int *) calloc(n_set01, sizeof(int));

   printf("\n%s", ANSI_BOLD);
   printf("rot_mat_ZYX ===================================================\n%s",
     ANSI_RESET);

   s_ori[0] = 60.0 * RAD;	// alpha, roll
   s_ori[1] = 45.0 * RAD;	// beta, yaw
   s_ori[2] = 90.0 * RAD;	// gamma, pitch
   printf("\n");
   printf("  Rotation around Z-axis: %.2lfº (alpha)\n", s_ori[0] * DEG);
   printf("  Rotation around Y-axis: %.2lfº (beta)\n", s_ori[1] * DEG);
   printf("  Rotation around X-axis: %.2lfº (gamma)\n", s_ori[2] * DEG);
   printf("\n");

   rot_mat_ZYX(M_ori, s_ori);
   mat_fprintf(stdout, (double const **) M_ori, 3, 3, 1.0, "R[ZYX]", 1);
   printf("\n");
   mat_transpose(M_tsp, (double const **) M_ori, 3, 3);
   mat_inv(M_inv, (double const **) M_ori, 3, 3);

   mat_fprintf(stdout, (double const **) M_tsp, 3, 3, 1.0, "TRANSPOSE R[ZYX]", 1);
   printf("\n");
   mat_fprintf(stdout, (double const **) M_inv, 3, 3, 1.0, "INVERSE R[ZYX]", 1);
   printf("\n");

   ci = 0;
   det = mat_det( (double const **) M_ori, 3, 3);
   res_set01[ci] = (1.0 - det) < TOLERANCE;
   printf ("  Determinant: %.3e %s %s %s\n", det, res_set01[ci] ? ANSI_BLUE : 
     ANSI_RED, res_set01[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

   ci = 1;
   res_set01[ci] = mat_equal((double const **) M_inv, (double const **) M_tsp, 
     3, 3, TOLERANCE);
   printf ("  A^-1 == A^T: %s %s %s\n\n", res_set01[ci] ? ANSI_BLUE : 
     ANSI_RED, res_set01[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

   ci = 2;
   u_rot[0] = M_ori[0][2];
   u_rot[1] = M_ori[1][2];
   u_rot[2] = M_ori[2][2];
   cos_to_sph_unit(s_rot, u_rot, random);
   res_set01[ci] = NUM_EQU(M_PI_2, s_rot[0], TOLERANCE) && 
   NUM_EQU(M_PI + M_PI_2 + s_ori[0], s_rot[1], TOLERANCE);
   printf("  New Z-axis: %.2lfº %.2lfº\n", s_rot[0] * DEG, s_rot[1] * DEG);
   printf("  Z-axis rotation: %s %s %s\n\n", res_set01[ci] ? ANSI_BLUE : 
     ANSI_RED, res_set01[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

   ci = 3;
   u_rot[0] = M_ori[0][1];
   u_rot[1] = M_ori[1][1];
   u_rot[2] = M_ori[2][1];
   cos_to_sph_unit(s_rot, u_rot, random);
   res_set01[ci] = NUM_EQU(s_ori[0], s_rot[1], TOLERANCE) && 
     NUM_EQU(s_ori[1], s_rot[0], TOLERANCE);
   printf("  New Y-axis: %.2lfº %.2lfº\n", s_rot[0] * DEG, s_rot[1] * DEG);
   printf("  Y-axis rotation: %s %s %s\n\n", res_set01[ci] ? ANSI_BLUE : 
     ANSI_RED, res_set01[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

   ci = 4;
   u_rot[0] = M_ori[0][0];
   u_rot[1] = M_ori[1][0];
   u_rot[2] = M_ori[2][0];
   cos_to_sph_unit(s_rot, u_rot, random);
   res_set01[ci] = NUM_EQU(M_PI - s_ori[1], s_rot[0], TOLERANCE) && 
     NUM_EQU(s_ori[0], s_rot[1], TOLERANCE);
   printf("  New X-axis: %.2lfº %.2lfº\n", s_rot[0] * DEG, s_rot[1] * DEG);
   printf("  X-axis rotation: %s %s %s\n\n", res_set01[ci] ? ANSI_BLUE : 
     ANSI_RED, res_set01[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

  /* Validation for rot_mat_ZY0 ***********************************************/
   int n_set02 = 5;
   int *res_set02 = (int *) calloc(n_set02, sizeof(int));

   printf("\n%s", ANSI_BOLD);
   printf("rot_mat_ZY0 ===================================================\n%s",
     ANSI_RESET);

   s_ori[0] = 45.0 * RAD;	// Theta
   s_ori[1] = 60.0 * RAD;	// phi
   s_ori[2] = 1.0;		// radius
   printf("\n");
   printf("  Rotation around Z-axis: %.2lfº (phi)\n", s_ori[1] * DEG);
   printf("  Rotation around Y-axis: %.2lfº (theta)\n", s_ori[0] * DEG);
   printf("\n");
   sph_to_cos_unit(u_ori, s_ori);
   vec_fprintf(stdout, u_ori, 3, 1.0, "Cosine direction of new Z-axis: X,Y,Z", 1);
   printf("\n");

   rot_mat_ZY0(M_ori, u_ori);
   mat_fprintf(stdout, (double const **) M_ori, 3, 3, 1.0, "R[ZY0]", 1);
   printf("\n");
   mat_transpose(M_tsp, (double const **) M_ori, 3, 3);
   mat_inv(M_inv, (double const **) M_ori, 3, 3);

   mat_fprintf(stdout, (double const **) M_tsp, 3, 3, 1.0, "TRANSPOSE R[ZY0]", 1);
   printf("\n");
   mat_fprintf(stdout, (double const **) M_inv, 3, 3, 1.0, "INVERSE R[ZY0]", 1);
   printf("\n");

   ci = 0;
   det = mat_det((double const **) M_ori , 3, 3);
   res_set02[ci] = (1.0 - det) < TOLERANCE;
   printf ("  Determinant: %.3e %s %s %s\n", det, res_set02[ci] ? ANSI_BLUE : 
     ANSI_RED, res_set02[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

   ci = 1;
   res_set02[ci] = mat_equal((double const **) M_inv, (double const **) M_tsp, 
     3, 3, TOLERANCE);
   printf ("  A^-1 == A^T: %s %s %s\n\n", res_set02[ci] ? ANSI_BLUE : 
     ANSI_RED, res_set02[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

   ci = 2;
   u_rot[0] = M_ori[0][2];
   u_rot[1] = M_ori[1][2];
   u_rot[2] = M_ori[2][2];
   cos_to_sph_unit(s_rot, u_rot, random);
   res_set02[ci] = NUM_EQU(s_ori[1], s_rot[1], TOLERANCE) && 
     NUM_EQU(s_ori[0], s_rot[0], TOLERANCE);
   printf("  New Z-axis: %.2lfº %.2lfº\n", s_rot[0] * DEG, s_rot[1] * DEG);
   printf("  Z-axis rotation: %s %s %s\n\n", res_set02[ci] ? ANSI_BLUE : 
     ANSI_RED, res_set02[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

   ci = 3;
   u_rot[0] = M_ori[0][1];
   u_rot[1] = M_ori[1][1];
   u_rot[2] = M_ori[2][1];
   cos_to_sph_unit(s_rot, u_rot, random);
   res_set02[ci] = NUM_EQU(M_PI_2, s_rot[0], TOLERANCE) && 
   NUM_EQU(M_PI_2 + s_ori[1], s_rot[1], TOLERANCE);
   printf("  New Y-axis: %.2lfº %.2lfº\n", s_rot[0] * DEG, s_rot[1] * DEG);
   printf("  Y-axis rotation: %s %s %s\n\n", res_set02[ci] ? ANSI_BLUE : 
     ANSI_RED, res_set02[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

   ci = 4;
   u_rot[0] = M_ori[0][0];
   u_rot[1] = M_ori[1][0];
   u_rot[2] = M_ori[2][0];
   cos_to_sph_unit(s_rot, u_rot, random);
   res_set02[ci] = NUM_EQU(M_PI - s_ori[0], s_rot[0], TOLERANCE) && 
     NUM_EQU(s_ori[1], s_rot[1], TOLERANCE);
   printf("  New X-axis: %.2lfº %.2lfº\n", s_rot[0] * DEG, s_rot[1] * DEG);
   printf("  X-axis rotation: %s %s %s\n\n", res_set02[ci] ? ANSI_BLUE : 
     ANSI_RED, res_set02[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

  /* Validation for rot_mat_ZYZ ***********************************************/

   int n_set03 = 5;
   int *res_set03 = (int *) calloc(n_set03, sizeof(int));

   printf("\n%s", ANSI_BOLD);
   printf("rot_mat_ZYZ ===================================================\n%s",
     ANSI_RESET);

   s_ori[0] = 45.0 * RAD;	// Theta
   s_ori[1] = 60.0 * RAD;	// Phi
   s_ori[2] = 90.0 * RAD;	// Alpha
   printf("\n");
   printf("  Rotation around Z-axis (1st): %.2lfº (phi)\n", s_ori[1] * DEG);
   printf("  Rotation around Y-axis:       %.2lfº (theta)\n", s_ori[0] * DEG);
   printf("  Rotation around Z-axis (2nd): %.2lfº (alpha)\n", s_ori[2] * DEG);
   printf("\n");
   sph_to_cos_unit(u_ori, s_ori);
   vec_fprintf(stdout, u_ori, 3, 1.0, "Cosine direction of new Z-axis: X,Y,Z", 1);
   printf("\n");

   rot_mat_ZYZ(M_ori, u_ori, s_ori, s_ori[2]);
   mat_fprintf(stdout, (double const **) M_ori, 3, 3, 1.0, "R[ZYZ]", 1);
   printf("\n");
   mat_transpose(M_tsp, (double const **) M_ori, 3, 3);
   mat_inv(M_inv, (double const **) M_ori, 3, 3);

   mat_fprintf(stdout, (double const **) M_tsp, 3, 3, 1.0, "TRANSPOSE R[ZYZ]", 1);
   printf("\n");
   mat_fprintf(stdout, (double const **) M_inv, 3, 3, 1.0, "INVERSE R[ZYZ]", 1);
   printf("\n");

   ci = 0;
   det = mat_det((double const **) M_ori, 3, 3);
   res_set03[ci] = (1.0 - det) < TOLERANCE;
   printf ("  Determinant: %.3e %s %s %s\n", det, res_set03[ci] ? ANSI_BLUE : 
     ANSI_RED, res_set03[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

   ci = 1;
   res_set03[ci] = mat_equal((double const **) M_inv, (double const **) M_tsp,
     3, 3, TOLERANCE);
   printf ("  A^-1 == A^T: %s %s %s\n\n", res_set03[ci] ? ANSI_BLUE : 
     ANSI_RED, res_set03[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

   ci = 2;
   u_rot[0] = M_ori[0][2];
   u_rot[1] = M_ori[1][2];
   u_rot[2] = M_ori[2][2];
   cos_to_sph_unit(s_rot, u_rot, random);
   res_set03[ci] = NUM_EQU(s_ori[1], s_rot[1], TOLERANCE) && 
     NUM_EQU(s_ori[0], s_rot[0], TOLERANCE);
   printf("  New Z-axis: %.2lfº %.2lfº\n", s_rot[0] * DEG, s_rot[1] * DEG);
   printf("  Z-axis rotation: %s %s %s\n\n", res_set03[ci] ? ANSI_BLUE : 
     ANSI_RED, res_set03[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

   ci = 3;
   u_rot[0] = M_ori[0][1];
   u_rot[1] = M_ori[1][1];
   u_rot[2] = M_ori[2][1];
   cos_to_sph_unit(s_rot, u_rot, random);
   res_set03[ci] = NUM_EQU(M_PI_4, s_rot[0], TOLERANCE) && 
     NUM_EQU(M_PI + s_ori[1], s_rot[1], TOLERANCE);
   printf("  New Y-axis: %.2lfº %.2lfº\n", s_rot[0] * DEG, s_rot[1] * DEG);
   printf("  Y-axis rotation: %s %s %s\n\n", res_set03[ci] ? ANSI_BLUE : 
     ANSI_RED, res_set03[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

   ci = 4;
   u_rot[0] = M_ori[0][0];
   u_rot[1] = M_ori[1][0];
   u_rot[2] = M_ori[2][0];
   cos_to_sph_unit(s_rot, u_rot, random);
   res_set03[ci] = NUM_EQU(M_PI_2, s_rot[0], TOLERANCE) && 
     NUM_EQU(M_PI_2 + s_ori[1], s_rot[1], TOLERANCE);
   printf("  New X-axis: %.2lfº %.2lfº\n", s_rot[0] * DEG, s_rot[1] * DEG);
   printf("  X-axis rotation: %s %s %s\n\n", res_set03[ci] ? ANSI_BLUE : 
     ANSI_RED, res_set03[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

  /* Validation for rot_vec ***************************************************/
   int n_set04 = 3;
   int *res_set04 = (int *) calloc(n_set04, sizeof(int));

   printf("\n%s", ANSI_BOLD);
   printf("rot_vec =======================================================\n%s",
     ANSI_RESET);

   s_ori[0] = 45.0 * RAD;	// Theta
   s_ori[1] = 60.0 * RAD;	// Phi
   s_ori[2] = 90.0 * RAD;	// Alpha
   printf("\n");
   printf("  Rotation around Z-axis (1st): %.2lfº (phi)\n", s_ori[1] * DEG);
   printf("  Rotation around Y-axis:       %.2lfº (theta)\n", s_ori[0] * DEG);
   printf("  Rotation around Z-axis (2nd): %.2lfº (alpha)\n", s_ori[2] * DEG);
   printf("\n");
   sph_to_cos_unit(u_ori, s_ori);
   vec_fprintf(stdout, u_ori, 3, 1.0, "Cosine direction of new Z-axis: X,Y,Z", 1);
   printf("\n");

   rot_mat_ZYZ(M_ori, u_ori, s_ori, s_ori[2]);
   mat_fprintf(stdout, (double const **) M_ori, 3, 3, 1.0, "R[ZYZ]", 1);
   printf("\n");

   ci = 0;
   v_ori[0] = 0.0;
   v_ori[1] = 0.0;
   v_ori[2] = 1.0;
   rot_vec_unit (v_rot, v_ori, (double const **) M_ori);
   res_set04[ci] =
     NUM_EQU(v_rot[0], M_ori[0][2], TOLERANCE) && 
     NUM_EQU(v_rot[1], M_ori[1][2], TOLERANCE) && 
     NUM_EQU(v_rot[2], M_ori[2][2], TOLERANCE); 
   printf("  New Z-axis: %.2e %.2e %.2e\n", v_rot[0], v_rot[1], v_rot[2]);
   printf("  Z-axis rotation: %s %s %s\n\n", res_set04[ci] ? ANSI_BLUE : 
     ANSI_RED, res_set04[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

   ci = 1;
   v_ori[0] = 0.0;
   v_ori[1] = 1.0;
   v_ori[2] = 0.0;
   rot_vec_unit (v_rot, v_ori, (double const **) M_ori);
   res_set04[ci] =
     NUM_EQU(v_rot[0], M_ori[0][1], TOLERANCE) && 
     NUM_EQU(v_rot[1], M_ori[1][1], TOLERANCE) && 
     NUM_EQU(v_rot[2], M_ori[2][1], TOLERANCE); 
   printf("  New Y-axis: %.2e %.2e %.2e\n", v_rot[0], v_rot[1], v_rot[2]);
   printf("  Y-axis rotation: %s %s %s\n\n", res_set04[ci] ? ANSI_BLUE : 
     ANSI_RED, res_set04[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

   ci = 2;
   v_ori[0] = 1.0;
   v_ori[1] = 0.0;
   v_ori[2] = 0.0;
   rot_vec_unit (v_rot, v_ori, (double const **) M_ori);
   res_set04[ci] =
     NUM_EQU(v_rot[0], M_ori[0][0], TOLERANCE) && 
     NUM_EQU(v_rot[1], M_ori[1][0], TOLERANCE) && 
     NUM_EQU(v_rot[2], M_ori[2][0], TOLERANCE); 
   printf("  New X-axis: %.2e %.2e %.2e\n", v_rot[0], v_rot[1], v_rot[2]);
   printf("  X-axis rotation: %s %s %s\n\n", res_set04[ci] ? ANSI_BLUE : 
     ANSI_RED, res_set04[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

  /* Summary ******************************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("Summary =======================================================\n%s",
     ANSI_RESET);
   int res = 0.0;
   for (size_t i = 0; i < n_set01; i++) res += res_set01[i];
   printf("\n  rot_mat_ZYX: %s %s %s\n", res == n_set01 ? ANSI_BLUE : 
     ANSI_RED, res == n_set01 ? "PASSED!" : "FAILED!", ANSI_RESET);
   res = 0.0;
   for (size_t i = 0; i < n_set02; i++) res += res_set02[i];
   printf("  rot_mat_ZY0: %s %s %s\n", res == n_set02 ? ANSI_BLUE : 
     ANSI_RED, res == n_set02 ? "PASSED!" : "FAILED!", ANSI_RESET);
   res = 0.0;
   for (size_t i = 0; i < n_set03; i++) res += res_set03[i];
   printf("  rot_mat_ZYZ: %s %s %s\n", res == n_set03 ? ANSI_BLUE : 
     ANSI_RED, res == n_set03 ? "PASSED!" : "FAILED!", ANSI_RESET);
   res = 0.0;
   for (size_t i = 0; i < n_set04; i++) res += res_set04[i];
   printf("  rot_vec:     %s %s %s\n\n", res == n_set04 ? ANSI_BLUE : 
     ANSI_RED, res == n_set04 ? "PASSED!" : "FAILED!", ANSI_RESET);

   printf("\nRunning R script val_rotation.R:\n");
   system("Rscript val_rotation.R");
   printf("Done!\n");
 }


