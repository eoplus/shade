
/*******************************************************************************
 val_ray.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 This is a driver program to test if the functions to alloc, free, initiate and
 move a ray work as expected.

 Compile instruction:
 gcc val_ray.c ../../../ray.c ../../../geometry.c ../../../memory.c ../../../rotation.c -lm -lgsl -O3 -o val_ray.o

 ./val_ray.o
*******************************************************************************/

 #include <stdio.h>
 #include <math.h>
 #include <gsl/gsl_blas.h>
 #include <gsl/gsl_rng.h>
 #include <sys/time.h>

 #include "../../../constants.h"
 #include "../../../memory.h"
 #include "../../../geometry.h"
 #include "../../../ray.h"

 void main (void)
 {
   double l = 10.0;

   double *u = calloc_1d(3, "u in val_ray.c");
   double *s = calloc_1d(3, "s in val_ray.c");
   double *src_pos = calloc_1d(3, "src_pos in val_ray.c");

   size_t iop_nw0 = 1; 
   size_t btt_nbr = 1;
   struct light_ray* ray = ray_alloc (1, iop_nw0, btt_nbr);

   // Initiate random number generator:
   struct timeval tv;
   unsigned long int seed = tv.tv_sec + tv.tv_usec;
   gettimeofday (&tv, 0);
   gsl_rng *random = gsl_rng_alloc (gsl_rng_ranlxs0);
   gsl_rng_set (random, seed);

  /* Validation for ray_ini ***************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("ray_ini =======================================================\n%s",
     ANSI_RESET);
   printf("  Ray starts at the origin, with spherical directions: \n"
     "  60.0º, 45.0º, 1 (Polar, Azimuth, Radius)\n\n");

   s[0] = 60.0 * RAD;
   s[1] = 45.0 * RAD;
   s[2] = 1.0;
   ray_ini(src_pos, s, ray);
   ray_printf(ray, 1);

  /* Validation for ray_mov ***************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("ray_mov =======================================================\n%s",
     ANSI_RESET);
   printf("  Ray moved by a distance of 1 meter. Terminal point should be:\n"
     "  6.12, 6.12, 5.00 (m)\n\n");

   ray_mov(l, ray);
   ray_printf(ray, 1);

  /* Validation for ray_scat **************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("ray_scat ======================================================\n%s",
     ANSI_RESET);
   printf("  Ray scattered by an angle of 180.0º. Final spherical directions\n"
     "  should be: 120.0º, 225.0º\n\n");

   s[0] = 180.0 * RAD;
   s[1] = 0.0 * RAD;
   s[2] = 1.0;
   ray_scat(s, ray, random);
   ray_printf(ray, 1);

  /* Validation for ray_free **************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("ray_free ======================================================\n%s",
     ANSI_RESET);
   printf("  Releasing memory of ray should also set pointer to NULL\n\n");
   ray_free(&ray);
   printf("  Is pointer NULL: %d\n\n", ray == NULL);

 }
