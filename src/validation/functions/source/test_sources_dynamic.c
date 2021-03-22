/*******************************************************************************
 test_sources_dynamic.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 This is a driver program to test if the speed when using using dynamic
 dispatching

 Compile instruction:
 gcc test_sources_dynamic.c ../../../src/rotation.c ../../../src/geometry.c ../../../src/memory.c ../../../src/sources.c -lm -lgsl -O3 -o test_sources_dynamic.o

 ./test_sources_dynamic.o
*******************************************************************************/

 #include <stdio.h>
 #include <string.h>
 #include <time.h>
 #include <sys/time.h>
 #include <gsl/gsl_rng.h>

 #include "../../../src/config.h"
 #include "../../../src/constants.h"
 #include "../../../src/sources.h"

 void main ( void )
 {
   double src_fov;
   double src_s[3];
   double src_o[3];
   double src_stks[1];
   char src_tp[STRMXLEN];

   struct source *src = src_alloc( );

   double prob[2];
   double s_scat[3];
   double stks[1];

   time_t start, end;

   // Initiate random number generator:
   struct timeval tv;
   unsigned long int seed;
   gettimeofday (&tv, 0);
   seed   = tv.tv_sec + tv.tv_usec;
   gsl_rng *random = gsl_rng_alloc (gsl_rng_ranlxs0);
   gsl_rng_set (random, seed);

  /* Val sir: *****************************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("PNT_SIR =======================================================\n%s",
     ANSI_RESET);

   src_fov = 180.0 * RAD;
   src_s[0] = 0.0 * RAD;
   src_s[1] = 0.0 * RAD;
   src_s[2] = 1.0;
   src_o[0] = 0.0;
   src_o[1] = 0.0;
   src_o[2] = 0.0;
   src_stks[0] = 1.0;

   FILE *fi = fopen("input_source.txt", "r");
   if(fi == NULL)
   {
     printf("\nERROR: Failed to open input file input_source.txt\n\n");
     exit(-1);
   }
   fscanf(fi, "%s", src_tp);

   src_setup(src, src_fov, src_s, src_o, src_stks, src_tp);
   src_printf(src);
   printf("\n");

   time(&start);
   for (size_t i = 0; i < 100000000; i++)
   {
     prob[0] = gsl_rng_uniform(random);
     prob[1] = gsl_rng_uniform(random);
     src_qtl(src, prob, s_scat);
   }
   time(&end);
   double total = (double) (end - start);
   printf("Dynamic execution time: %.2e\n\n", total);

   time(&start);
   for (size_t i = 0; i < 100000000; i++)
   {
     prob[0] = gsl_rng_uniform(random);
     prob[1] = gsl_rng_uniform(random);
     src_qtl_pnt_pir(src, prob, s_scat);
   }
   time(&end);
   total = (double) (end - start);
   printf("Direct execution time: %.2e\n\n", total);
 }

