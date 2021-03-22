
/*******************************************************************************
 val_reflectance.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 This is a driver program to test if the BRDF functions are properly
 implemented.

 Compile instruction:
 gcc val_reflectance.c ../../../src/reflectance.c ../../../src/memory.c -lm -lgsl -O3 -o val_reflectance.o

 ./val_reflectance.o
*******************************************************************************/

 #include <stdio.h>
 #include <string.h>
 #include <sys/time.h>

 #include "../../../src/config.h"
 #include "../../../src/constants.h"
 #include "../../../src/reflectance.h"

 #include <gsl/gsl_spline.h>
 #include <gsl/gsl_rng.h>

 void main ( void )
 {
   struct reflectance *refl = refl_alloc();
   char   refl_tp[STRMXLEN];
   double k;
   double bhr[1] = {1.0};
   int    nbhr = 1;
   double step1, step2;

   FILE   *fpo;
   int N = 1000000;

   double s_scat_i[2];
   double s_scat_r[2];
   double dens[2];
   double prob[2];
   double brdf[1];

   // Initiate random number generator:
   struct timeval tv;
   unsigned long int seed;
   gettimeofday (&tv, 0);
   seed   = tv.tv_sec + tv.tv_usec;
   gsl_rng *random = gsl_rng_alloc (gsl_rng_ranlxs0);
   gsl_rng_set (random, seed);

  /* Lambertian: **************************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("Lambertian ====================================================\n%s",
     ANSI_RESET);

   strncpy(refl_tp, "lambert", STRMXLEN);
   refl_setup(refl, refl_tp, nbhr, bhr, k);
   refl_printf(refl);
   printf("\n");

   fpo = fopen ("plots/lambert_cdf_i_00_000.txt", "w");
   fprintf(fpo, "Psi_r\tPhi_r\tPDF_Psi_r\tCDF_Psi_r\tPDF_Phi_r\tCDF_Phi_r\n");
   s_scat_i[0] = cos(0.0 * RAD);
   s_scat_i[1] = 0.0;
   s_scat_r[0] = cos(0.0);
   s_scat_r[1] = 0.0;
   step1 = M_PI_2 / 360.0;
   step2 = K_2PI  / 360.0;
   for (size_t i = 0; i < 361; i++)
   {
     refl_pdf(refl, s_scat_i, s_scat_r, dens);
     refl_cdf(refl, s_scat_i, s_scat_r, prob);

     fprintf(fpo, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n", 
       acos(s_scat_r[0]) * DEG,
       s_scat_r[1] * DEG,
       dens[0], 
       prob[0],
       dens[1],
       prob[1]);
     s_scat_r[0] = cos(acos(s_scat_r[0]) + step1);
     s_scat_r[1] += step2;
   }
   fclose(fpo);

   fpo = fopen ("plots/lambert_brdf_i_00_000_smp.txt", "w");
   fprintf(fpo, "Psi_r\tPhi_r\n");
   s_scat_i[0] = cos(0.0 * RAD);
   s_scat_i[1] = 0.0;
   for (size_t i = 0; i < N; i++)
   {
     prob[0] = gsl_rng_uniform (random);
     prob[1] = gsl_rng_uniform (random);
     refl_qtl(refl, s_scat_i, s_scat_r, prob);
     fprintf(fpo, "%.6E\t%.6E\n", acos(s_scat_r[0]) * DEG, s_scat_r[1] * DEG);
   }
   fclose(fpo);

   fpo = fopen ("plots/lambert_cdf_i_60_000.txt", "w");
   fprintf(fpo, "Psi_r\tPhi_r\tPDF_Psi_r\tCDF_Psi_r\tPDF_Phi_r\tCDF_Phi_r\n");
   s_scat_i[0] = cos(60.0 * RAD);
   s_scat_i[1] = 0.0;
   s_scat_r[0] = cos(0.0);
   s_scat_r[1] = 0.0;
   step1 = M_PI_2 / 360.0;
   step2 = K_2PI  / 360.0;
   for (size_t i = 0; i < 361; i++)
   {
     refl_pdf(refl, s_scat_i, s_scat_r, dens);
     refl_cdf(refl, s_scat_i, s_scat_r, prob);

     fprintf(fpo, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n", 
       acos(s_scat_r[0]) * DEG,
       s_scat_r[1] * DEG,
       dens[0], 
       prob[0],
       dens[1],
       prob[1]);
     s_scat_r[0] = cos(acos(s_scat_r[0]) + step1);
     s_scat_r[1] += step2;
   }
   fclose(fpo);

   fpo = fopen ("plots/lambert_brdf_i_60_000_smp.txt", "w");
   fprintf(fpo, "Psi_r\tPhi_r\n");
   s_scat_i[0] = cos(60.0 * RAD);
   s_scat_i[1] = 0.0;
   for (size_t i = 0; i < N; i++)
   {
     prob[0] = gsl_rng_uniform (random);
     prob[1] = gsl_rng_uniform (random);
     refl_qtl(refl, s_scat_i, s_scat_r, prob);
     fprintf(fpo, "%.6E\t%.6E\n", acos(s_scat_r[0]) * DEG, s_scat_r[1] * DEG);
   }
   fclose(fpo);

  /* Minnaert: **************************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("Minnaert ======================================================\n%s",
     ANSI_RESET);

   k = 1.0;
   strncpy(refl_tp, "minnaert", STRMXLEN);
   refl_setup(refl, refl_tp, nbhr, bhr, k);
   refl_printf(refl);
   printf("\n");

   fpo = fopen ("plots/minnaert_cdf_i_00_000.txt", "w");
   fprintf(fpo, "Psi_r\tPhi_r\tPDF_Psi_r\tCDF_Psi_r\tPDF_Phi_r\tCDF_Phi_r\n");
   s_scat_i[0] = cos(0.0 * RAD);
   s_scat_i[1] = 0.0;
   s_scat_r[0] = cos(0.0);
   s_scat_r[1] = 0.0;
   step1 = M_PI_2 / 360.0;
   step2 = K_2PI  / 360.0;
   for (size_t i = 0; i < 361; i++)
   {
     refl_pdf(refl, s_scat_i, s_scat_r, dens);
     refl_cdf(refl, s_scat_i, s_scat_r, prob);

     fprintf(fpo, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n", 
       acos(s_scat_r[0]) * DEG,
       s_scat_r[1] * DEG,
       dens[0], 
       prob[0],
       dens[1],
       prob[1]);
     s_scat_r[0] = cos(acos(s_scat_r[0]) + step1);
     s_scat_r[1] += step2;
   }
   fclose(fpo);

   fpo = fopen ("plots/minnaert_brdf_i_00_000_smp.txt", "w");
   fprintf(fpo, "Psi_r\tPhi_r\n");
   s_scat_i[0] = cos(0.0 * RAD);
   s_scat_i[1] = 0.0;
   for (size_t i = 0; i < N; i++)
   {
     prob[0] = gsl_rng_uniform (random);
     prob[1] = gsl_rng_uniform (random);
     refl_qtl(refl, s_scat_i, s_scat_r, prob);
     fprintf(fpo, "%.6E\t%.6E\n", acos(s_scat_r[0]) * DEG, s_scat_r[1] * DEG);
   }
   fclose(fpo);

   fpo = fopen ("plots/minnaert_cdf_i_60_000.txt", "w");
   fprintf(fpo, "Psi_r\tPhi_r\tPDF_Psi_r\tCDF_Psi_r\tPDF_Phi_r\tCDF_Phi_r\n");
   s_scat_i[0] = cos(60.0 * RAD);
   s_scat_i[1] = 0.0;
   s_scat_r[0] = cos(0.0);
   s_scat_r[1] = 0.0;
   step1 = M_PI_2 / 360.0;
   step2 = K_2PI  / 360.0;
   for (size_t i = 0; i < 361; i++)
   {
     refl_pdf(refl, s_scat_i, s_scat_r, dens);
     refl_cdf(refl, s_scat_i, s_scat_r, prob);

     fprintf(fpo, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n", 
       acos(s_scat_r[0]) * DEG,
       s_scat_r[1] * DEG,
       dens[0], 
       prob[0],
       dens[1],
       prob[1]);
     s_scat_r[0] = cos(acos(s_scat_r[0]) + step1);
     s_scat_r[1] += step2;
   }
   fclose(fpo);

   fpo = fopen ("plots/minnaert_brdf_i_60_000_smp.txt", "w");
   fprintf(fpo, "Psi_r\tPhi_r\n");
   s_scat_i[0] = cos(60.0 * RAD);
   s_scat_i[1] = 0.0;
   for (size_t i = 0; i < N; i++)
   {
     prob[0] = gsl_rng_uniform (random);
     prob[1] = gsl_rng_uniform (random);
     refl_qtl(refl, s_scat_i, s_scat_r, prob);
     fprintf(fpo, "%.6E\t%.6E\n", acos(s_scat_r[0]) * DEG, s_scat_r[1] * DEG);
   }
   fclose(fpo);

   printf("\nRunning R script val_reflectance.R:\n");
   system("Rscript val_reflectance.R");
   printf("Done!\n");
 }
