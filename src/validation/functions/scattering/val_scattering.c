
/*******************************************************************************
 val_scattering.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 This is a driver program to test if the scattering functions are properly
 implemented.

 Compile instruction:
 gcc val_scattering.c ../../../scattering.c ../../../scattering_iso.c ../../../scattering_rlg.c ../../../scattering_hg.c ../../../scattering_ff.c ../../../scattering_truncate.c -lm -lgsl -O3 -o val_scattering.o

 ./val_scattering.o
*******************************************************************************/

 #include <stdio.h>
 #include <string.h>
 #include <sys/time.h>

 #include "../../../config.h"
 #include "../../../constants.h"
 #include "../../../scattering.h"

 #include <gsl/gsl_spline.h>
 #include <gsl/gsl_rng.h>

 void main ( void )
 {
   struct scattering *scat = scat_alloc();
   char   scat_tp[STRMXLEN];
   double depolr;
   double g;
   double fbb;
   char   method[STRMXLEN];
   int    truncate;
   double step;

   FILE   *fpo;
   int N = 1000000;

   double s_scat[2];
   double dens[2];
   double prob[2];
   double pf[1];

   // Initiate random number generator:
   struct timeval tv;
   unsigned long int seed;
   gettimeofday (&tv, 0);
   seed   = tv.tv_sec + tv.tv_usec;
   gsl_rng *random = gsl_rng_alloc (gsl_rng_ranlxs0);
   gsl_rng_set (random, seed);

  /* Isotropic: ***************************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("Isotropic =====================================================\n%s",
     ANSI_RESET);

   strncpy(scat_tp, "isotropic", STRMXLEN);

   scat_setup(scat, scat_tp, depolr, g, fbb, method, truncate);
   scat_fprintf(stdout, scat, 1);
   printf("\n");

   printf("Average cosine: %.4e\n\n", scat_av_mu(scat));

   fpo = fopen ("plots/isotropic.txt", "w");
   fprintf(fpo, "Angle\tPF\tPDF\tCDF\n");
   step = M_PI / 720.0;
   s_scat[0] = 0.0;
   for (size_t i = 0; i < 721; i++)
   {
     scat_pf(scat, s_scat, pf); 
     scat_pdf(scat, s_scat, dens);
     scat_cdf(scat, s_scat, prob);

     fprintf(fpo, "%.6E\t%.6E\t%.6E\t%.6E\n", 
       s_scat[0] * DEG,
       pf[0], 
       dens[0],
       prob[0]);
     s_scat[0] += step;
   }
   fclose(fpo);

   fpo = fopen ("plots/isotropic_smp.txt", "w");
   fprintf(fpo, "psi\n");
   for (size_t i = 0; i < N; i++)
   {
     prob[0] = gsl_rng_uniform(random);
     prob[1] = gsl_rng_uniform(random);
     scat_qtl(scat, s_scat, prob);
     fprintf(fpo, "%06.2f\n", s_scat[0] * DEG, s_scat[1] * DEG);
   }
   fclose(fpo);

  /* Rayleigh air: ************************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("Rayleigh (air) ================================================\n%s",
     ANSI_RESET);

   strncpy(scat_tp, "rayleigh", STRMXLEN);
   depolr = 1.0;

   scat_setup(scat, scat_tp, depolr, g, fbb, method, truncate);
   scat_fprintf(stdout, scat, 1);
   printf("\n");

   printf("Average cosine: %.4e\n\n", scat_av_mu(scat));

   fpo = fopen ("plots/rayleigh_air.txt", "w");
   fprintf(fpo, "Angle\tPF\tPDF\tCDF\n");
   step = M_PI / 720.0;
   s_scat[0] = 0.0;
   for (size_t i = 0; i < 721; i++)
   {
     scat_pf(scat, s_scat, pf); 
     scat_pdf(scat, s_scat, dens);
     scat_cdf(scat, s_scat, prob);

     fprintf(fpo, "%.6E\t%.6E\t%.6E\t%.6E\n", 
       s_scat[0] * DEG,
       pf[0], 
       dens[0],
       prob[0]);
     s_scat[0] += step;
   }
   fclose(fpo);

   fpo = fopen ("plots/rayleigh_air_smp.txt", "w");
   fprintf(fpo, "psi\n");
   for (size_t i = 0; i < N; i++)
   {
     prob[0] = gsl_rng_uniform(random);
     prob[1] = gsl_rng_uniform(random);
     scat_qtl(scat, s_scat, prob);
     fprintf(fpo, "%06.2f\n", s_scat[0] * DEG, s_scat[1] * DEG);
   }
   fclose(fpo);

  /* Rayleigh water: **********************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("Rayleigh (water) ==============================================\n%s",
     ANSI_RESET);

   strncpy(scat_tp, "rayleigh", STRMXLEN);
   depolr = WATER_DEPOLR;

   scat_setup(scat, scat_tp, depolr, g, fbb, method, truncate);
   scat_fprintf(stdout, scat, 1);
   printf("\n");

   printf("Average cosine: %.4e\n\n", scat_av_mu(scat));

   fpo = fopen ("plots/rayleigh_water.txt", "w");
   fprintf(fpo, "Angle\tPF\tPDF\tCDF\n");
   step = M_PI / 720.0;
   s_scat[0] = 0.0;
   for (size_t i = 0; i < 721; i++)
   {
     scat_pf(scat, s_scat, pf); 
     scat_pdf(scat, s_scat, dens);
     scat_cdf(scat, s_scat, prob);

     fprintf(fpo, "%.6E\t%.6E\t%.6E\t%.6E\n", 
       s_scat[0] * DEG,
       pf[0], 
       dens[0],
       prob[0]);
     s_scat[0] += step;
   }
   fclose(fpo);

   fpo = fopen ("plots/rayleigh_water_smp.txt", "w");
   fprintf(fpo, "psi\n");
   for (size_t i = 0; i < N; i++)
   {
     prob[0] = gsl_rng_uniform(random);
     prob[1] = gsl_rng_uniform(random);
     scat_qtl(scat, s_scat, prob);
     fprintf(fpo, "%06.2f\n", s_scat[0] * DEG, s_scat[1] * DEG);
   }
   fclose(fpo);

  /* Henyey-Greestein: ********************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("Henyey-Greestein ==============================================\n%s",
     ANSI_RESET);

   strncpy(scat_tp, "hg", STRMXLEN);
   g = 0.924;

   scat_setup(scat, scat_tp, depolr, g, fbb, method, truncate);
   scat_fprintf(stdout, scat, 1);
   printf("\n");

   printf("Average cosine: %.4e\n\n", scat_av_mu(scat));

   fpo = fopen ("plots/hg.txt", "w");
   fprintf(fpo, "Angle\tPF\tPDF\tCDF\n");
   step = M_PI / 720.0;
   s_scat[0] = 0.0;
   for (size_t i = 0; i < 721; i++)
   {
     scat_pf(scat, s_scat, pf); 
     scat_pdf(scat, s_scat, dens);
     scat_cdf(scat, s_scat, prob);

     fprintf(fpo, "%.6E\t%.6E\t%.6E\t%.6E\n", 
       s_scat[0] * DEG,
       pf[0], 
       dens[0],
       prob[0]);
     s_scat[0] += step;
   }
   fclose(fpo);

   fpo = fopen ("plots/hg_smp.txt", "w");
   fprintf(fpo, "psi\n");
   for (size_t i = 0; i < N; i++)
   {
     prob[0] = gsl_rng_uniform(random);
     prob[1] = gsl_rng_uniform(random);
     scat_qtl(scat, s_scat, prob);
     fprintf(fpo, "%06.2f\n", s_scat[0] * DEG, s_scat[1] * DEG);
   }
   fclose(fpo);

  /* Fournier-Forand: *********************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("Fournier-Forand ===============================================\n%s",
     ANSI_RESET);

   strncpy(scat_tp, "ff", STRMXLEN);
   fbb = 0.0183;
   strcpy(method, "itp");
   truncate = 0;

   scat_setup(scat, scat_tp, depolr, g, fbb, method, truncate);
   scat_fprintf(stdout, scat, 1);
   printf("\n");

   printf("Average cosine: %.4e\n\n", scat_av_mu(scat));

   fpo = fopen ("plots/ff_full.txt", "w");
   fprintf(fpo, "Angle\tPF\tPDF\tCDF\n");
   step = M_PI / 720.0;
   s_scat[0] = 0.0;
   for (size_t i = 0; i < 720; i++)
   {
     scat_pf(scat, s_scat, pf);
     scat_pdf(scat, s_scat, dens);
     scat_cdf(scat, s_scat, prob);

     fprintf(fpo, "%.6E\t%.6E\t%.6E\t%.6E\n", 
       s_scat[0] * DEG,
       pf[0], 
       dens[0],
       prob[0]);
     s_scat[0] += step;
   }
   fclose(fpo);

   fpo = fopen ("plots/ff_full_smp.txt", "w");
   fprintf(fpo, "psi\n");
   for (size_t i = 0; i < N; i++)
   {
     prob[0] = gsl_rng_uniform(random);
     prob[1] = gsl_rng_uniform(random);
     scat_qtl(scat, s_scat, prob);
     fprintf(fpo, "%06.2f\n", s_scat[0] * DEG, s_scat[1] * DEG);
   }
   fclose(fpo);

  /* Fournier-Forand (truncated): *********************************************/
   printf("\n%s", ANSI_BOLD);
   printf("Fournier-Forand (truncated) ===================================\n%s",
     ANSI_RESET);

   strncpy(scat_tp, "ff", STRMXLEN);
   fbb = 0.0183;
   strcpy(method, "itp");
   truncate = 1;

   scat_setup(scat, scat_tp, depolr, g, fbb, method, truncate);
   scat_fprintf(stdout, scat, 1);
   printf("\n");

   printf("Average cosine: %.4e\n\n", scat_av_mu(scat));

   fpo = fopen ("plots/ff_trunc.txt", "w");
   fprintf(fpo, "Angle\tPF\tPDF\tCDF\n");
   step = M_PI / 720.0;
   s_scat[0] = 0.0;
   for (size_t i = 0; i < 720; i++)
   {
     scat_pf(scat, s_scat, pf);
     scat_pdf(scat, s_scat, dens);
     scat_cdf(scat, s_scat, prob);

     fprintf(fpo, "%.6E\t%.6E\t%.6E\t%.6E\n", 
       s_scat[0] * DEG,
       pf[0], 
       dens[0],
       prob[0]);
     s_scat[0] += step;
   }
   fclose(fpo);

   fpo = fopen ("plots/ff_trunc_smp.txt", "w");
   fprintf(fpo, "psi\n");
   for (size_t i = 0; i < N; i++)
   {
     prob[0] = gsl_rng_uniform(random);
     prob[1] = gsl_rng_uniform(random);
     scat_qtl(scat, s_scat, prob);
     fprintf(fpo, "%06.2f\n", s_scat[0] * DEG, s_scat[1] * DEG);
   }
   fclose(fpo);


   printf("\nRunning R script val_scattering.R:\n");
   system("Rscript val_scattering.R");
   printf("Done!\n");
 }
