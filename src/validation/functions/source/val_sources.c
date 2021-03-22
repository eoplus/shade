
/*******************************************************************************
 val_sources.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 This is a driver program to test if the sources are properly implemented.

 Compile instruction:
 gcc val_sources.c ../../../rotation.c ../../../geometry.c ../../../memory.c ../../../sources.c -lm -lgsl -O3 -o val_sources.o

 ./val_sources.o
*******************************************************************************/

 #include <stdio.h>
 #include <string.h>
 #include <sys/time.h>
 #include <gsl/gsl_rng.h>

 #include "../../../config.h"
 #include "../../../constants.h"
 #include "../../../sources.h"

 void main ( void )
 {
   FILE   *fpo;
   double step;
   int N = 1000000;
   int i;

   double src_fov;
   double src_s[3];
   double src_ref_o[3];
   double src_rel_o[3] = {0.0};
   double src_stks[1];
   char src_tp[STRMXLEN];

   struct source *src = src_alloc( );

   double prob[2];
   double prob_1[2];
   double prob_2[2];
   double s_scat[3];
   double stks[1];

   // Initiate random number generator:
   struct timeval tv;
   unsigned long int seed = tv.tv_sec + tv.tv_usec;
   gettimeofday (&tv, 0);
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
   src_ref_o[0] = 0.0;
   src_ref_o[1] = 0.0;
   src_ref_o[2] = 0.0;
   src_stks[0] = 1.0;
   strncpy(src_tp, "sir", STRMXLEN);

   src_setup(src, src_fov, src_s, src_ref_o, src_rel_o, src_stks, src_tp);
   src_fprintf(stdout, src, 1);
   printf("\n");

   prob[0] = 0.5;
   prob[1] = 0.5;
   src_qtl(src, prob, s_scat);
   printf("QTL:\n");
   printf("Prob: %.2lf %.2lf\n", prob[0], prob[1]);
   printf("Scat: %.2lf %.2lf\n\n", s_scat[0] * DEG, s_scat[1] * DEG);

   s_scat[0] = 60.0 * RAD;
   s_scat[1] = 180.0 * RAD;
   src_cdf(src, prob, s_scat);
   printf("CDF:\n");
   printf("Scat: %.2lf %.2lf\n", s_scat[0] * DEG, s_scat[1] * DEG);
   printf("CDF: %.2lf %.2lf\n\n", prob[0], prob[1]);

   src_pdf(src, prob, s_scat);
   printf("PDF:\n");
   printf("Scat: %.2lf %.2lf\n", s_scat[0] * DEG, s_scat[1] * DEG);
   printf("PDF: %.2lf %.2lf\n\n", prob[0], prob[1]);

   src_ef(src, stks, s_scat);
   printf("EF:\n");
   printf("Scat: %.2lf %.2lf\n", s_scat[0] * DEG, s_scat[1] * DEG);
   printf("STKS: %.2lf\n\n", stks[0]);

   fpo = fopen ("plots/ed0_eu0.txt", "w");
   fprintf(fpo, "Angle\tEF\tPDF\tCDF\n");
   src_fov = 180.0 * RAD;
   src_setup(src, src_fov, src_s, src_ref_o, src_rel_o, src_stks, src_tp);
   step = M_PI_2 / 720.0;
   s_scat[0] = 0.0;
   for (i = 0; i < 721; i++)
   {
     src_ef(src, stks, s_scat); 
     src_pdf(src, prob_1, s_scat);
     src_cdf(src, prob_2, s_scat);

     fprintf(fpo, "%.6E\t%.6E\t%.6E\t%.6E\n", 
       s_scat[0] * DEG,
       stks[0], 
       prob_1[0],
       prob_2[0]);
     s_scat[0] += step;
   }
   fclose(fpo);

   fpo = fopen ("plots/ed0_eu0_smp.txt", "w");
   fprintf(fpo, "psi\n");
   for (i = 0; i < N; i++)
   {
     prob[0] = gsl_rng_uniform(random);
     prob[1] = gsl_rng_uniform(random);
     src_qtl(src, prob, s_scat);
     fprintf(fpo, "%06.2f\n", s_scat[0] * DEG, s_scat[1] * DEG);
   }
   fclose(fpo);

  /* Val pir: *****************************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("PNT_PIR =======================================================\n%s",
     ANSI_RESET);

   src_fov = 180.0 * RAD;
   src_s[0] = 0.0 * RAD;
   src_s[1] = 0.0 * RAD;
   src_s[2] = 1.0;
   src_ref_o[0] = 0.0;
   src_ref_o[1] = 0.0;
   src_ref_o[2] = 0.0;
   src_stks[0] = 1.0;
   strncpy(src_tp, "pir", STRMXLEN);

   src_setup(src, src_fov, src_s, src_ref_o, src_rel_o, src_stks, src_tp);
   src_fprintf(stdout, src, 1);
   printf("\n");

   prob[0] = 0.5;
   prob[1] = 0.5;
   src_qtl(src, prob, s_scat);
   printf("QTL:\n");
   printf("Prob: %.2lf %.2lf\n", prob[0], prob[1]);
   printf("Scat: %.2lf %.2lf\n\n", s_scat[0] * DEG, s_scat[1] * DEG);

   s_scat[0] = 60.0 * RAD;
   s_scat[1] = 180.0 * RAD;
   src_cdf(src, prob, s_scat);
   printf("CDF:\n");
   printf("Scat: %.2lf %.2lf\n", s_scat[0] * DEG, s_scat[1] * DEG);
   printf("CDF: %.2lf %.2lf\n\n", prob[0], prob[1]);

   src_pdf(src, prob, s_scat);
   printf("PDF:\n");
   printf("Scat: %.2lf %.2lf\n", s_scat[0] * DEG, s_scat[1] * DEG);
   printf("PDF: %.2lf %.2lf\n\n", prob[0], prob[1]);

   src_ef(src, stks, s_scat);
   printf("EF:\n");
   printf("Scat: %.2lf %.2lf\n", s_scat[0] * DEG, s_scat[1] * DEG);
   printf("STKS: %.2lf\n\n", stks[0]);

   fpo = fopen ("plots/ed_eu.txt", "w");
   fprintf(fpo, "Angle\tEF\tPDF\tCDF\n");
   src_fov = 180.0 * RAD;
   src_setup(src, src_fov, src_s, src_ref_o, src_rel_o, src_stks, src_tp);
   step = M_PI_2 / 720.0;
   s_scat[0] = 0.0;
   for (i = 0; i < 721; i++)
   {
     src_ef(src, stks, s_scat); 
     src_pdf(src, prob_1, s_scat);
     src_cdf(src, prob_2, s_scat);

     fprintf(fpo, "%.6E\t%.6E\t%.6E\t%.6E\n", 
       s_scat[0] * DEG,
       stks[0], 
       prob_1[0],
       prob_2[0]);
     s_scat[0] += step;
   }
   fclose(fpo);

   fpo = fopen ("plots/ed_eu_smp.txt", "w");
   fprintf(fpo, "psi\n");
   for (i = 0; i < N; i++)
   {
     prob[0] = gsl_rng_uniform(random);
     prob[1] = gsl_rng_uniform(random);
     src_qtl(src, prob, s_scat);
     fprintf(fpo, "%06.2f\n", s_scat[0] * DEG, s_scat[1] * DEG);
   }
   fclose(fpo);
  /* Val rad: *****************************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("PNT_RAD =======================================================\n%s",
     ANSI_RESET);

   src_fov = 60.0 * RAD;
   src_s[0] = 0.0 * RAD;
   src_s[1] = 0.0 * RAD;
   src_s[2] = 1.0;
   src_ref_o[0] = 0.0;
   src_ref_o[1] = 0.0;
   src_ref_o[2] = 0.0;
   src_stks[0] = 1.0;
   strncpy(src_tp, "rad", STRMXLEN);

   src_setup(src, src_fov, src_s, src_ref_o, src_rel_o, src_stks, src_tp);
   src_fprintf(stdout, src, 1);
   printf("\n");

   prob[0] = 0.5;
   prob[1] = 0.5;
   src_qtl(src, prob, s_scat);
   printf("QTL:\n");
   printf("Prob: %.2lf %.2lf\n", prob[0], prob[1]);
   printf("Scat: %.2lf %.2lf\n\n", s_scat[0] * DEG, s_scat[1] * DEG);

   s_scat[0] = 30.0 * RAD;
   s_scat[1] = 180.0 * RAD;
   src_cdf(src, prob, s_scat);
   printf("CDF:\n");
   printf("Scat: %.2lf %.2lf\n", s_scat[0] * DEG, s_scat[1] * DEG);
   printf("CDF: %.2lf %.2lf\n\n", prob[0], prob[1]);

   src_pdf(src, prob, s_scat);
   printf("PDF:\n");
   printf("Scat: %.2lf %.2lf\n", s_scat[0] * DEG, s_scat[1] * DEG);
   printf("PDF: %.2lf %.2lf\n\n", prob[0], prob[1]);

   src_ef(src, stks, s_scat);
   printf("EF:\n");
   printf("Scat: %.2lf %.2lf\n", s_scat[0] * DEG, s_scat[1] * DEG);
   printf("STKS: %.2lf\n\n", stks[0]);

   fpo = fopen ("plots/lu_10.txt", "w");
   fprintf(fpo, "Angle\tEF\tPDF\tCDF\n");
   src_fov = 10.0 * RAD;
   src_setup(src, src_fov, src_s, src_ref_o, src_rel_o, src_stks, src_tp);
   step = (5.0 * RAD) / 720.0;
   s_scat[0] = 0.0;
   for (i = 0; i < 721; i++)
   {
     src_ef(src, stks, s_scat); 
     src_pdf(src, prob_1, s_scat);
     src_cdf(src, prob_2, s_scat);

     fprintf(fpo, "%.6E\t%.6E\t%.6E\t%.6E\n", 
       s_scat[0] * DEG,
       stks[0], 
       prob_1[0],
       prob_2[0]);
     s_scat[0] += step;
   }
   fclose(fpo);

   fpo = fopen ("plots/lu_10_smp.txt", "w");
   fprintf(fpo, "psi\n");
   for (i = 0; i < N; i++)
   {
     prob[0] = gsl_rng_uniform(random);
     prob[1] = gsl_rng_uniform(random);
     src_qtl(src, prob, s_scat);
     fprintf(fpo, "%06.2f\n", s_scat[0] * DEG, s_scat[1] * DEG);
   }
   fclose(fpo);

   printf("\nRunning R script val_source.R:\n");
   system("Rscript val_sources.R");
   printf("Done!\n");
 }
