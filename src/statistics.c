
/*******************************************************************************
 statistics.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 Version: 1.6
 Date: 2021-03-25
 License: GPL-3.0

*******************************************************************************/

 #include <stdio.h>	// Functions:  fprintf, printf,
 #include <string.h>	// Functions:  strcat, strcpy, sprintf
 #include <stdlib.h>	// Functions:  exit
 #include <math.h>      // Functions:  pow, sqrt

 #include "config.h"		// STRMXLEN
 #include "constants.h"
 #include "accumulators.h"	// Structures: str_acm

/*
 * At the moment, only calculates the SE for the total spatial integral. It is
 * simple to extend to do the same for the spatially resolved integral
 *
 */

 void
 calc_se
 (
   struct accumulator_bmc *accm_se,
   struct accumulator_bmc *accm_mn,
   double *btt_bhr,
   char   *sufx,
   char   *ofbn
 )
 {
   int cw, cb, cz, cj;
   char bf[10];
   char ifn[STRMXLEN], pre[STRMXLEN];
   double val;
   FILE *fpi;

   for (cj = 0; cj < SUBDV; cj++)
   {
     sprintf(pre, "_S%02d", cj);
     for (cb = 0; cb < accm_se->nbr; cb++)
     {
       sprintf(bf, "_B%05.3f", btt_bhr[cb]);
       strncpy(ifn, ofbn, STRMXLEN - 1);
       strcat(ifn, pre);
       strcat(ifn, sufx);
       strcat(ifn, bf);
       strcat(ifn, ".txt");
 
       fpi = fopen(ifn, "r");
       if(fpi == NULL)
       {
         printf("Failed to open intermediate result file %s\n", ifn);
         exit(-1);
       }

       fscanf(fpi, "%*[^\n]\n");
       for (cz = 0; cz < accm_se->ns; cz++)
       {
         fscanf(fpi, "%*s\t");
         for (cw = 0; cw < accm_se->nw0; cw++)
         {
           fscanf(fpi, "%lf", &val);
           accm_se->vector[cz][cw][cb] += 
             pow(val - accm_mn->vector[cz][cw][cb], 2.0) / (SUBDV - 1);
         }
       }
       fclose(fpi);
       remove(ifn);
     }
   }

   for (cz = 0; cz < accm_se->ns; cz++)
   {
     for (cw = 0; cw < accm_se->nw0; cw++)
     {
       for (cb = 0; cb < accm_se->nbr; cb++)
       {
         accm_se->vector[cz][cw][cb] = sqrt(accm_se->vector[cz][cw][cb] / SUBDV);
       }
     }
   }
 }


