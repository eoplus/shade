
/*******************************************************************************
 scattering_hg.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 Functions for Henyey-Greestein scattering.
 
 In the text below, equations Ax refers to equations in the scattering.c source
 file.

 ...

*******************************************************************************/

 #include <stdio.h>
 #include <string.h>
 #include <math.h>

 #include "config.h"
 #include "constants.h"
 #include "scattering.h"
 #include "scattering_hg.h"

 void
 scat_setup_hg
 ( struct scattering * scat,
   double const g
 )
 {
   scat->fun.qtl = &scat_qtl_hg;
   scat->fun.cdf = &scat_cdf_hg;
   scat->fun.pdf = &scat_pdf_hg;
   scat->fun.pf  = &scat_pf_hg;

//   scat->hg = (struct scattering_hg*) malloc(sizeof(struct scattering_hg));

   scat->g        = g;
   scat->hg.k_2g = 2.0 * scat->g;
   scat->hg.k_g2 = scat->g * scat->g;
   scat->fbb      = ((1.0 - g) / scat->hg.k_2g) * 
                    (((1.0 + g) / sqrt(1.0 + scat->hg.k_g2)) - 1.0);
   scat->pfn      = (1.0 / (4.0 * M_PI)) * (1.0 - scat->hg.k_g2);
   scat->truncate = 0;
   strncpy(scat->method, "analytical", STRMXLEN);
 }

 void
 scat_qtl_hg
 (
   struct scattering const * scat,
   double * s_scat,
   double const * prob 
 )
 {
   s_scat[0] = acos((1.0 + scat->hg.k_g2 - pow((1.0 - scat->hg.k_g2) / 
     (1.0 + scat->g - scat->hg.k_2g * prob[0]), 2.0)) / scat->hg.k_2g );
   s_scat[1] = K_2PI * prob[1];
 }

 void
 scat_cdf_hg
 (
   struct scattering const * scat,
   double const * s_scat,
   double * prob
 )
 {
   prob[0] = ((1.0 - scat->hg.k_g2) / scat->hg.k_2g) * 
     ( (1.0 / (1.0 - scat->g)) - 
     (1.0 / sqrt(1.0 + scat->hg.k_g2 - scat->hg.k_2g * cos(s_scat[0]))) );
   prob[1] = K_1_2PI * s_scat[1];
 }

 void
 scat_pdf_hg
 (
   struct scattering const * scat,
   double const * s_scat,
   double * dens
 )
 {
   scat_pf_hg(scat, s_scat, &dens[0]);
   dens[0] *= sin(s_scat[0]) * K_2PI;
   dens[1] = K_1_2PI;
 }

 void
 scat_pf_hg
 (
   struct scattering const * scat,
   double const * s_scat,
   double * pf
 )
 {
   pf[0] = scat->pfn / 
     pow(1.0 + scat->hg.k_g2 - scat->hg.k_2g * cos(s_scat[0]), 1.5);
 }

