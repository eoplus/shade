
/*******************************************************************************
 scattering_iso.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 Functions for isotropic scattering.
 
 In the text below, equations Ax refers to equations in the scattering.c source
 file.

 Isotropic scattering results in equal scattering per solid angle in all 
 directions, therefore the phase function is constant and integrating Eq. A7 
 over all directions in the unit sphere results in:

 1 = 2 PI Integrate{dpsi B(psi,phi) sin(psi)} ->
 1 = 2 PI B Integral{dpsi sin(psi)} ->
 1 = 4 PI B ->
 B = 1 / (4 PI)                                                             (B1)

 By substituting B into Eq. A6 and A7 we have:

 PDF(psi) = 2 PI B(psi,phi) sin(psi) ->
          = 2 PI B sin(psi)
          = 0.5 sin(psi),                                                   (B2)

 CDF(psi) = Integral(0 to psi){dpsi PDF(psi)} ->
          = 0.5 Integral(0 to psi){dpsi sin(psi)} ->
          = 0.5 (1 - cos(psi)).                                             (B3)

 The quantile function solves Eq. B3 for psi as:

 CDF(psi) = 0.5 (1 - cos(psi)) ->
 2 CDF(psi) = 1 - cos(psi) ->
 psi = acos(1 - 2 CDF(psi))                                                 (B4)

 All angles are in the reference frame of the ray's direction and need to be 
 rotated to the system's reference frame. This is not performed by the 
 scattering functions. Random sampling from the PDF is made by providing a
 positive random uniform deviate as CDF(psi) (termed prob, for probability).

*******************************************************************************/

 #include <stdio.h>
 #include <string.h>
 #include <math.h>

 #include "config.h"
 #include "constants.h"
 #include "scattering.h"
 #include "scattering_iso.h"

 #define ISOTROPIC_PF  (1.0 / (4.0 * M_PI))

 void
 scat_setup_iso
 ( struct scattering * scat )
 {
   scat->fun.qtl = &scat_qtl_iso;
   scat->fun.cdf = &scat_cdf_iso;
   scat->fun.pdf = &scat_pdf_iso;
   scat->fun.pf  = &scat_pf_iso;

   scat->pfn = 1.0;
   scat->fbb = 0.5;
   scat->g   = 0.0;
   scat->truncate = 0;
   strncpy(scat->method, "analytical", STRMXLEN);
 }

 void
 scat_qtl_iso
 (
   struct scattering const * scat,
   double * s_scat,
   double const * prob 
 )
 {
   s_scat[0] = acos(1.0 - 2.0 * prob[0]);
   s_scat[1] = K_2PI * prob[1];
 }

 void
 scat_cdf_iso
 (
   struct scattering const * scat,
   double const * s_scat,
   double * prob
 )
 {
   prob[0] = 0.5 * (1.0 - cos(s_scat[0]));
   prob[1] = K_1_2PI * s_scat[1];
 }

 void
 scat_pdf_iso
 (
   struct scattering const * scat,
   double const * s_scat,
   double * dens
 ) 
 { 
   dens[0] = 0.5 * sin(s_scat[0]);
   dens[1] = K_1_2PI;
 }

 void
 scat_pf_iso
 (
   struct scattering const * scat,
   double const * s_scat,
   double * pf
 )
 { pf[0] = ISOTROPIC_PF; }


