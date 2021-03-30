
/*******************************************************************************
 scattering_rlg.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 Version: 1.6
 Date: 2021-03-25
 License: GPL-3.0

 The rlgleigh phase function describes molecular scattering as proportional to 
 1 + f cos(psi)^2 (Zhang et al., 2019), where the value of f is given by 
 (1 - d) / (1 + d) and d is the depolarization ratio of the medium. Following
 from the equations on the 'scattering.c' source file, we can calculate the
 normalization factor necessary for the phase function to integrate to 1 (Eq. A7
 integrated over the full range of angles):

 1 = 2 PI Integral{dpsi B(psi) sin(psi)} ->
 1 = 2 PI Integral{dpsi K (1 + f cos(psi)^2) sin(psi)} ->
 1 = K 4 PI (3 + f) / 3 ->
 K = 3 / (4 PI (3 + f))                                                     (B1)
 
 therefore, B(psi) is:

 B(psi) = K (1 + f cos(psi)^2) ->
        = 3 (1 + f cos(psi)^2) / (4 PI (3 + f))                             (B2)

 Substituting B(psi) in Eq. A6 and A7, the probability density function (PDF)
 and cumulative distribution function (CDF) are then:
 
 PDF(psi) = 2 PI B(psi) sin(psi) ->
          = 2 PI sin(psi) 3 (1 + f cos(psi)^2) / (4 PI (3 + f)) ->
          = (1 + f cos(psi)^2) sin(psi) 3 / (6 + 2 f),                      (B3)

 CDF(psi) = Integrate(0 to psi){dpsi PDF(psi)} ->
          = Integrate(0 to psi){dpsi (1 + f cos(psi)^2) sin(psi) 3 / (6 + 2 f)}
          = (3 / (6 + 2 f))
              Integrate(0 to psi){dpsi sin(psi) + f sin(psi) cos(psi)^2}->
          = -(f / (6 + 2 f)) cos(psi)^3 - (3 / (6 + 2 f)) cos(psi) + 0.5    (B4)

 The quantile function solves Eq. B4 for psi:

 CDF(psi) = -(f / (6 + 2 f)) cos(psi)^3 - (3 / (6 + 2 f)) cos(psi) + 0.5 ->
 psi      = acos(delta / (3 cbrt(2) f) - 3 cbrt(2) / delta),                (B5)

 delta    = cbrt( 54 f^3 (0.5 - CDF(psi)) + 162 f^2 (0.5 - CDF(psi)) + 
              sqrt(2916 f^3 + (54 f^3 (0.5 - CDF(psi)) + 
              162 f^2 (0.5 - CDF(psi)))^2)),                                (B6)

 where cbrt is the cubit root.
 
 To improve efficiency of the calculations, those equations can be simplified
 for given f values. The depolarization ratio of air is 0, such that f = 1. 
 The depolarization ratio of pure water (Farinato & Rowel, 1976) and pure 
 (artificial) seawater (Zhang et al., 2019) is 0.039, such that f = 0.9249278.
 There is some evidence that the depolarization ratio of water has a small
 dependency on salinity and some uncertain dependency on wavelength (Zhang et 
 al., 2019). To accomodate variable depolarization ratios, as may be of 
 interest in special simulations or be included in future version as a model
 for water, the general Rayleigh equations are provided below. However, since 
 the value of 0.039 is in general recommended for water this will be used, but
 no specific function is written and the water scattering is called with the 
 generic function for full precision.

 REFERENCES:

 Farinato, R. S; Rowell, R. L. 1976. New values of the light scattering
   depolarization and anisotropy of water. The Journal of Chemical Physics 65,
   2, 593-595.
   DOI: http://doi.org/10.1063/1.433115

 Zhang, X.; Stramski, D.; Reynolds, R. A.; Blocker, E. R. 2019. Light 
   scattering by pure water and seawater: the depolarization ratio and its
   variation with salinity. Applied Optics 58, 4, 991-1004. 
   DOI: https://doi.org/10.1364/AO.58.000991 

*******************************************************************************/

 #include <stdio.h>
 #include <string.h>
 #include <stdlib.h>
 #include <math.h>

 #include "config.h"
 #include "aux.h"
 #include "constants.h"
 #include "scattering.h"
 #include "scattering_rlg.h"

 #define RAYLEIGH_A_NORM (3.0 / (16.0 * M_PI))

/* Setup Rayleigh scattering: **************************************************

 scat_setup_rlg

 Setup Rayleigh scattering functions.

 INPUT:
 scat   - Scattering parameters
          Pointer to scattering struct;
 depolr - Linear depolarization ratio
          Range: (0,1], unitless
          Double.

 OUTPUT:
 None. Updates the values pointed by scat.

*******************************************************************************/

 void
 scat_setup_rlg
 (
   struct scattering * scat,
   double const depolr
 )
 {
//   scat->rlg = (struct scattering_rlg*) malloc (sizeof(struct scattering_rlg));
   scat->rlg.depolr = depolr;
   scat->rlg.f      = (1.0 - depolr) / (1.0 + depolr);

   scat->rlg.f2 = scat->rlg.f * scat->rlg.f;
   scat->rlg.f3 = scat->rlg.f2 * scat->rlg.f;
   scat->rlg.c0 = 3.0 / (4.0 * M_PI * (3.0 + scat->rlg.f));
   scat->rlg.c1 = 3.0 / (2.0 * (3.0 + scat->rlg.f));
   scat->rlg.c2 = scat->rlg.f / (2.0 * (3.0 + scat->rlg.f));
   scat->rlg.c3 = 54.0 * scat->rlg.f3 + 162.0 * scat->rlg.f2;
   scat->rlg.c4 = 2916.0 * scat->rlg.f3;
   scat->rlg.c5 = 1.0 / (3.0 * cbrt(2.0) * scat->rlg.f);
   scat->rlg.c6 = 3.0 * cbrt(2.0);

   if ( ABS(1.0 - depolr) < TOLERANCE )
   {
     // depolr = 1.0 (Air)
     scat->fun.qtl = &scat_qtl_rlg_air;
     scat->fun.cdf = &scat_cdf_rlg_air;
     scat->fun.pdf = &scat_pdf_rlg_air;
     scat->fun.pf  = &scat_pf_rlg_air;
   }
   else
   {
     scat->fun.qtl = &scat_qtl_rlg;
     scat->fun.cdf = &scat_cdf_rlg;
     scat->fun.pdf = &scat_pdf_rlg;
     scat->fun.pf  = &scat_pf_rlg;
   }

   scat->pfn = RAYLEIGH_A_NORM;
   scat->fbb = 0.5;
   scat->g   = 0.0;
   scat->truncate = 0;
   strncpy(scat->method, "analytical", STRMXLEN);
 }

/* Rayleigh scattering functions: **********************************************

 scat_qtl_rlg, scat_cdf_rlg, scat_pdf_rlg, scat_ef_rlg

 Functions to retrieve the quantile, CDF, PDF and phase function for Rayleigh
 scattering with arbitrary depolarization ratio.

 INPUT:
 scat   - Scattering parameters
          Pointer to constant scattering struct; 
 s_scat - Spherical direction
          [0]Theta, [1]Phi, [2]Radius
          Range: [0,PI] radians, [0,2PI] radians, 1.0 unitless
          Pointer to constant double or pointer to double;
 prob   - Cumulative probability associated with a polar and a azimuth angle
          [0]Theta, [1]Phi
          Range: [0,1] unitless
          Pointer to constant double or pointer to double;
 dens   - Probability density associated with a polar and a azimuth angle
          [0]Theta, [1]Phi
          Range: [0,Inf) unitless
          Pointer to double
 pf     - Phase function at a particular angle
          Range: [0,Inf), sr
          Pointer to double;

 OUTPUT:
 None. Updates values pointed by s_scat, prob or pf.

*******************************************************************************/

 void
 scat_qtl_rlg
 (
   struct scattering const * scat,
   double * s_scat,
   double const * prob
 )
 {
   double const b0 = 0.5 - prob[0];
   double const b1 = scat->rlg.c3 * b0;
   double const b2 = scat->rlg.c4 + b1 * b1;
   double const d0 = cbrt(b1 + sqrt(b2));
   s_scat[0] = acos(scat->rlg.c5 * d0 - scat->rlg.c6 / d0);
   s_scat[1] = K_2PI * prob[1];
 }

 void
 scat_cdf_rlg
 (
   struct scattering const * scat,
   double const * s_scat,
   double * prob
 )
 {
   double const mu = cos(s_scat[0]);
   prob[0] = -scat->rlg.c2 * mu * mu * mu - scat->rlg.c1 * mu + 0.5;
   prob[1] = K_1_2PI * s_scat[1];
 }

 void
 scat_pdf_rlg
 (
   struct scattering const * scat,
   double const * s_scat,
   double * dens
 ) 
 {
   double const mu = cos(s_scat[0]);
   dens[0] = (1.0 + scat->rlg.f * mu * mu) * sin(s_scat[0]) * scat->rlg.c1;
   dens[1] = K_1_2PI;
 }

 void
 scat_pf_rlg
 (
   struct scattering const * scat,
   double const * s_scat,
   double * pf
 )
 {
   double const mu = cos(s_scat[0]);
   pf[0] = (1.0 + scat->rlg.f * mu * mu) * scat->rlg.c0;
 }

/* Rayleigh for standard air: **************************************************

 scat_qtl_rlg_air, scat_cdf_rlg_air, scat_pdf_rlg_air, scat_ef_rlg_air

 Functions to retrieve the quantile, CDF, PDF and phase function for Rayleigh
 scattering for air (depolarization ratio = 1).

 Since f = 1, Eq. B2 simplifies to:

 B(psi) = (1 + cos(psi)^2) 3 / (16 PI).                                     (C1)

 And Eq. B3 and B4 simplify to:

 PDF(psi) = 2 PI B(psi) sin(psi) ->
          = 0.375 (1 + cos(psi)^2) sin(psi),                                (C2)

 CDF(psi) = Integral(0 to psi){dpsi PDF(psi)} ->
          = -0.125 cos(psi)^3 - 0.375 cos(psi) + 0.5                        (C3)

 And the quantile function of Eq. B5 simplifies to:

 psi   =  acos(delta - 1 / delta),                                          (C4)

 delta = cbrt(sqrt(16 CDF(psi)^2 - 16 CDF(psi) + 5) - 4 CDF(psi) + 2).      (C5)

 INPUT:
 scat   - Scattering parameters
          Pointer to constant scattering struct; 
 s_scat - Spherical direction
          [0]Theta, [1]Phi, [2]Radius
          Range: [0,PI] radians, [0,2PI] radians, 1.0 unitless
          Pointer to constant double or pointer to double;
 prob   - Cumulative probability associated with a polar and a azimuth angle
          [0]Theta, [1]Phi
          Range: [0,1] unitless
          Pointer to constant double or pointer to double;
 dens   - Probability density associated with a polar and a azimuth angle
          [0]Theta, [1]Phi
          Range: [0,Inf) unitless
          Pointer to double
 pf     - Phase function at a particular angle
          Range: [0,Inf), sr
          Pointer to double;

 OUTPUT:
 None. Updates values pointed by s_scat, prob or pf.

*******************************************************************************/

 void
 scat_qtl_rlg_air
 (
   struct scattering const * scat,
   double * s_scat,
   double const * prob
 )
 {
   double delta = cbrt(sqrt(16.0 * prob[0] * prob[0] - 16.0 * prob[0] + 5.0) - 
     4.0 * prob[0] + 2.0);
   s_scat[0] = acos(delta - 1.0 / delta);
   s_scat[1] = K_2PI * prob[1];
 }

 void
 scat_cdf_rlg_air
 (
   struct scattering const * scat,
   double const * s_scat,
   double * prob
 )
 {
   double const mu = cos(s_scat[0]);
   prob[0] = -(0.125 * mu * mu * mu) - (0.375 * mu) + 0.5;
   prob[1] = K_1_2PI * s_scat[1];
 }

 void
 scat_pdf_rlg_air
 (
   struct scattering const * scat,
   double const * s_scat,
   double * dens
 ) 
 {
   double const mu = cos(s_scat[0]);
   dens[0] = 0.375 * (1.0 + mu * mu) * sin(s_scat[0]);
   dens[1] = K_1_2PI;
 }

 void
 scat_pf_rlg_air
 (
   struct scattering const * scat,
   double const * s_scat,
   double * pf
 )
 {
   double const mu = cos(s_scat[0]);
   pf[0] = RAYLEIGH_A_NORM * (1.0 + mu * mu);
 }


