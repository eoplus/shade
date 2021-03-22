/********************************************************************************
 *										*
 * trunc.c									* 
 *										*
 *	Implements geometrical truncation of a phase function (PF). Higly 	*
 * anisoropic phase functions can cause high variability with the local 	*
 * estimate variance reduction technique, not only requiring increased number 	*
 * of simulations for convergency but also causing random spikes in the the 	*
 * convergency series (c.f. Buras & Mayer 2011). Amoung the methods to reduce 	*
 * this effect is the truncation of the diffraction peak (small angle 		*
 * scattering; e.g. Iwabuchi & Suzuki, 2009). The function presented here 	*
 * perform a geometrical truncation of the PF at an angle of 5º. The PF at 	*
 * angles < 5º is substituted by the average sine weighted PF between 0 and 5º:	*
 *										*
 * PF(psi < 5º) = k								*
 * k = integral PF(psi)*sin(psi) dpsi / integral sin(psi) dpsi			*
 * CDF(psi < 5º) = k * 2 * PI * integral sin(psi) dpsi				*
 *										*
 *	With those definitions, the normalization condition for PF is satisfied *
 * and the RTE is unchanged. The truncation might introduce a bias since the 	*
 * angular distribution is now altered. However, for thick optical depth medium *
 * as the water, where multiple scattering is important, the scpecific angular 	*
 * shape of the PF loses importance for the final solution (Nakajima & Tanaka, 	*
 * 1988; Gordon 1993; Mobley et al. 2002). Also, the intermediate and 		*
 * backscatter angles are left unchanged, so reducing much of any possible bias *
 * (Mobley et al. 2002). The truncation was validated for Lu between 0 and 80º 	*
 * SZA and for diffuse light for the Petzold PF and agrees with Hydrolight 	*
 * (version 4.2) without any significant bias. It is expected but not garanteed *
 * to be the same outside those conditions.					*
 *										*
 * References:									*
 * Buras, R. & Mayer, B. 2011. Efficient unbiased variance reduction		*
 *	techniques for Monte Carlo simulations of radiative transfer in 	*
 *	cloudy atmospheres: The solution. J. Quant. Spectrosc. Radiat. Transfer *
 *	112, 3, 434-447.							*
 *	DOI: 10.1016/j.jqsrt.2010.10.005					*
 * Iwabuchi, H. & Suzuki, T. 2009. Fast and accurate radiance calculations 	*
 *	using truncation approximation for anisotropic scattering phase 	*
 *	functions. J. Quant. Spectrosc. Radiat. Transfer 110, 17, 1926-1939.	*
 *	DOI: 10.1016/j.jqsrt.2009.04.006					*
 * Gordon, H. R. 1993. Sensitivity of radiative transfer to small-angle 	*
 *	scattering in the ocean: Quantitative assessment. Applied Optics 32, 	*
 *	36, 7505-7511.								*
 *	DOI: 10.1364/AO.32.007505						* 
 * Mobley, C. D.; Sundman, L. K. & Boss, E. S. 2002. Phase function effects on 	*
 *	oceanic light fields. Applied Optics 41, 6, 1035-50.			*
 *	DOI: 10.1364/AO.41.001035						*
 * Nakajima, T. & Tanaka, M. 1988. Algorithms for radiative intensity 		*
 *	calculations in moderately thick atmospheres using a truncation 	*
 *	approximation. J. Quant. Spectrosc. Radiat. Transfer, 40, 51–69.	*
 *	DOI: 10.1016/0022-4073(88)90031-3					*
 *										*
 ********************************************************************************/

 #include <stdio.h>	// Functions: printf
 #include <stdlib.h>	// Functions: exit
 #include <math.h>	// Constants: M_PI

 #include "constants.h"	// Constants: RAD, DM_PI
 #include "scattering.h"
 #include "scattering_ff.h"	// Functions: cdf_ff, pf_ff

/* Truncation of PF and CDF for Fournier-Forand: *******************************
 
 scat_trunc_ff 

 INPUT:

 OUTPUT:
 Update the pf and cdf vectors to be used in GSL interpolation functions.
 
*******************************************************************************/

 void
 scat_trunc_ff
 ( struct scattering *scat )
 {
   int ci = 0;									// General counting variable
   int ci_cpsi = 0;								// Index of the critical angle for truncation
   double cpsi = 5 * RAD;							// Critical angle for truncation
   double k = 0.0;								// Constant for angles < cpsi
   double step_ff = M_PI / (double) (scat->lut.npsi - 1);				// Step increase of the psi points
   double s_scat[2] = {0.0};
   double prob[2] = {0.0};

   s_scat[0] = cpsi;
   scat_cdf_ff (scat, s_scat, prob);
   k = (prob[0] * K_1_2PI) / (1.0 - cos(cpsi));		// 1.0 - cos(psi) is the analytical integral of sin(psi) dpsi
  
   ci_cpsi = (int) trunc (cpsi / step_ff);
  
   for (ci = 0; ci <= ci_cpsi; ci++)
   {
     scat->lut.pf[ci] = k;
     scat->lut.cdf[ci] = k * (1.0 - cos(scat->lut.psi[ci])) * K_2PI;
   }
 }

/* Old version: 
 void trunc_ff (int np_ff, double* psi_ff, double* pfv_ff, double* itg_ff, 
                double j, double nh) {
  int ci = 0;									// General counting variable
  int ci_cpsi = 0;								// Index of the critical angle for truncation
  double cpsi = 5 * RAD;							// Critical angle for truncation
  double k = 0.0;								// Constant for angles < cpsi
  double step_ff = M_PI / (double) (np_ff - 1);					// Step increase of the psi points
  
  k = (cdf_ff (cpsi, j, nh) / (2.0 * M_PI)) / (1.0 - cos(cpsi));		// 1.0 - cos(psi) is the analytical integral of sin(psi) dpsi
  
  ci_cpsi = (int) trunc (cpsi / step_ff);
  
  for (ci = 0; ci <= ci_cpsi; ci++) {
      pfv_ff[ci] = k;
      itg_ff[ci] = k * (1.0 - cos(psi_ff[ci])) * 2.0 * M_PI;
  }
 }
*/

/*
 * Truncation of PF and CDF for tabulated function
 * 
 * NOT IMPLEMENTED IN THIS VERSION!
 *
 */

 void scat_trunc_pf (void) {
  printf("Truncation for tabulated phase functions not implemented in \
         this version\n");
  exit(-1);
 }
 

