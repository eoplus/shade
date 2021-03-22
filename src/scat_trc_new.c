/* Truncation of a phase function (PF)
 *
 * Large particles have strong forward scattering due to the large diffraction 
 * signal and this effect can create large variability in Monte Carlo 
 * simulations when the local estimate method is used (e.g., Burras & Mayer, 
 * 2011). The local estimate method creates a virtual ray out of the current ray 
 * position and calculates the probability of reaching a given spatially and / 
 * or directionally finite target. The extreme example of such is a spatially 
 * finite radiance sensor with a small FOV in a 3D radiative transfer problem. 
 * This is analogous to the Sun, in backward Monte Carlo simulations, where due 
 * to its distance from Earth, there is a very narrow range of directions that 
 * that point directly from the surface to the Sun.
 * 
 * A common solution to avoid the large variability associated with the 
 * diffraction peak is the truncation of the PF on the diffraction peak, while
 * adjusting the IOPs accordingly. However, since the PF is altered even in a 
 * small set of near forward angles, some specific effects (like glory) may be 
 * affected. However, sophisticated truncation methods, such as the adaptative 
 * truncation (Iwabuchi & Suzuki, 2009) have very low bias. Unbiased methods
 * also exist (Burras & Mayer, 2011) but are computationally more intensive. 
 * Additionally, the magnitude of the problem and the bias is proportional to 
 * how peaked the PF is, and this is really important for computations for 
 * cloudy atmospheres.
 *
 * Here we implement the method of adaptative PF truncation.
 * 
 * Lets recall that in this code, the phase functions (beta) are normalized such 
 * that:
 *
 * 2 PI Integral(0 to PI){beta(psi) sin(psi) dpsi} = 1.                      (1)
 *
 * If the PF has a strong forward scattering peak, it can be useful to consider 
 * that all rays scattered at a near 0 angle were not scattered, since they have 
 * approximately the same direction. In this case, the PF can be represented as 
 * the following sum:
 * 
 * beta(psi) = beta'(psi) + Dirac(0 - psi) (1 - f) / (2 PI),                 (2)
 *
 * where Dirac is the Dirac delta function and with:
 *
 * beta'(psi) = { k          , if psi <  psi_trunc
 *              { beta(psi)  , if psi >= psi_trunc,                          (3)
 *
 * and:
 *
 * f = 2 PI Integral(0 to PI){beta'(psi) sin(psi) dpsi}.                     (4)
 *
 * Note that f, as defined by Eq. 4, will be < 1 if:
 * 
 * k 2 PI Integral(0 to psi_trunc){sin(psi) dpsi} < 
 *   2 PI Integral(0 to psi_trunc){beta(psi) sin(psi) dpsi},                 (5)
 *
 * and the difference 1 - f is the scattering that was removed from the 
 * diffraction peak. Therefore, in Eq. 2, Dirac(0 - psi) (1 - f) / (2 PI) 
 * attributes the sine weighted beta for psi < psi_trunc to psi = 0. Since 
 * sin(0) = 0, the normalization condition in Eq. 1 is no longer fulfilled as 
 * this fraction of the total normalized scattering will not contribute to the 
 * integral. In this situation, the IOPs must compensate for the missing 
 * fractional scattering, and the RTE remains unchanged with the following 
 * substitutions (Gordon, 1985):
 *
 * c'  = c (1 - w0 (1 - f)),                                                 (6)
 * w0' = w0 (1 - (1 - f)) / (1 - w0 (1 - f)).                                (7)
 *
 * That is, all scattering that was removed, does not contribute to the beam 
 * attenuation, reducing it and reducing the single scattering albedo. We can 
 * show with a simple example that this is the case. Following from the two flow 
 * equations and with a series of simplifications that are appropriate for 
 * water, the hemispherical-directional reflectance can be approximated as:
 *
 * K b bbr / (a + b bbr),                                                    (8) 
 * 
 * K is an unspecified constant, b is the scattering coefficient, bbr is the 
 * backscattering ratio and a is the absorption coefficient. Then:
 *
 * K b' bbr' / (a + b' bbr') ->
 * K (b f) (bbr / f) / (a + (b f) (bbr / f)) ->
 * K b bbr / (a + b bbr).
 *
 * Note that bbr' = bbr / f because with the removal of the diffraction peak, 
 * the backscatter represents a higher proportion of the total scattering.
 *
 * The constant value k in Eq. 3 may be defined in different ways. Here we will
 * define it in a way to preserve the 0th and the 1st moments of the phase 
 * function (that is, the the integral over all directions is 1 and the 
 * expectation (average angle) is unchanged). The 0th moment is preserved if
 * Eq. 5 is an equality and the 1st moment if:
 *
 * k 2 PI Integral(0 to psi_trunc){psi sin(psi) dpsi} = 
 *   2 PI Integral(0 to psi_trunc){psi beta(psi) sin(psi) dpsi}.             (9)
 *
 * Then b, c and w0 remain unchanged and we minimize changes of the 
 * directionality of the PF.
 *
 * The specific method of truncation is the adaptative method of Iwabuchi & 
 * Suzuki (2009). In this method, the PF is smoothed progressively, with the
 * order and direction of scatter. The principle is that with multiple scatter
 * the details of the PF are not important. For the first scatter, the real beta
 * beta is used. If the scattering angles is higher that psi_trc_1, then the 
 * beta is updated for the next scattering event, by making a small truncation,
 * resulting in beta_1. And so it continues until a limit. Here 4 scattering 
 * events at angles > psi_trc_n are necessary to reach the full truncation of 
 * the diffraction peak. The maximum truncation angle is defined a priori to be 
 * 7 degrees.
 * 
 */

 #include "constants.h"	// Constants: RAD

 #define PSI_TRC_MAX 7 * RAD

 struct scat_trc_par
 {
   double psi_trc;
   double beta_k[3];
   int    id;
 };

 scat
 {
   struct scat_trc_par* par = (struct scat_trc_par*) scat_par;
   if ( psi < par->psi_trunc )
   {
     return par->beta_k[par->id];
   } else {
     par->id += 1;  //This is not good as it will always be increasing... should stop at 3 
     par->psi_trunc += par->step;

     
   }
 }
