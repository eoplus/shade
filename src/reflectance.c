
/*******************************************************************************
 reflectance.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 This source file provides a representation of Bidirectional Reflectance
 Distribution Function (BRDF) and functions to be used with the Monte Carlo
 code.

 The BRDF can be used as bi-variate probability distribution function (PDF) to
 sample reflected directions dependent on the incident directions. The BRDF is
 defined as (Nikodemus, 1977):

 BRDF(t_i, p_i, t_r, p_r) =
   L(t_r, p_r) / (L(t_i, p_i) cos(t_i) omega(t_i, p_i)),                     (1)
 
 where L is radiance, omega is a solid angle centered on the directions of 
 incidence, t is the polar angle and p is the azimuth, with the subscript i for 
 the incident ray and r for the reflected one. 

 The directional-hemispherical reflectance is the integral of the BRDF over
 all reflected directions and represent the probability that a ray incident in 
 a given direction will be reflected:
 
 rho_dh(t_i, p_i) = 
   Integral{domega(t_r, p_r) BRDF(t_i, p_i, t_r, p_r) cos(t_r)}.             (2)

 By virtue of the Helmholtz reciprocity principle, this is the same as the 
 hemispherical-directional reflectance, if the incident light distribution is 
 uniform. To isolate directional aspects from the reflectivity, it is common 
 to present a normalized BDRF, defined as:
 
 BRDF_n(t_i, p_i, t_r, p_r) = BRDF(t_i, p_i, t_r, p_r) / rho_dh(t_i, p_i).   (3)

 The bi-hemispherical reflectance integrates over all incident and reflected 
 directions, and therefore is an apparent optical property as it depends on 
 the incident light distribution. It can be written as:

 rho = Eu / Ed ->
     = Integral{domega(t_r, p_r) Lu cos(t_r)} / 
         Integral{domega(t_i, p_i) Ld cos(t_i)} ->
     = Integral{domega(t_r, p_r) Integral{domega(t_i, p_i) Ld 
         BRDF(t_i, p_i, t_r, p_r) cos(t_i)} cos(t_r) } /
         Integral{domega(t_i, p_i) Ld cos(t_i)},                             (4)
  
 where the integrals are over the incidence / reflectance hemisphere.
 Simplification of the above equation is only possible in some simple conditions
 of illumination and BRDF.

 The Cumulative Distribution Function (CDF) for the azimuthal angle is given 
 by integrating Eq. 2 from 0 up to the azimuth of interest and normalizing by 
 the full integral (rho_dh(t_i, p_i)). The CDF for the polar angle for a given
 incidence and reflected azimuth is given by the integral from 0 to t_r,
 normalized by the full integral from 0 to PI/2, in both cases the azimuth is
 fixed.
  
 For each BRDF model, there are 8 functions, though only 5 are actually needed
 in the Monte Carlo solver. Those functions have the following name 
 convention:
 
 brdf_build_XXX   - Setup the parameter structure of the BRDF model;
 brdf_XXX         - Returns a pointer to the BRDF vector at the requested 
                     geometry;
 brdf_dhr_XXX     - Returns a pointer to the DHR vector at the requested 
                     incidence angle;
 brdf_pdf_phi_XXX - Returns the PDF for phi_r at the requested geometry;
 brdf_cdf_phi_XXX - Returns the CDF for phi_r at the requested geometry;
 brdf_pdf_psi_XXX - Returns the PDF for psi_r at the requested geometry;
 brdf_cdf_psi_XXX - Returns the CDF for psi_r at the requested geometry;
 brdf_qtl_XXX     - Updates the values of phi_r and psi_r for the requested 
                    quantiles of phi_r and psi_r.

 As elsewhere in the code, the quantile function can be used to perform random
 sampling of the PDFs for phi_r and psi_r at a given incidence angle, by 
 passing random uniform deviates as the probabilities.
  
 The parameters for each BRDF function are made generic by passing a void 
 pointer and internally defining the parameter structure which is unique to 
 each BRDF model. All parameter structures are named brdf_par_XXXX, where XXXX
 is the name of the BRDF model and they are defined in the brdf.h header file.
 
 The inputs of the BRDF functions, excluding the setup functions, are 
 described below. The outputs were described above.
  
 INPUT:
 mut_i      - (double)  Incident cosine polar angle;
 phi_i      - (double)  Incident azimuth angle;
 mut_r      - (double)  Reflected cosine polar angle;
 phi_r      - (double)  Reflected azimuth angle;
 p_mut_r    - (double*) Pointer to reflected cosine polar angle;
 p_phi_r    - (double*) Pointer to reflected azimuth angle;
 p_brdf_par - (void*)   Pointer to BRDF parameter structure; 
 prob_psi   - (double)  CDF probability for psi;
 prob_phi   - (double)  CDF probability for phi. 
  
*******************************************************************************/

 #include <stdlib.h>
 #include <stdio.h>
 #include <string.h>
 #include <math.h>

 #include "config.h"
 #include "memory.h"
 #include "aux.h"
 #include "constants.h"
 #include "reflectance.h"


 struct reflectance*
 refl_alloc
 ( void )
 {
   struct reflectance *refl =
     (struct reflectance*) malloc(sizeof(struct reflectance));
   #ifdef CHCK_MEM_ALLOC
   if ( !refl )
   {
     printf("\nERROR: Could not allocate memory for refl in refl_alloc"
       " (reflectance.c)\n");
     exit(-1);
   }
   #endif // CHCK_MEM_ALLOC

   // Allocated if necessary by the setup functions:
//   refl->minn = NULL;
//   refl->lamb = NULL;

   return refl;
 }

 void
 refl_free
 ( struct reflectance** refl )
 {
//   free( (*refl)->fun );
//   if ( (*refl)->minn ) refl_free_minn( &(*refl)->minn );
//   if ( (*refl)->lamb ) refl_free_lamb( &(*refl)->lamb );
   free( *refl );
   *refl = NULL;
 }

 void
 refl_setup
 (
   struct reflectance *refl,
   char const   *tp,
   int const    nbhr,
   double const *bhr,
   double const k
 )
 {
   strncpy(refl->tp, tp, STRMXLEN);
   if ( strcmp(tp, "minnaert") == 0 )
   {
     // Minnaert:
     refl_setup_minn(refl, nbhr, bhr, k);
     return;
   }
   else if (strcmp(tp, "lambert") == 0 )
   {
     // Lambert:
     refl_setup_lamb(refl, nbhr, bhr);
     return;
   }
 }

/* Print reflectance properties: ***********************************************

 refl_fprintf

 Prints scattering parameters.

 INPUT:
 refl - Reflectance parameters
        Pointer to pointer to constant reflectance struct;

 OUTPUT:
 None. Prints the parameters of the reflectance function.

*******************************************************************************/

 void
 refl_fprintf
 (
   FILE *odv,
   struct reflectance const *refl,
   int const indent
 )
 {
   char pre_0[STRMXLEN] = "";
   char pre_1[STRMXLEN] = "";
   for (size_t i = 0; i < indent; i++)
     strcat(pre_0, "  ");
   strncpy(pre_1, pre_0, STRMXLEN);
   strcat(pre_1, "  ");

   fprintf(odv, "%sReflectance parameters:\n", pre_0);
   if ( strcmp(refl->tp, "minnaert") == 0 )
   {
     fprintf(odv, "%sType: %s (k = %.3lf)\n", pre_1, refl->tp, refl->minn.k);
     fprintf(odv, "%sN.:   %d\n", pre_1, refl->minn.nbhr);
     fprintf(odv, "%sBHR:  ", pre_1);
     for (size_t i = 0; i < (refl->minn.nbhr - 1); i++)
     {
       if ( i && !(i % 5) )
         fprintf(odv, "\n%s      ", pre_1);
       fprintf(odv, "%.4lf, ", refl->minn.bhr[i]);
     }
     fprintf(odv, "%.4lf (unitless)\n", refl->minn.bhr[(refl->minn.nbhr - 1)]);
   }
   else if ( strcmp(refl->tp, "lambert") == 0 )
   {
     fprintf(odv, "%sType: %s (Minnaert with k = 0)\n", pre_1, refl->tp);
     fprintf(odv, "%sN.:   %d\n", pre_1, refl->lamb.nbhr);
     fprintf(odv, "%sBHR:  ", pre_1);
     for (size_t i = 0; i < (refl->lamb.nbhr - 1); i++)
     {
       if ( i && !(i % 5) )
         fprintf(odv, "\n%s      ", pre_1);
       fprintf(odv, "%.4lf, ", refl->lamb.bhr[i]);
     }
     fprintf(odv, "%.4lf (unitless)\n", refl->lamb.bhr[(refl->lamb.nbhr - 1)]);
   }
 }


/* Minnaert BRDF model: ********************************************************

 refl_alloc_minn, refl_free_minn, refl_setup_minn, refl_qtl_minn, refl_cdf_minn,
 refl_pdf_minn, refl_dhr_minn, refl_brdf_minn

 The Minnaert BRDF model is defined as:

 BRDF(ti,pi,tr,pr) = (BHR / PI) (cos(ti) cos(tr))^k                          (1)

 With the integral for tr from 0 to PI/2 and the integral for pr from 0 to 
 2 PI, the DHR(ti,pi) is then:

 DHR(ti,pi) = Integrate{Integrate{BRDF(ti,pi,tr,pr) cos(tr) sin(tr) dtr} dpr} ->
 = Integrate{Integrate{(BHR / PI) cos(ti)^k cos(tr)^k cos(tr) sin(tr) dtr} dpr} ->
 = (BHR / PI) cos(ti)^k Integrate{Integrate{cos(tr)^k cos(tr) sin(tr) dtr} dpr} ->
 = 2 BHR cos(ti)^k Integrate{cos(tr)^k cos(tr) sin(tr) dtr} ->
 = 2 BHR cos(ti)^k / (k + 2).                                                (2)

 With the integral from 0 to PI/2, the PDF for tr is:

 PDF(ti,pi,tr) = BRDF(ti,pi,tr,pr) cos(tr) sin(tr) / Integrate{BRDF(ti,pi,tr,pr) cos(tr) sin(tr) dtr} ->
 = (R / PI) (cos(ti) cos(tr))^k cos(tr) sin(tr) / Integrate{(R / PI) (cos(ti) cos(tr))^k cos(tr) sin(tr) dtr} ->
 = (R / PI) cos(ti)^k cos(tr)^(k+1) sin(tr) / (R / PI) cos(ti)^k Integrate{cos(tr)^(k+1) sin(tr) dtr} ->
 = cos(tr)^(k+1) sin(tr) / (1 / (k + 2)) ->
 = (k + 2) cos(tr)^(k+1) sin(tr).                                          (3)

 The CDF for tr is:
 
 CDF(ti,pi,tr) = Integrate(from 0 to tr){PDF(ti,pi,pr) dtr} ->
 = Integrate(from 0 to pr){(k + 2) cos(tr)^(k+1) sin(tr) dtr} ->
 = (k + 2) Integrate(from 0 to pr){cos(tr)^(k+1) sin(tr) dtr} ->
 = 1 - ((k + 2) cos(tr)^(k+2) / (k + 2)) ->
 = 1 - cos(tr)^(k+2)                                                       (4)

 With the integral from 0 to PI/2 for tr, the PDF for pr is:

 PDF(ti,pi,pr) = (1 / DHR(ti,pi)) Integrate{BRDF(ti,pi,tr,pr) cos(tr) sin(tr) dtr} ->
 = (1 / DHR(ti,pi)) Integrate{BRDF(ti,pi,tr,pr) cos(tr) sin(tr) dtr} ->
 = (1 / DHR(ti,pi)) Integrate{(R / PI) (cos(ti) cos(tr))^k cos(tr) sin(tr) dtr} ->
 = (1 / DHR(ti,pi)) (R / PI) cos(ti)^k Integrate{cos(tr)^(k+1) sin(tr) dtr} ->
 = (1 / DHR(ti,pi)) (R / PI) cos(ti)^k / (k + 2) ->
 = ((k + 2) / (2 R cos(ti)^k)) (R / PI) cos(ti)^k / (k + 2) ->
 = ((k + 2) (R / PI) cos(ti)^k) / ((k + 2) 2 R cos(ti)^k) ->
 = 1 / (2 PI)                                                              (5)
 
 The CDF for pr is then:

 CDF(ti,pi,pr) = Integrate(from 0 to pr){PDF(ti,pi,pr) dpr} ->
 = Integrate(from 0 to pr){dpr / (2 PI)} ->
 = pr / (2 PI)                                                             (6)

 The quantile function for tr and pr is found by solving the CDFs for tr and 
 pr:

 CDF(ti,pi,tr) = 1 - cos(tr)^(k+2) ->
 cos(tr)^(k+2) = 1 - CDF(ti,pi,tr) ->
 cos(tr) = (1 - CDF(ti,pi,tr))^(1/(k + 2)) ->
 tr = acos((1 - CDF(ti,pi,tr))^(1 / (k + 2)))                              (7)

 CDF(ti,pi,pr) = pr / (2 PI) ->
 pr = 2 PI CDF(ti,pi,pr)                                                   (8)

*******************************************************************************/
 
 struct reflectance_minn*
 refl_alloc_minn
 ( void )
 {
   return (struct reflectance_minn*) malloc(sizeof(struct reflectance_minn));
 }

 void
 refl_free_minn
 ( struct reflectance_minn **minn )
 {
   free( *minn );
   *minn = NULL;
 }

 void
 refl_setup_minn
 (
   struct reflectance *refl,
   int const    nbhr,
   double const *bhr,
   double const k
 )
 {
   refl->minn.nbhr = nbhr;
   refl->minn.bhr  = calloc_1d (nbhr, 
     "refl->minn->bhr in refl_setup_minn (reflectance.c)");
   refl->minn.bdr  = calloc_1d (nbhr, 
     "refl->minn->bdr in refl_setup_minn (reflectance.c)");

   for (size_t i = 0; i < nbhr; i++)
   {
     refl->minn.bhr[i]  = bhr[i];
     refl->minn.bdr[i]  = bhr[i] * M_1_PI;
   }

   refl->minn.k      = k;
   refl->minn.kp1    = k + 1.0;
   refl->minn.kp2    = k + 2.0;
   refl->minn.invkp2 = 1.0 / (k + 2.0);

   refl->fun.qtl  = refl_qtl_minn;
   refl->fun.cdf  = refl_cdf_minn;
   refl->fun.pdf  = refl_pdf_minn;
   refl->fun.dhr  = refl_dhr_minn;
   refl->fun.brdf = refl_brdf_minn;
 }

 void
 refl_qtl_minn
 (
   struct reflectance const *refl,
   double const *s_scat_i,
   double       *s_scat_r,
   double const *prob
 )
 {
   s_scat_r[0] = pow(1.0 - prob[0], refl->minn.invkp2); // COSINE!
   s_scat_r[1] = prob[1] * K_2PI;
 }

 void
 refl_cdf_minn
 (
   struct reflectance const *refl,
   double const *s_scat_i,
   double const *s_scat_r,
   double       *prob
 )
 {
   prob[0] = 1.0 - pow(s_scat_r[0], refl->minn.kp2);
   prob[1] = s_scat_r[1] * K_1_2PI;
 }

 void
 refl_pdf_minn
 (
   struct reflectance const *refl,
   double const *s_scat_i,
   double const *s_scat_r,
   double       *dens
 )
 {
   dens[0] = refl->minn.kp2 * pow(s_scat_r[0], refl->minn.kp1) * 
     sqrt(1.0 - s_scat_r[0] * s_scat_r[0]);
   dens[1] = K_1_2PI;
 }

 void
 refl_dhr_minn
 (
   struct reflectance const *refl,
   double const *s_scat_i,
   double       *dhr
 )
 {
   double scl = 2.0 * pow(s_scat_i[0], refl->minn.k) * refl->minn.invkp2;
   for (size_t i = 0; i < refl->minn.nbhr; i++)
   {
     dhr[i] = refl->minn.bhr[i] * scl;
   }
 }

 void
 refl_brdf_minn
 (
   struct reflectance const *refl,
   double const *s_scat_i,
   double const *s_scat_r,
   double       *bdr
 )
 {
   double scl = pow((s_scat_i[0] * s_scat_i[1]), refl->minn.k);
   for (size_t i = 0; i < refl->minn.nbhr; i++)
   {
     bdr[i] = refl->minn.bdr[i] * scl;
   }
 }

/* Lambert BRDF model: *********************************************************

 Though a special case of the Minnaert BRDF when the k parameter is equal to
 0, it is more representative as a model surface and a direct implementation 
 is more efficient. The Lambertian model results in a constant BRDF by having 
 a reflectance in a element of surface that is proportionally decreasing with 
 the cosine of the angle, compensating the proportionally increasing area with 
 the cosine of the angle.
 
 rho_dh = Integral{BRDF(ti,pi,tr,pr) cos(theta) sin(theta) dtheta dphi} -> 
        = BRDF * pi.

 Since the BRDF is constant, rho_dh does not dependent on direction either and
 is equal to rho. The equations below are simplifications arising when k = 0 
 and need not an explicit derivation.

*******************************************************************************/

 struct reflectance_lamb*
 refl_alloc_lamb
 ( void )
 {
   return (struct reflectance_lamb*) malloc(sizeof(struct reflectance_lamb));
 }

 void
 refl_free_lamb
 ( struct reflectance_lamb **lamb )
 {
   free( *lamb );
   *lamb = NULL;
 }

 void
 refl_setup_lamb
 (
   struct reflectance *refl,
   int const    nbhr,
   double const *bhr
 )
 {
   refl->lamb.nbhr = nbhr;
   refl->lamb.bhr  = calloc_1d(nbhr,
     "refl->lamb.bhr in refl_setup_lamb (reflectance.c)");
   refl->lamb.bdr  = calloc_1d(nbhr,
     "refl->lamb.bdr in refl_setup_lamb (reflectance.c)");
   for (size_t i = 0; i < nbhr; i++)
   {
     refl->lamb.bhr[i] = bhr[i];
     refl->lamb.bdr[i] = bhr[i] * M_1_PI;
   }

   refl->fun.qtl  = refl_qtl_lamb;
   refl->fun.cdf  = refl_cdf_lamb;
   refl->fun.pdf  = refl_pdf_lamb;
   refl->fun.dhr  = refl_dhr_lamb;
   refl->fun.brdf = refl_brdf_lamb;
 }

 void
 refl_qtl_lamb
 (
   struct reflectance const *refl,
   double const *s_scat_i,
   double       *s_scat_r,
   double const *prob
 )
 {
   s_scat_r[0] = sqrt(prob[0]); // COSINE!
   s_scat_r[1] = prob[1] * K_2PI;
 }

 void
 refl_cdf_lamb
 (
   struct reflectance const *refl,
   double const *s_scat_i,
   double const *s_scat_r,
   double       *prob
 )
 {
   prob[0] = 1.0 - s_scat_r[0] * s_scat_r[0];
   prob[1] = s_scat_r[1] * K_1_2PI;
 }

 void
 refl_pdf_lamb
 (
   struct reflectance const *refl,
   double const *s_scat_i,
   double const *s_scat_r,
   double       *dens
 )
 {
   dens[0] = 2.0 * s_scat_r[0] * sqrt(1.0 - s_scat_r[0] * s_scat_r[0]);
   dens[1] = K_1_2PI;
 }

 void
 refl_dhr_lamb
 (
   struct reflectance const *refl,
   double const *s_scat_i,
   double       *dhr
 )
 {
   for (size_t i = 0; i < refl->lamb.nbhr; i++)
     dhr[i] = refl->lamb.bhr[i];
 }

 void
 refl_brdf_lamb
 (
   struct reflectance const *refl,
   double const *s_scat_i,
   double const *s_scat_r,
   double       *bdr
 )
 {
   for (size_t i = 0; i < refl->lamb.nbhr; i++)
     bdr[i] = refl->lamb.bdr[i];
 }



/* Unidirectional BRDF model
 *
 * This is a simple construct to allow specific test calculations, not intended 
 * for general simulation. It will always reflect in a given (fixed) direction, 
 * independent of the incident direction. Additionally, its BRDF from all 
 * directions to that direction is constant.
 * 
 */
/*
 void brdf_build_unid (void  **pp_brdf_par,
                       const int    nbhr,
                       const double *p_bhr,
                       const double thetaf,
                       const double phif)
 {
   int i;
   *pp_brdf_par = (void*) calloc(1, sizeof(struct brdf_par_unid));

   struct brdf_par_unid *p_par = (struct brdf_par_unid*) *pp_brdf_par;
   p_par->nbhr    = nbhr;
   p_par->p_bhr   = (double*) calloc (nbhr, sizeof(double));
   p_par->p_dhr   = (double*) calloc (nbhr, sizeof(double));
   p_par->p_brdf  = (double*) calloc (nbhr, sizeof(double));
   p_par->p_brdf0 = (double*) calloc (nbhr, sizeof(double));
   for (i = 0; i < nbhr; i++)
   {
     p_par->p_bhr[i]  = p_bhr[i];
     p_par->p_dhr[i]  = p_bhr[i];
     p_par->p_brdf[i] = p_bhr[i] / M_PI;
   }
   p_par->mutf   = cos(thetaf);
   p_par->phif   = phif;
 } 

 double* brdf_unid (const double mut_i,
                    const double phi_i,
                    const double mut_r,
                    const double phi_r,
                    void  *p_brdf_par)
 {
   struct brdf_par_unid *p_par = (struct brdf_par_unid*) p_brdf_par;
   if ( NUM_EQU(mut_r, p_par->mutf, 1E-12) && 
        NUM_EQU(phi_r, p_par->phif, 1E-12) )
   {
     return p_par->p_brdf;
   } else {
     return p_par->p_brdf0;
   }
 }

 double* brdf_dhr_unid (const double mut_i,
                        const double phi_i,
                        void  *p_brdf_par)
 {
   return ((struct brdf_par_unid*) p_brdf_par)->p_dhr;
 }

 double brdf_cdf_phi_unid (const double mut_i,
                           const double phi_i,
                           const double mut_r,
                           const double phi_r,
                           const void   *p_brdf_par)
 {
   return (phi_r >= ((struct brdf_par_unid*) p_brdf_par)->phif) ? 1.0 : 0.0;
 }

 double brdf_cdf_psi_unid (const double mut_i,
                           const double phi_i,
                           const double mut_r,
                           const double phi_r,
                           const void   *p_brdf_par)
 {
   return (mut_r >= ((struct brdf_par_unid*) p_brdf_par)->mutf) ? 1.0 : 0.0;
 }

 void brdf_qtl_unid (const double mut_i,
                     const double phi_i,
                     double       *p_mut_r,
                     double       *p_phi_r,
                     const void   *p_brdf_par,
                     const double prob_psi,
                     const double prob_phi)
 {
   *p_mut_r = ((const struct brdf_par_unid*) p_brdf_par)->mutf;
   *p_phi_r = ((const struct brdf_par_unid*) p_brdf_par)->phif;
 }

*/




