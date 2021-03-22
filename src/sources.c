
/*******************************************************************************
 sources.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 Provides functions describing light sources and to take random samples of their
 probability distribution function.

 Light sources can be of various types. The Sun as seen from Earth is a source
 of nearly colimated unpolarized light. Lasers are colimated sources,
 transpectral processes can be taken as point isotropic unpolarized sources and
 in Backward Monte Carlo, sensors (radiance, plane irradiance, scalar
 irradiance) are taken as sources.

 All sources are described by a set of properties presented in detail in the
 sources.h header file. In summary, they are described by an origin (center),
 radius (for its area, currently = 0) and orientation axis. 

 THEORY:
 For reasons discussed later, the theory here deals with point sources
 (infinitesimal area extent) only. As for the scattering phase functions, we
 normalize the emission functions such that the angular integral equals unit:

 Integral{ dphi
   Integral{ dpsi
     S(psi,phi) sin(psi) } } = 1,                                            (1)

 where the subscripts l and u on the integration limits stand for lower and
 upper polar angle of emission. We use the convention that if an integral does
 not specifies bounds, the integration is to be carried over the full range of
 the variable. Those ranges are [0,2PI] for the azimuth angle phi and [0,PI] for
 the polar angle psi. The normalization condition results in the Emission
 Function (EF), S(psi,phi), being a probability function, with Eq. 1 being the
 Cumulative Distribution Function (CDF).

 It is possible to separate this bi-dimensional probability function into the 
 probabilities of its components for the polar and azimuth angles:

 PDF(psi) = S(psi,phi) sin(psi) / Integrate{S(psi,phi) sin(psi) dpsi},       (2)

 CDF(psi) = Integrate(0 to psi){PDF(psi) dpsi},                              (3)

 PDF(phi) = Integrate{S(psi,phi) sin(psi) dpsi} / CDF(psi,phi),              (4)

 CDF(phi) = Integrate(0 to phi){PDF(phi) dphi}.                              (5)

 Those general equations will be the starting point for development of the
 functions of each source. Without specific information for a real source,
 sources are assumed to be ideal in the sense that S(psi,phi) has no directional
 dependency (is a constant, when compensating for strictly geometrical factors,
 such as cosine weight for plane irradiance sensors as sources) within its
 surface area and emission solid angle (= 0 elsewhere).

 In the current version, all sources as treated as with infinitesimal area
 (point). Sources with finite spatial extent can be implemented in the future, 
 but add complexity as angular distribution is dependent on position within the
 surface of the source (e.g. a hemispherical scalar irradiance sensor as a
 source has restricted range of emission angles everywhere except on the top, 
 since hemispherical irradiance sensors have a circular guard on the bottom to
 prevent light rays from the bottom hemisphere).

 IMPLEMENTATION:
 The implementation of these functions follow the paradigm of dynamic
 dispatching, as used throughout the program. This flexibility might result in
 penalty performances due to additional function pointer deferencing and unknown
 functions at compile time preventing inlining. However, limited execution time
 tests showed no loss in performance against direct call of specific functions.

 The group of functions for each source are:
 src_qtl_x_y - Emission quantile function.
 src_cdf_x_y - Emission cumulative distribution function;
 src_pdf_x_y - Emission probability distribution function;
 src_ef_x_y  - Emission function;

 where x can be:
 pnt - Point source (infinitesimal area);
 srf - Area source (finite surface area);

 and y can be:
 sun - The sun; 					// NOT IMPLEMENTED YET
 las - A laser with Gaussian profile; 			// NOT IMPLEMENTED YET
 sir - A scalar irradiance sensor as a source;
 rad - A radiance sensor as a source;
 pir - A plane irradiance sensor as a source;

 An isotropic point source can be implemented as a "pnt_sir". A radiance sensor 
 as a point source can also be implemented with "pnt_sir". The Sun is treated as 
 a point source only and lasers are obligatory area sources.

 It is possible to sample a random direction by using the qtl function with
 random uniform (0,1] deviates for probabilities of psi and phi.

*******************************************************************************/

 #include <stdlib.h>
 #include <string.h>
 #include <math.h>

 #include "config.h"
 #include "memory.h"
 #include "constants.h"
 #include "geometry.h"
 #include "rotation.h"
 #include "sources.h"

/* src_alloc: ******************************************************************

 Allocates memory to a source struct.

 INPUT:
 None.

 OUTPUT:
 Pointer to a source struct.

*******************************************************************************/

 struct source*
 src_alloc
 ( void )
 {
   struct source* src = (struct source*) malloc(sizeof(struct source));
   #ifdef CHCK_MEM_ALLOC
   if ( !src )
   {
     printf("\nERROR: Failed to allocate memory for src in src_alloc "
       "(source.c)\n");
     exit(-1);
   }
   #endif // CHCK_MEM_ALLOC
 
   src->M = calloc_2d(3, 3, "src->M in src_alloc (sources.c)");
   return src;
 }

/* src_free: *******************************************************************

 Deallocate memory of a source struct.

 INPUT:
 src - Source parameters
       Pointer to pointer to source struct;

 OUTPUT:
 None. Releases memory and updates pointed address to NULL.

*******************************************************************************/

 void
 src_free
 ( struct source **src )
 {
   free_2d(3, &(*src)->M);
   free(*src);
   *src = NULL;
 }

/* src_setup: ******************************************************************

 Setups the parameters of the source.

 INPUT:
 src_fov - Field-of-view of the source;
 src_s   - Pointer to spherical directions of the source axis;
 src_o   - Pointer to center position (origin) of the source;
 src_tp  - Pointer to string with source type name;
 src_par - Pointer to source parameter pointer.

 OUTPUT:
 None. Allocate and initialize source parameters.

*******************************************************************************/

 void
 src_setup
 (
   struct source *src,
   double const src_fov,
   double const *src_s,
   double const *src_ref_o,
   double const *src_rel_o,
   double const *src_stks,
   char const   *src_tp
 )
 {
   src->fov     = src_fov;
   src->hfov    = src_fov * 0.5;
   src->lnf     = (1.0 - cos(src->hfov));
   src->lnf_inv = 1.0 / src->lnf;

   strncpy(src->tp, src_tp, STRMXLEN);
   if ( strcmp (src_tp, "pir") == 0 )
   {
     src->stks[0] = M_1_PI;
     src->fun.qtl = &src_qtl_pnt_pir; 
     src->fun.cdf = &src_cdf_pnt_pir; 
     src->fun.pdf = &src_pdf_pnt_pir; 
     src->fun.ef  = &src_ef_pnt_pir; 
   }
   else if ( strcmp(src_tp, "sir") == 0 )
   {
     src->stks[0] = 1.0 / (K_2PI * src->lnf);
     src->fun.qtl = &src_qtl_pnt_sir; 
     src->fun.cdf = &src_cdf_pnt_sir; 
     src->fun.pdf = &src_pdf_pnt_sir; 
     src->fun.ef  = &src_ef_pnt_sir; 
   }
   else if ( strcmp(src_tp, "rad") == 0 )
   {
     src->stks[0] = 1.0 / (K_2PI * src->lnf);
     src->fun.qtl = &src_qtl_pnt_rad; 
     src->fun.cdf = &src_cdf_pnt_rad; 
     src->fun.pdf = &src_pdf_pnt_rad; 
     src->fun.ef  = &src_ef_pnt_rad; 
   }

   #ifdef VECTOR_RT
   for (size_t i = 0; i < 4; i++)
     src->stks[i] = src_stks[i];
   #else
   src->stks[0] = src_stks[0];
   #endif // VECTOR_RT

   src->s[0] = src_s[0];
   src->s[1] = src_s[1];
   src->s[2] = 1.0;
   sph_to_cos_unit(src->u, src->s);

   if ( (1.0 - ABS(src->u[2]) < TOLERANCE) )
   {
     src->u[0] = 0.0;
     src->u[1] = 0.0;
     src->u[2] = 1.0 * SIGN(src->u[2]);
     src->M[0][0] = 1.0 * SIGN(src->u[2]);
     src->M[1][1] = 1.0 * SIGN(src->u[2]);
     src->M[2][2] = 1.0 * SIGN(src->u[2]);
     src->rotate_f = (src->u[2] > 0.0) ? 0: 1;
   }
   else
   {
     rot_mat_ZY0(src->M, src->u);
     src->rotate_f = 1;
   }

   double src_rel_o_rot[3];   
   rot_vec_unit(src_rel_o_rot, src_rel_o, (double const **) src->M);
   src->o[0] = src_ref_o[0] + src_rel_o_rot[0];
   src->o[1] = src_ref_o[1] + src_rel_o_rot[1];
   src->o[2] = src_ref_o[2] + src_rel_o_rot[2];
   if ( (1.0 / src_ref_o[2]) == -INFINITY  && ABS( src_rel_o_rot[2] ) < TOLERANCE)
     src->o[2] = -0.0; // The sign of -0.0 is loss with -0.0+0.0 or 0.0-0.0.

//mat_fprintf(stdout, (double const **)  src->M, 3, 3, 1.0, "M", 1);
//vec_fprintf(stdout, src_rel_o, 3, 1.0, "REL", 1);
//vec_fprintf(stdout, src_rel_o_rot, 3, 1.0, "ROT_REL", 1);
//vec_fprintf(stdout, src_ref_o, 3, 1.0, "REF", 1);
//vec_fprintf(stdout, src->o, 3, 1.0, "FINAL", 1);
 }

/* src_printf: *****************************************************************

 Prints parameters of the source.

 INPUT:
 odv    - Output device. Should be a file stream (e.g. 'stdout' to print to
          default);
 src    - Source parameters
          Pointer to pointer to constant source struct;
 indent - Number of indentations at the header level (each indent = 2 spaces)
          Constant int;

 OUTPUT:
 None. Prints the parameters of the source.

*******************************************************************************/

 void
 src_fprintf
 ( 
   FILE *odv,
   struct source const *src,
   int const indent
 )
 {
   // Set indentation:
   char pre_0[STRMXLEN] = "";
   char pre_1[STRMXLEN] = "";
   for (size_t i = 0; i < indent; i++)
     strcat(pre_0, "  ");
   strncpy(pre_1, pre_0, STRMXLEN);
   strcat(pre_1, "  ");

   char src_lnm[STRMXLEN];
   if ( strcmp (src->tp, "pir") == 0 )
     strncpy(src_lnm, "Plane irradiance sensor (as source)", STRMXLEN);
   if ( strcmp (src->tp, "sir") == 0 )
     strncpy(src_lnm, "Scalar irradiance sensor (as source)", STRMXLEN);
   if ( strcmp (src->tp, "rad") == 0 )
     strncpy(src_lnm, "Radiance sensor (as source)", STRMXLEN);

   fprintf(odv, "%sSource parameters:\n", pre_0);
   fprintf(odv, "%sType:     %s\n", pre_1, src_lnm);
   fprintf(odv, "%sFOV:      %.2lfº\n", pre_1, src->fov * DEG);   
   fprintf(odv, "%sOrigin:  % .2e, % .2e, % .2e (m)\n", pre_1, 
     src->o[0], 
     src->o[1], 
     src->o[2]);
   fprintf(odv, "%sAxis:      %6.2lfº, %6.2lfº (% .2e, % .2e, % .2e)\n", pre_1, 
     src->s[0] * DEG, 
     src->s[1] * DEG, 
     src->u[0], 
     src->u[1],
     src->u[2]);
   #ifdef VECTOR_RT
   fprintf(odv, "%sStokes:  % .2e, % .2e, % .2e, % .2e (W m^-2 sr^-1)\n", pre_1, 
     src->stks[0], src->stks[1], src->stks[2], src->stks[3]);
   #else
   fprintf(odv, "%sStokes:  % .2e (W m^-2 sr^-1)\n", pre_1, src->stks[0]);
   #endif // VECTOR_RT
   fprintf(odv, "%sRotate:   %d\n", pre_1, src->rotate_f);
   mat_fprintf(odv, (double const **) src->M, 3, 3, 1.0, "  Rotation matrix",
     indent+1);
 }

/* Point sources, scalar irradiance sensors as sources: ************************

 src_qtl_sir, src_cdf_sir, src_pdf_sir, src_ef_sir

 Isotropic point sources (area -> 0, FOV = 2 PI) and scalar irradiance sensors
 as a source.

 The FOV in this case in principle is variable, though typically is PI 
 (hemispherical sensor), or 2 PI (spherical sensor, isotropic point source). It
 is assumed that the sensor is equally sensitive to all directions in its FOV,
 resulting in no directional dependency of the emission. We will take psi_u,
 the integration upper limit, to be FOV / 2, then Eq. A1 becomes:

 Integral{ dphi Integral{ dpsi S(psi,phi) sin(psi) } } = 1 ->
 Integral{ dphi Integral(0 to psi_u){ dpsi S(psi,phi) sin(psi) } } = 1 ->
 2 PI S Integral(0 to psi_u){ dpsi sin(psi)} = 1 ->
 2 PI S (1 - cos(psi_u)) = 1 ->
 S = 1 / (2 PI (1 - cos(psi_u)))                                            (B1)

 Substituting this S into Eq. A2 and A3, the PDF and CDF for psi are then:

 PDF(psi) = S(psi,phi) sin(psi) / Integrate{S(psi,phi) sin(psi) dpsi} ->
          = S sin(psi) / Integrate(0 to psi_u){S sin(psi) dpsi} ->
          = S sin(psi) / S (1 - cos(psi_u)) ->
          = sin(psi) / (1 - cos(psi_u)),                                    (B2)

 CDF(psi) = Integrate(0 to psi){dpsi PDF(psi)} ->
          = Integrate(0 to psi){dpsi sin(psi) / (1 - cos(psi_u))} ->
          = (1 / (1 - cos(psi_u))) Integrate(0 to psi){dpsi sin(psi)} ->
          = (1 - cos(psi)) / (1 - cos(psi_u))                               (B3)

 The quantile can be found by solving Eq. B3 for psi:

 CDF(psi) = (1 - cos(psi)) / (1 - cos(psi_u)) ->
 CDF(psi) (1 - cos(psi_u)) = 1 - cos(psi) ->
 cos(psi) = 1 - CDF(psi) (1 - cos(psi_u)) ->
 psi      = acos(1 - CDF(psi) (1 - cos(psi_u)))                             (B4) 
 
 Similarly, substituting this S into Eq. A4 and A5, the PDF and CDF for phi are
 then:

 PDF(phi) = Integrate{dpsi S(psi,phi) sin(psi)} / CDF(psi,phi) ->
          = Integrate(0 to psi_u){dpsi S sin(psi)} / CDF(psi,phi) ->
          = Integrate(0 to psi_u){dpsi S sin(psi)} / 
              (2 PI Integral(0 to psi_u){dpsi S sin(psi)}) ->
          = 1 / (2 PI)                                                      (B5)

 CDF(phi) = Integrate(0 to phi){dphi PDF(phi)} ->
          = Integrate(0 to phi){dphi / (2 PI)} ->
          = phi / (2 PI)                                                    (B6)

 The quantile can be found by solving Eq. B6 for phi:

 CDF(phi) = phi / (2 PI) ->
 phi      = 2 PI CDF(phi)                                                   (B7) 

 All angles are in the reference frame of the sensor. If the source axis is not
 parallel and with the same direction as the Z-axis, then the emission angle
 needs to be rotated to the system's reference frame. This is not performed by
 the source functions. Random sampling from the PDF is made by providing a
 positive random uniform deviate as CDF(psi) (termed prob, for probability).

*******************************************************************************/

 void
 src_qtl_pnt_sir
 (
   struct source const *src,
   double const *prob,
   double *s_src
 )
 {
   s_src[0] = acos(1.0 - prob[0] * src->lnf);
   s_src[1] = K_2PI * prob[1];
 }

 void
 src_cdf_pnt_sir
 (
   struct source const *src,
   double *prob,
   double const *s_src
 )
 {
   prob[1] = s_src[1] * K_1_2PI;
   if ( s_src[0] < src->hfov )
   {
     prob[0] = (1.0 - cos(s_src[0])) * src->lnf_inv;
   } else {
     prob[0] = 1.0;
   }
 }

 void
 src_pdf_pnt_sir
 (
   struct source const *src,
   double *prob,
   double const *s_src
 )
 {
   if ( s_src[0] < src->hfov )
   {
     prob[0] = sin(s_src[0]) * src->lnf_inv;
     prob[1] = K_1_2PI;
   } else {
     prob[0] = 0.0;
     prob[1] = 0.0;
   }
 }

 void
 src_ef_pnt_sir
 (
   struct source const *src,
   double *stks,
   double const *s_src
 )
 { 
   if ( s_src[0] <= src->hfov )
   {
     #ifdef VECTOR_RT
     for (size_t i = 0; i < 4; i++)
       stks[i] = src->stks[i];
     #else
     stks[0] = src->stks[0];
     #endif // VECTOR_RT
   } else {
     #ifdef VECTOR_RT
     for (size_t i = 0; i < 4; i++)
       stks[i] = 0.0;
     #else
     stks[0] = 0.0;
     #endif // VECTOR_RT
   }
 }

/* Radiance sensors as sources: ************************************************

 src_qtl_rad, src_cdf_rad, src_pdf_rad, src_ef_rad

 For infinitesimal area source, the equations and functions for a radiance
 sensor as a source are exactly equal to those of a scalar irradiance sensor as
 a source.

 However, since in future versions finite area sources might be implemented,
 those functions are repeated here. All equations appearing in this block in the
 future will have prefix 'C'.

*******************************************************************************/

 void
 src_qtl_pnt_rad
 (
   struct source const *src,
   double const *prob,
   double *s_src
 )
 {
   s_src[0] = acos(1.0 - prob[0] * src->lnf);
   s_src[1] = K_2PI * prob[1];
 }

 void
 src_cdf_pnt_rad
 (
   struct source const *src,
   double *prob,
   double const *s_src
 )
 {
   prob[1] = s_src[1] * K_1_2PI;
   if ( s_src[0] < src->hfov )
   {
     prob[0] = (1.0 - cos(s_src[0])) * src->lnf_inv;
   } else {
     prob[0] = 1.0;
   }
 }

 void
 src_pdf_pnt_rad
 (
   struct source const *src,
   double *prob,
   double const *s_src
 )
 {
   if ( s_src[0] < src->hfov )
   {
     prob[0] = sin(s_src[0]) * src->lnf_inv;
     prob[1] = K_1_2PI;
   } else {
     prob[0] = 0.0;
     prob[1] = 0.0;
   }
 }

 void
 src_ef_pnt_rad
 (
   struct source const *src,
   double *stks,
   double const *s_src
 )
 {
   if ( s_src[0] <= src->hfov )
   {
     #ifdef VECTOR_RT
     for (size_t i = 0; i < 4; i++)
       stks[i] = src->stks[i];
     #else
     stks[0] = src->stks[0];
     #endif // VECTOR_RT
   } else {
     #ifdef VECTOR_RT
     for (size_t i = 0; i < 4; i++)
       stks[i] = 0.0;
     #else
     stks[0] = 0.0;
     #endif // VECTOR_RT
   }
 }

/* Plane sources (e.g., irradiance sensors as sources): ************************

 src_qtl_pir, src_cdf_pir, src_pdf_pir, src_ef_pir

 Plane sources, such as plane irradiance sensors as a source.

 In this case, the FOV is fixed to PI as the receiving surface must be flat. Now
 the sensor has a cosine angular response (that is assumed to be well corrected
 to remove cosine response errors). With similar development as for the
 scalar irradiance sensors, we can solve Eq. A1 as:

 Integral{ dphi Integral{ dpsi S(psi,phi) sin(psi) } } = 1 ->
 2 PI S Integral{ dpsi cos(psi) sin(psi)} = 1 ->
 2 PI S 0.5 = 1 ->
 S = 1 / PI,                                                                (D1)
 S(psi) = S cos(psi) = cos(psi) / PI                                        (D2) 

 Substituting this S(psi) into Eq. A2 and A3, the PDF and CDF for psi are then:

 PDF(psi) = S(psi,phi) sin(psi) / Integrate{dpsi S(psi,phi) sin(psi)} ->
          = S cos(psi) sin(psi) / Integrate{dpsi S cos(psi) sin(psi)} ->
          = S cos(psi) sin(psi) / (S 0.5) ->
          = 2 cos(psi) sin(psi),                                            (D3)

 CDF(psi) = Integrate(0 to psi){dpsi PDF(psi)} ->
          = Integrate(0 to psi){dpsi 2 cos(psi) sin(psi)} ->
          = 2 Integrate(0 to psi){dpsi cos(psi) sin(psi)} ->
          = sin(psi)^2
          = 1 - cos(psi)^2                                                  (D4)

 The quantile can be found by solving Eq. D4 for psi:

 CDF(psi) = sin(psi)^2 ->
 sqrt(CDF(psi)) = sin(psi) ->
 psi = asin( sqrt( CDF(psi) ) ) = acos( sqrt( 1 - CDF(psi) ) ).             (D5)

 Similarly, substituting this S into Eq. A4 and A5, the PDF and CDF for phi are
 then:

 PDF(phi) = Integrate{dpsi S(psi,phi) sin(psi)} / CDF(psi,phi) ->
          = Integrate{dpsi S cos(psi) sin(psi)} / CDF(psi,phi) ->
          = Integrate{dpsi S cos(psi) sin(psi)} / 
              (2 PI Integral{dpsi S cos(psi) sin(psi)}) ->
          = 1 / (2 PI)                                                      (D6)

 CDF(phi) = Integrate(0 to phi){dphi PDF(phi)} ->
          = Integrate(0 to phi){dphi / (2 PI)} ->
          = phi / (2 PI)                                                    (D7)

 Those are exactly the same as Eq. B5 and Eq. B6. The quantile can be found by
 solving Eq. B6 for phi:

 CDF(phi) = phi / (2 PI) ->
 phi      = 2 PI CDF(phi)                                                   (D8) 

 All angles are in the reference frame of the sensor. If the source axis is not
 parallel and with the same direction as the Z-axis, then the emission angle
 needs to be rotated to the system's reference frame. This is not performed by
 the source functions. Random sampling from the PDF is made by providing a
 positive random uniform deviate as CDF(psi) (termed prob, for probability).

*******************************************************************************/

 void
 src_qtl_pnt_pir
 (
   struct source const *src,
   double const *prob,
   double *s_src
 )
 {
   s_src[0] = acos(sqrt(1.0 - prob[0]));
   s_src[1] = K_2PI * prob[1];
 }

 void
 src_cdf_pnt_pir
 (
   struct source const *src,
   double *prob,
   double const *s_src
 )
 {
   prob[1] = s_src[1] * K_1_2PI;
   if ( s_src[0] < src->hfov )
   {
     double mu = cos(s_src[0]);
     prob[0] = 1.0 - (mu * mu);
   } else {
     prob[0] = 1.0;
   }
 }

 void
 src_pdf_pnt_pir
 (
   struct source const *src,
   double *prob,
   double const *s_src
 )
 {
   if ( s_src[0] < src->hfov )
   {
     prob[0] = 2.0 * cos(s_src[0]) * sin(s_src[0]);
     prob[1] = K_1_2PI;
   } else {
     prob[0] = 0.0;
     prob[1] = 0.0;
   }
 }

 void
 src_ef_pnt_pir
 (
   struct source const *src,
   double *stks,
   double const *s_src
 )
 {
   if ( s_src[0] <= src->hfov )
   {
     double mu = cos(s_src[0]);
     #ifdef VECTOR_RT
     for (size_t i = 0; i < 4; i++)
       stks[i] = src->stks[i] * mu;
     #else
     stks[0] = src->stks[0] * mu;
     #endif // VECTOR_RT
   } else {
     #ifdef VECTOR_RT
     for (size_t i = 0; i < 4; i++)
       stks[i] = 0.0;
     #else
     stks[0] = 0.0;
     #endif // VECTOR_RT
   }
 }


