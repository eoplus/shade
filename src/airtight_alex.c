// SPECIFIC ALEX:
//
// For out air-tight system, the air-water boundary is at the exit of the 
// cylinder and therefore the photon package should be propagated directly 
// there. Comment out the next line when not simulating this sensor and when 
// simulating the reference condition without the sensor. Uncomment the line 
// above for teh general condition.

// SPECIFIC ALEX:
//
// For out air-tight system, the air-water boundary is at the exit of the 
// cylinder and therefore it is important to check this surface interaction.
// Comment out when not simulating this sensor and when simulating the reference
// condition without the sensor.
/*

     */

 #include <math.h>

 #include "config.h"
 #include "constants.h"
 #include "ray.h"
 #include "structures.h"
 #include "sources.h"
 #include "tr_interface.h"
 #include "geometry.h"
 #include "rotation.h"
 #include "intersect.h"

 #include <gsl/gsl_rng.h>

/* Perform reflection at the opening: ******************************************

 check_interface_up_airtight_alex

 INPUT:
 ray      - Light ray parameters
 l        - Pathlength to transverse (= Inf since is in air and air RTE is not 
            tracked in version 1.6 of the code
 cyln     - Cylinder that contains the sensor lens (acting as source)
 src      - Source parameters, from which we use the rotation matrix from source
            axis (= cylinder axis) to the system's reference frame.
 muc      - Cosine of the critical angle for total internal reflection (NaN for
            air to water)
 n_rat_aw - Ratio of air to water refractive indexes
            [0]rat, [1]rat^2, [2]rat^3
            Pointer to double;
 iop_na   - Refractive index of air
            Double;
 iop_na   - Refractive index of water
            Double;
 Fmat     - Fresnel reflectance and transmittance matrix
            Struct Fm;
 btt_nbr  - Number of bottom reflectances;
 iop_nw0  - Number of single scattering albedos;
 ref_o_1  - Origin of the system from the perspective of the cylinder's origin;
 ref_o_2  - Origin of the cylinder from the perspective of the cylinder's top;
 M_open   - Back-rotation matrix from the cylinder's reference frame to the
            cylinder's reference frame
 random   - GSL random number generator;

 OUPUT:
 None. Updates the position and direction of the ray.

*******************************************************************************/

 void
 check_interface_up_airtight_alex
 (
   struct light_ray * ray,
   double const l,
   str_cyln const * cyln,
   struct source const * src,
   double const muc, 
   double const *n_rat_wa,
   double const iop_na,
   double const iop_nw,
   struct Fm Fmat,
   int const btt_nbr,
   int const iop_nw0,
   double const * ref_o_1,
   double const * ref_o_2,
   double const ** M_open,
   gsl_rng * random
 )
 {
   double a_rot_1[3], a_rot_2[3];
   double b_rot_1[3], b_rot_2[3];
   double u_rot_1[3], u_rot_2[3];
   double c_vec[3], u_refrac[3];
   double s_sin;
   int intrs_f;
   double intrs_l;
   double intrs_p[3];
   FRESNEL_MC_VAR

   if ( ray->u[2] < 0.0 )
   {
     /*
     Center to cylinder origin and rotate if necessary
     */
     if ( cyln->rotate_f[0] )
     {
       ROT_VEC(a_rot_1, ray->a, cyln->o, cyln->M, c_vec);
       ROT_VEC_UNIT(u_rot_1, ray->u, cyln->M);
     }
     else
     {
       for (size_t i = 0; i < 3; i++)
       {
         a_rot_1[i] = ray->a[i] - cyln->o[i];
         u_rot_1[i] = ray->u[i];
       }
     }

     /*
     Center to cylinder top and rotate if necessary
     */
     a_rot_1[2] -= cyln->h;
     if ( cyln->rotate_f[2] )
     {
       ROT_VEC_UNIT(a_rot_2, a_rot_1, cyln->top->M);
       ROT_VEC_UNIT(u_rot_2, u_rot_1, cyln->top->M);
     }
     else
     {
       for (size_t i = 0; i < 3; i++)
       {
         a_rot_2[i] = a_rot_1[i];
         u_rot_2[i] = u_rot_1[i];
       }
     }

     // Check intersection with plane and the ellipse
     intrs_plane(&intrs_f, &intrs_l, intrs_p, a_rot_2, u_rot_2, l);
     if ( intrs_f )
     {
       intrs_f = ( (intrs_p[0] * intrs_p[0]) * cyln->top->rsq_inv[0] + 
                   (intrs_p[1] * intrs_p[1]) * cyln->top->rsq_inv[1] - 1.0 ) < 
                 TOLERANCE;
     }

     if ( intrs_f )
     {
       // Specular reflection
       u_rot_2[2] *= -1.0;

       // Account for partial transmittance loss
       if ( u_rot_2[2] > muc )
       {
         SNELL_U(u_rot_2[2], u_refrac[2], n_rat_wa, s_sin);
         FRESNEL_MC_U(u_rot_2[2], u_refrac[2], iop_nw, iop_na, Fmat);
         for (size_t cw = 0; cw < iop_nw0; cw++)
         {
           for (size_t cb = 0; cb < btt_nbr; cb++)
           {
             #ifdef VECTOR_RT
             printf("\nERROR: VECTOR TRANSMITTANCE NOT IMPLEMENTED\n");
             exit(-1);
             #else
             ray->stks[cw][cb][0] *= Fmat.R[0];
             #endif // VECTOR_RT
           }
         }
       }

       // Update ray's position:
       a_rot_2[0] = intrs_p[0];
       a_rot_2[1] = intrs_p[1];
       a_rot_2[2] = intrs_p[2];
       b_rot_2[0] = a_rot_2[0] + u_rot_2[0] * (l - intrs_l);
       b_rot_2[1] = a_rot_2[1] + u_rot_2[1] * (l - intrs_l);
       b_rot_2[2] = a_rot_2[2] + u_rot_2[2] * (l - intrs_l);

       // Rotate back:
       ROT_VEC(a_rot_1, a_rot_2, ref_o_2, M_open, c_vec);
       ROT_VEC(b_rot_1, b_rot_2, ref_o_2, M_open, c_vec);
       ROT_VEC_UNIT(u_rot_1, u_rot_2, M_open);

       ROT_VEC(ray->a, a_rot_1, ref_o_1, src->M, c_vec);
       ROT_VEC(ray->b, b_rot_1, ref_o_1, src->M, c_vec);
       ROT_VEC_UNIT(ray->u, u_rot_1, src->M);

       COS_TO_SPH_UNIT(ray->s, ray->u, random);
     }
   }
 }

/* Perform transmission to water: **********************************************

 check_interface_down_airtight_alex

 For this very specific simulation condition, it is understood that the source
 has the same inclination axis as the first cylinder, such that the ray does not
 need to be projected if it was just created. 

 If the top opening of the cylinder is inclined relative to the cylinder's axis,
 the ray position and direction are rotated to the top's reference frame. With 
 this position and direction the point of intersection with the water surface is
 checked and the ray propagated to that position.

 Note that in the code below a test is made for intersection to the ellipse
 formed by the opening, that is the water surface. If the ray does not intersect
 this surface, it intersected the side walls and is flagged.

 The pseudo-code is:
 (1) Center ray's initial point to the cylinder's origin and if necessary, 
     rotate point and direction to its reference frame;
 (2) Center ray's initial position to cylinder's top origin and if necessary, 
     rotate point and direction to its reference frame;
 (3) Get intersection point on the supporting plane of the opening (cylinder's
     top);
 (4) Calculate Snell refraction and Fresnel transmittance;
 (5) Update the Stokes vector of the transmitted ray;
 (6) Back-center (and back-rotate if necessary) from the cylinder's top
     reference frame to the cylinder's reference frame;
 (7) Back-center (and back-rotate if necessary) from the cylinder's reference 
     frame to the system's reference frame;
 (8) If intersection point with the interface was not within the ellipse of the 
     opening (top of the cylinder), ray necessarily intersected the side walls
     and is flagged.

 INPUT:
 ray      - Light ray parameters
 l        - Pathlength to transverse (= Inf since is in air and air RTE is not 
            tracked in version 1.6 of the code
 cyln     - Cylinder that contains the sensor lens (acting as source)
 src      - Source parameters, from which we use the rotation matrix from source
            axis (= cylinder axis) to the system's reference frame.
 muc      - Cosine of the critical angle for total internal reflection (NaN for
            air to water)
 n_rat_aw - Ratio of air to water refractive indexes
            [0]rat, [1]rat^2, [2]rat^3
            Pointer to double;
 iop_na   - Refractive index of air
            Double;
 iop_na   - Refractive index of water
            Double;
 Fmat     - Fresnel reflectance and transmittance matrix
            Struct Fm;
 btt_nbr  - Number of bottom reflectances;
 iop_nw0  - Number of single scattering albedos;
 ref_o_1  - Origin of the system from the perspective of the cylinder's origin;
 ref_o_2  - Origin of the cylinder from the perspective of the cylinder's top;
 M_open   - Back-rotation matrix from the cylinder's reference frame to the
            cylinder's reference frame
 random   - GSL random number generator;

 OUPUT:
 None. Updates the position, direction and shadow flag of the ray.

*******************************************************************************/

 void
 check_interface_down_airtight_alex
 (
   struct light_ray * ray,
   double const l,
   str_cyln const * cyln,
   struct source const * src,
   double const muc, 
   double const *n_rat_aw,
   double const iop_na,
   double const iop_nw,
   struct Fm Fmat,
   int const btt_nbr,
   int const iop_nw0,
   double const * ref_o_1,
   double const * ref_o_2,
   double const ** M_open,
   gsl_rng * random
 )
 {
   double a_rot_1[3], a_rot_2[3];
   double b_rot_1[3], b_rot_2[3];
   double u_rot[3], v_rot[3];
   double c_vec[3], u_refrac[3], s_refrac[3];
   double s_sin, sin_t;
   int intrs_f;
   double intrs_l;
   double intrs_p[3];
   FRESNEL_MC_VAR

   if ( cyln->rotate_f[2] )
   {
     ROT_VEC(u_rot, ray->a, cyln->o, cyln->M, c_vec); // NOW ALSO ROTATE FROM SYSTEM TO CYLINDER FIRST!
     u_rot[2] -= cyln->h;
     ROT_VEC_UNIT(v_rot, u_rot, cyln->top->M);
     ROT_VEC_UNIT(u_rot, ray->u, cyln->top->M);
   }
   else
   {
     for (size_t ci = 0; ci < 3; ci++)
       u_rot[ci] = ray->u[ci];
     ROT_VEC(v_rot, ray->a, cyln->top->o, cyln->M, c_vec);
   }

   intrs_plane(&intrs_f, &intrs_l, intrs_p, v_rot, u_rot, INFINITY);
   intrs_f = ( (intrs_p[0] * intrs_p[0]) * cyln->top->rsq_inv[0] + 
               (intrs_p[1] * intrs_p[1]) * cyln->top->rsq_inv[1] - 1.0 ) < 
             TOLERANCE;
   SNELL_U(u_rot[2], u_refrac[2], n_rat_aw, s_sin);
   FRESNEL_MC_U(u_rot[2], u_refrac[2], iop_na, iop_nw, Fmat);

   for (size_t cw = 0; cw < iop_nw0; cw++)
   {
     for (size_t cb = 0; cb < btt_nbr; cb++)
     {
       #ifdef VECTOR_RT
       printf("\nERROR: VECTOR TRANSMITTANCE NOT IMPLEMENTED\n");
       exit(-1);
       #else
       ray->stks[cw][cb][0] *= Fmat.T[0];
       #endif // VECTOR_RT
     }
   }
   s_refrac[0] = CMPLX_F(creal)(CMPLX_P(acos)(u_refrac[2]));
   s_refrac[1] = ray->s[1];
   SPH_TO_COS_UNIT(u_rot, s_refrac, sin_t);

   if ( cyln->rotate_f[2] )
   {
     ROT_VEC(v_rot, intrs_p, ref_o_2, M_open, c_vec);
     ROT_VEC(ray->b, v_rot, ref_o_1, src->M, c_vec);
     ROT_VEC_UNIT(v_rot, u_rot, M_open);
     ROT_VEC_UNIT(ray->u, v_rot, src->M);
   }
   else 
   {
     for (size_t ci = 0; ci < 3; ci++)
       ray->b[ci] = intrs_p[ci] + cyln->top->o[ci];
     ROT_VEC_UNIT(ray->u, u_rot, src->M);
   }
   COS_TO_SPH_UNIT(ray->s, ray->u, random);

   if ( !intrs_f )
     ray->nsdw_f = 0;
 }

