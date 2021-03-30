
/*******************************************************************************
 airtight_alex.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 Version: 1.6
 Date: 2021-03-25
 License: GPL-3.0

 This source code provide functions to resolve specific interface interactions 
 for a airtight skylight-blocking cylinder, as used in the PONDER project. It is
 a specific condition but can be generalized for other skylight-blocking
 structures (e.g., cones).

 The code allows for both right-circular and oblique-circular cylinders, with a
 maximum inclination axis of the cylinder to the Z-axis of 40 degrees. It should 
 represent measurements where an airtight structure is used, with low
 inclination angles, when the water level does not rise inside the tube to its
 external level. The interface has the same position and angle of the opening of 
 cylinder.

 Two functions are provided to resolve surface interactions for an emitted ray 
 (downwards, air to water) and a ray in the system reaching the opening of the
 cylinder (upward, water to air).

 This simulation requires the conditional compilation flag AIRTIGHT_ALEX to
 be defined in the config.h header file of passed as a gcc compiler option 
 -DAIRTIGHT_ALEX. Note that because the rays will be reflected downwards when 
 reaching the opening of the cylinder, the result for the FREE condition cannot
 realistically represent the absence of such a sensor since the interface is
 discontinuous. To calculate shadowing is therefore necessary to run the
 simulation twice one with AIRTIGHT_ALEX defined and one without, normalizing
 the shadowing values of the AIRTIGHT_ALEX to the FREE values of the simulation
 without AIRTIGHT_ALEX.

*******************************************************************************/

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

 The pseudo-code is:
 (1) Center ray's initial point to the cylinder's origin and if necessary, 
     rotate point and direction to its reference frame;
 (2) Center ray's initial position to cylinder's top origin and if necessary, 
     rotate point and direction to its reference frame;
 (3) Get intersection point on the supporting plane of the opening (cylinder's
     top);
 (4) Calculate Snell refraction and Fresnel transmittance;
 (5) If incidence angle lower than the critical angle for total internal
     reflection, update the Stokes vector of the reflected ray;
 (6) Back-center (and back-rotate if necessary) from the cylinder's top
     reference frame to the cylinder's reference frame;
 (7) Back-center (and back-rotate if necessary) from the cylinder's reference 
     frame to the system's reference frame;

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
   struct str_cyln const * cyln,
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
       for (int i = 0; i < 3; i++)
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
       for (int i = 0; i < 3; i++)
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
         for (int cw = 0; cw < iop_nw0; cw++)
         {
           for (int cb = 0; cb < btt_nbr; cb++)
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
       for (int i = 0; i < 3; i++)
       {
         a_rot_2[i] = intrs_p[i];
         b_rot_2[i] = a_rot_2[i] + u_rot_2[i] * (l - intrs_l);
       }

       // Rotate back:
       if ( cyln->rotate_f[2] )
       {
         ROT_VEC(a_rot_1, a_rot_2, ref_o_2, M_open, c_vec);
         ROT_VEC(b_rot_1, b_rot_2, ref_o_2, M_open, c_vec);
         ROT_VEC_UNIT(u_rot_1, u_rot_2, M_open);
       }
       else 
       {
         for (int ci = 0; ci < 3; ci++)
         {
           a_rot_1[ci] = a_rot_2[ci] - ref_o_2[ci];
           b_rot_1[ci] = b_rot_2[ci] - ref_o_2[ci];
           u_rot_1[ci] = u_rot_2[ci];
         }
       }

       if ( cyln->rotate_f[0] )
       {
         ROT_VEC(ray->a, a_rot_1, ref_o_1, src->M, c_vec);
         ROT_VEC(ray->b, b_rot_1, ref_o_1, src->M, c_vec);
         ROT_VEC_UNIT(ray->u, u_rot_1, src->M);
       }
       else 
       {
         for (int ci = 0; ci < 3; ci++)
         {
           ray->a[ci] = a_rot_1[ci] - ref_o_1[ci];
           ray->b[ci] = b_rot_1[ci] - ref_o_1[ci];
           ray->u[ci] = u_rot_1[ci];
         }
       }

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
   struct str_cyln const * cyln,
   struct source const * src,
   double const muc, 
   double CMPLX_T const *n_rat_aw,
   double CMPLX_T const iop_na,
   double CMPLX_T const iop_nw,
   struct Fm Fmat,
   int const btt_nbr,
   int const iop_nw0,
   double const * ref_o_1,
   double const * ref_o_2,
   double const ** M_open,
   gsl_rng * random
 )
 {
   double a_rot_1[3], a_rot_2[3], b_rot_1[3];
   double u_rot_1[3], u_rot_2[3];
   double c_vec[3], s_refrac[3];
   double CMPLX_T u_refrac[3];
   double CMPLX_T s_sin;
   double sin_t;
   int intrs_f;
   double intrs_l;
   double intrs_p[3];
   FRESNEL_MC_VAR

   if ( cyln->rotate_f[0] )
   {
     ROT_VEC(a_rot_1, ray->a, cyln->o, cyln->M, c_vec);
     ROT_VEC_UNIT(u_rot_1, ray->u, cyln->M);
   }
   else
   {
     for (int i = 0; i < 3; i++)
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
     for (int i = 0; i < 3; i++)
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

   SNELL_U(u_rot_2[2], u_refrac[2], n_rat_aw, s_sin);
   FRESNEL_MC_U(u_rot_2[2], u_refrac[2], iop_na, iop_nw, Fmat);

   for (int cw = 0; cw < iop_nw0; cw++)
   {
     for (int cb = 0; cb < btt_nbr; cb++)
     {
       #ifdef VECTOR_RT
       printf("\nERROR: VECTOR TRANSMITTANCE NOT IMPLEMENTED\n");
       exit(-1);
       #else
       ray->stks[cw][cb][0] *= Fmat.T[0];
       #endif // VECTOR_RT
     }
   }
   COS_TO_SPH_UNIT(s_refrac, u_rot_2, random);
   s_refrac[0] = CMPLX_F(creal)(CMPLX_P(acos)(u_refrac[2]));
   SPH_TO_COS_UNIT(u_rot_2, s_refrac, sin_t);

   if ( cyln->rotate_f[2] )
   {
     ROT_VEC(b_rot_1, intrs_p, ref_o_2, M_open, c_vec);
     ROT_VEC_UNIT(u_rot_1, u_rot_2, M_open);
   }
   else 
   {
     for (int ci = 0; ci < 3; ci++)
     {
       b_rot_1[ci] = intrs_p[ci] - ref_o_2[ci];
       u_rot_1[ci] = u_rot_2[ci];
     }
   }

   if ( cyln->rotate_f[0] )
   {
     ROT_VEC(ray->b, b_rot_1, ref_o_1, src->M, c_vec);
     ROT_VEC_UNIT(ray->u, u_rot_1, src->M);
   }
   else 
   {
     for (int ci = 0; ci < 3; ci++)
     {
       ray->b[ci] = b_rot_1[ci] - ref_o_1[ci];
       ray->u[ci] = u_rot_1[ci];
     }
   }

   COS_TO_SPH_UNIT(ray->s, ray->u, random);

   if ( !intrs_f )
     ray->nsdw_f = 0;
 }

