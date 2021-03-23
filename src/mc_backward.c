

 #include <stdlib.h>		// calloc, free, exit
 #include <string.h>		// strcmp
 #include <math.h>		// M_PI, M_1_PI, INFINITY, sqrt, asin, sin, log
 #include <complex.h>		// creal, csqrt, cacos
 #include <sys/time.h>     	// timeval, gettimeofday
 #include <gsl/gsl_rng.h>  	// gsl_rng, gsl_rng_alloc, gsl_rng_free, gsl_rng_set, gsl_rng_uniform, gsl_rng_uniform_pos

 #include "config.h"		// SUN_OMG_W, SHADOWING, AIRTIGHT_ALEX, VECTOR_RT 
 #include "aux.h"		// FINDBIN
 #include "constants.h"		// K_2PI
 #include "memory.h"		// calloc_1d, calloc_2d, free_1d, free_2d
 #include "intersect.h"		// chk_intrs
 #include "geometry.h"		// SPH_TO_COS_UNIT, COS_TO_SPH_UNIT, mat_transpose
 #include "ray.h"		// RAY_INI_S, RAY_INI_U, RAY_MOV, struct light_ray, ray_alloc, ray_free
 #include "sources.h"		// SRC_QTL, struct source
 #include "tr_interface.h"	// Fm, SNELL_U, FRESNEL_MC_U, FRESNEL_MC_VAR
 #include "skyrad.h"		// struct skyradiance
 #include "accumulators.h"	// ACCM_B_ADD, accm_b_norm, accm_b_add_k, struct accumulator_bmc
 #include "bottom.h"		// struct bottom
 #include "reflectance.h"	// REFL_BRDF, REFL_DHR, REFL_QTL
 #include "scattering.h"	// SCAT_PF, SCAT_QTL, struct scattering
 #include "rotation.h"		// ROT_VEC_UNIT, ROT_MAT_ZY0
 #include "structures.h"	// str_cyln, str_cone, str_cubd
 #include "airtight_alex.h"	// check_interface_down_airtight_alex, check_interface_up_airtight_alex

 int
 bmc
 (
   int const sim_nr,
   int const sim_ns, 
   int const iop_nw0,
   int const btt_nbr,
   double const * sim_sza, 
   double const * sim_saa, 
   double const sim_f0,
   double CMPLX_T const iop_na, 
   double CMPLX_T const iop_nw, 
   double const iop_c, 
   double const * iop_w0,
   int const str_ncl,
   int const str_ncn,
   int const str_ncb,
   struct source const * src,
   struct scattering const * scat,
   struct bottom const * btt,
   struct skyradiance const * skr,
   struct str_cyln const ** cylns,
   struct str_cone const ** cones,
   str_cubd const ** cubds,
   struct accumulator_bmc * accm_df_f,
   struct accumulator_bmc * accm_df_s,
   struct accumulator_bmc * accm_dr_f,
   struct accumulator_bmc * accm_dr_s
 )
 {
   int    ci  = 0;								// General counting variable
   int    cs  = 0;								// Sun angle counting variable
   int    cr  = 0; 								// Ray simulation number counting variable
   int    cw  = 0;								// Single scattering albedo counting variable
   int    cb  = 0;								// Bottom reflectance counting variable

   double rttaz  = 0.0;								// Rotated azimuth (needed for the current sky radiance definition)
   double sclw   = 0.0;								// Scaling vector for ray power
   double sclwdr = 0.0;								// Temporary scaling for direct component stokes
   double psi = 0.0; 								// Polar scattering angle
   double phi = 0.0;								// Azimuthal scattering angle
   double * dif_n = calloc_1d (sim_ns + 2, "dif_n (mc.c)");			// Pointer to normalization multiplier for diffuse component
   double * dir_n = calloc_1d (sim_ns, "dir_n (mc.c)");				// Pointer to normalization multiplier vector for the direct component

   // Variables for the accumulator:
   int rid, cid;

   // Variables for light ray:
   struct light_ray ray;							// Ray
   struct light_ray ray_vrt; 							// Virtual ray
   ray_alloc(1, iop_nw0, btt_nbr, &ray);
   ray_alloc(1, iop_nw0, btt_nbr, &ray_vrt);
   double l = 0.0;								// Metric pathlength to transverse
   double l_vrt = 0.0;								// Temporary pathlength for surfaces and for sun

   // Variables for light source:
   double inv_src_zpos = 1.0 / src->o[2];
   double s_em[3] = {1.0};							// Spherical directions of emitted ray
   double u_rot[3];

   // Variables for scattering:
   double pfv[1];								// Scattering phase function value at an given angle
   double s_scat[3] = {1.0};							// Spherical direction of scattering
   double u_scat[3];								// Cartesian direction of scattering
   double ** M_scat = calloc_2d(3, 3, "M_scat (mc_backward.c)");		// Rotation matrix from the ray's reference frame to the system's reference frame

   // Variables for interaction with the interface:
   struct Fm Fmat; 								// Fresnel reflectance / transmittance
   double sin_t;
   double CMPLX_T ra_mu = 0.0;							// Cosine of the refracted angle polar angle
   double CMPLX_T s_sin;							// Intermediate calculation when using macro SNELL_U
   double s_refrac[3];								// Refracted angle in spherical directions.
   double u_refrac[3];								// Refracted angle in Cartesian directions.
   double CMPLX_T n_rat_wa[3];							// Refractive index ratio water to air
   double CMPLX_T n_rat_aw[3];							// Refractive index ratio air to water
   n_rat_wa[0] = iop_nw / iop_na;						// Ratio of water to air refractive indexes
   n_rat_wa[1] = n_rat_wa[0] * n_rat_wa[0];					// n_rat_wa^2
   n_rat_wa[2] = n_rat_wa[1] * n_rat_wa[0];					// n_rat_wa^3
   n_rat_aw[0] = 1.0 / n_rat_wa[0];						// Ratio of air to water refractive indexes
   n_rat_aw[1] = 1.0 / n_rat_wa[1];						// n_rat_aw^2
   n_rat_aw[2] = 1.0 / n_rat_wa[2];						// n_rat_aw^3
   double muc = CMPLX_F(creal)(CMPLX_P(sqrt)(1.0 - n_rat_aw[1]));		// Critical angle for total internal reflection
   FRESNEL_MC_VAR								// Group of variables needed to use the macro FRESNEL_MC_U

   // Variables for random sampling:
   double prob[2];
   struct timeval tv;
   gettimeofday (&tv, 0);
   unsigned long int seed = tv.tv_sec + tv.tv_usec;
   gsl_rng * random = gsl_rng_alloc (gsl_rng_ranlxs0);
   gsl_rng_set (random, seed);

   // Variables for interaction with the bottom:
   double s_scat_i[3] = {1.0};							// Incident spherical direction
   double s_scat_r[3] = {1.0};							// Reflected spherical direction
   double * btt_dhr = calloc_1d(btt_nbr, "btt_dhr (mc_backward.c)");		// Directional-hemispherical reflectance of the bottom at a given incidence angle
   double * btt_bdr = calloc_1d(btt_nbr, "btt_bdr (mc_backward.c)");		// Bi-directional reflectance of the bottom for a given geometry

   // Variables for interaction from bottom to Sun: 
   struct Fm * sun_Fmat = (struct Fm*) calloc(sim_ns, sizeof(struct Fm));	// Fresnel transmittance / reflectance matrix
   double ** sun_sw = calloc_2d (sim_ns, 3, "sun_sw (mc.c)");			// Underwater Sun nadir angle, with dimensions [sim_ns][2]: [][0] Polar angle; [][1] Azimuth
   double ** sun_sa = calloc_2d (sim_ns, 3, "sun_sa (mc.c)");			// In-air Sun nadir angle, with dimensions [sim_ns][2]: [][0] Polar angle; [][1] Azimuth
   double ** sun_uw = calloc_2d (sim_ns, 3, "sun_uw (mc.c)");			// Underwater Sun Cartesian directions, with dimensions [sim_ns][3]: [][0] X; [][1] Y; [][2] Z
   double ** sun_ua = calloc_2d (sim_ns, 3, "sun_ua (mc.c)");			// In-air Sun Cartesian directions, with dimensions [sim_ns][3]: [][0] X; [][1] Y; [][2] Z
   double * sun_rra = calloc_1d (sim_ns, "sun_rra (mc.c)");			// Sun relative azimuth rotation for the sky radiance distribution reference, with dimension [sim_ns]
   double * sun_btt_l   = calloc_1d (sim_ns, "sun_btt_l (mc.c)");		// Pathlengths from bottom to surface in the direction of Sun, with dimension [sim_ns]
   double * sun_btt_scl = calloc_1d (sim_ns, "sun_btt_scl (mc.c)");		// Partial scaling for reflectance at bottom to Sun, with dimension [sim_ns]

   for (cs = 0; cs < sim_ns; cs++) 
   {
     sun_sa[cs][0] = M_PI - sim_sza[cs];
     sun_sa[cs][1] = sim_saa[cs];
     SPH_TO_COS_UNIT(sun_ua[cs], sun_sa[cs], sin_t)

     SNELL_U(cos(sim_sza[cs]), ra_mu, n_rat_aw, s_sin);
     sun_sw[cs][0] = M_PI - CMPLX_F(creal)(CMPLX_P(acos)(ra_mu));
     sun_sw[cs][1] = sim_saa[cs];

     // M_PI is used since the sky rad used here has the sun at M_PI...
     sun_rra[cs] = M_PI - sim_saa[cs];
     SPH_TO_COS_UNIT(sun_uw[cs], sun_sw[cs], sin_t);

     FRESNEL_MC_U(cos(M_PI - sun_sw[cs][0]), cos(sim_sza[cs]), iop_na, iop_nw,
       sun_Fmat[cs]);

     sun_btt_l[cs]   = btt->depth / -sun_uw[cs][2];
     sun_btt_scl[cs] = SUN_OMG_W * exp(-sun_btt_l[cs] * iop_c) * 
       sun_Fmat[cs].T[0];
   }

   #ifdef AIRTIGHT_ALEX
   // Variables necessary to propagate from source to exit of airtight tube:
   double ** M_open = calloc_2d(3, 3, "M_open (mc_backward.c)");		// Reverse rotation from opening axis to cylinder/source axis
   mat_transpose(M_open, (double const **) cylns[0]->top->M, 3, 3);		// M_open is the inverse (= transpose for rotation matrices) of the rotation from cylinder to top
   double ref_o_1[3];								// System origin from the reference frame of the cylinder
   double ref_o_2[3];								// Cylinder origin from the perspective of the reference frame of the top
   for (ci = 0; ci < 3; ci++)
   {
     ref_o_2[ci] = -cylns[0]->h * cylns[0]->top->M[ci][2];			
     ref_o_1[ci] = -cylns[0]->o[2] * cylns[0]->M[ci][2];			
   }
   #endif // AIRTIGHT_ALEX

   // Main loop of simulation: *************************************************

   for(cr = 0; cr < sim_nr; cr++) 
   {
     /*
     Start a new ray. The pseudo-code is:
     (1) Sample the emission direction from the source's PDF;
     (2) Initiate ray at the source position and sampled emission direction;
     (3) If the source axis does not have the same direction as the system's 
         Z-azis, rotate the emission direction from the source's reference frame
         to the system's reference frame.

     Note that because the source's origin is described on the system's
     reference frame, it is not necessary to rotate the initial point of the
     ray.
     */ 
     prob[0] = gsl_rng_uniform (random);
     prob[1] = gsl_rng_uniform (random);
     SRC_QTL(src, prob, s_em);
     RAY_INI_S(src->o, s_em, src->stks, ray, sin_t)
     if ( src->rotate_f )
     {
       ROT_VEC_UNIT(u_rot, ray.u, src->M);
       for (ci = 0; ci < 3; ci++)
         ray.u[ci] = u_rot[ci];
       COS_TO_SPH_UNIT(ray.s, ray.u, random);
     }

     /*
     If the source is above the air-water interface, propagate the ray to and
     across surface. The pseudo-code is:
     (1) Calculate the pathlength in air to the interface and move the ray
         terminal point to the interface;
     (2) Calculate the angle of refraction and the Fresnel transmittance from 
         air to water;
     (3) Update the ray's direction to the refracted direction;
     (4) Update the ray's Stokes vector for the transmittance across the
         interface.

     Note that sources above water must have negative depths. -0.0 can be used 
     to indicate the the sensor acting as source is just above the surface. 
     However, since -0.0 cannot be checked directly for its sign, the inverse
     position is checked since the inverse is guaranteed to keep the sign (e.g., 
     1.0 / 0.0 = INFINITY and 1.0 / -0.0 = -INFINITY).

     Note also that the refractive index change also changes the solid angle and
     alters the magnitude of the diffuse Stokes vector. However, in this code, 
     we do not track the magnitude of the Stokes vector, but its relative value 
     to the square of the refractive index: L / n^2, where n is the refractive 
     index of the current medium. Therefore, the transmittance is simply 
     1.0 - Reflectance. At the end of the simulation, refractive index effects
     on the solid angle are taken into account.
     */
     #ifndef AIRTIGHT_ALEX
     if ( inv_src_zpos < 0.0 )
     {
       l = -(src->o[2]) / ray.u[2];
       RAY_MOV(l, ray);
       ray.b[2] = round(ray.b[2] * 1E6) / 1E6;
       SNELL_U(ray.u[2], u_refrac[2], n_rat_aw, s_sin);
       FRESNEL_MC_U(ray.u[2], u_refrac[2], iop_na, iop_nw, Fmat);
       ray.s[0] = CMPLX_F(creal)(CMPLX_P(acos)(u_refrac[2]));
       SPH_TO_COS_UNIT(ray.u, ray.s, sin_t);
       for (cw = 0; cw < iop_nw0; cw++)
       {
         for (cb = 0; cb < btt_nbr; cb++)
         {
           #ifdef VECTOR_RT
           printf("\nERROR: VECTOR TRANSMITTANCE NOT IMPLEMENTED\n");
           exit(-1);
           #else
           ray.stks[cw][cb][0] *= Fmat.T[0];
           #endif // VECTOR_RT
         }
       }
     }
     #else
     /*
     If AIRTIGHT_ALEX flag is defined, the sensor acting as source is inside an
     airtight tube represented by the first cylinder and ray has to come from 
     the source to the opening of the tube, traveling in air, to intersect the
     air-water interface. Due to being airtight, having the small radius, and
     small maximum allowed inclination, the water does not go up in the tube and
     the interface is at the same depth and with the same inclination as the 
     opening (top of the cylinder).
   
     The steps are presented in detail in the airtight_alex.c source file.
     */
     check_interface_down_airtight_alex(&ray, INFINITY, cylns[0], src, muc, 
       n_rat_aw, iop_na, iop_nw, Fmat, btt_nbr, iop_nw0, ref_o_1, ref_o_2,
       (double const **) M_open, random);
     #endif // AIRTIGHT_ALEX

     /*
     Loop through scattering events:

     The limit is defined by the Stokes intensity magnitude. When the intensity
     of the ray is below the critical weight (CRTCW), the ray is extinguished. 
     This test looks at the ray for the highest single scattering albedo and 
     highest bi-hemispherical reflectance of the bottom.
     */
     while ( ray.stks[iop_nw0 - 1][btt_nbr - 1][0] > CRTCW )
     {
       // Sample free optical path length and move ray:
       l = -log(gsl_rng_uniform_pos (random)) / iop_c;
       RAY_MOV(l, ray);

       #ifdef AIRTIGHT_ALEX
       /*
       For Alex's airtight system, check (if in path up) if the ray reached
       the air-water interface on the tip and reflect it.

       The steps are presented in detail in the airtight_alex.c source file.
       */
       check_interface_up_airtight_alex(&ray, INFINITY, cylns[0], src, muc, 
         n_rat_wa, iop_na, iop_nw, Fmat, btt_nbr, iop_nw0, ref_o_1, ref_o_2,
         (double const **) M_open, random);
       #endif // AIRTIGHT_ALEX

       /*
       Loop over multiple interaction between surface and bottom:

       Such interaction could occur when water are very transparent, bottom at a
       very shallow depth and/or bottom highly reflective.
       */
       while ( ( (ray.b[2] <= 0) || (ray.b[2] > btt->depth) ) && 
               (ray.stks[iop_nw0 - 1][btt_nbr - 1][0] > CRTCW) )
       {
         // If ray interacts with surface:
         if ( ray.b[2] <= 0.0 )
         {
           /*
           Set virtual ray to be at the position of the intersection with the
           surface. To avoid rounding errors, the Z-axis position is set
           directly.
           */
           l_vrt = ray.b[2] / ray.u[2];
           ray_vrt.b[0] =  ray.b[0] + (-ray.u[0] * l_vrt);
           ray_vrt.b[1] =  ray.b[1] + (-ray.u[1] * l_vrt);
           ray_vrt.b[2] =  0.0; // to avoid rounding errors

           #ifdef SHADOWING
           /*
           Check for intersection with shadowing structures in the path to 
           surface.
           */
           if (ray.nsdw_f)
           {
             ray.nsdw_f = !chk_intrs (ray.a, ray_vrt.b, ray.u, l-l_vrt, 
               str_ncl, cylns, str_ncn, cones, str_ncb, cubds);
           }
           #endif // SHADOWING

           // Update ray after surface interaction. Note that l_vrt is the 
           // remainder path of the total pathlength.
           l = l_vrt;
           ray.a[0] = ray_vrt.b[0];
           ray.a[1] = ray_vrt.b[1];
           ray.a[2] = ray_vrt.b[2];
           ray.b[2] = -ray.b[2];
           ray.u[2] = -ray.u[2];
           ray.s[0] = M_PI - ray.s[0];
           ray.s[2] = 1.0;

           /*
           If incidence angle is lower than the critical angle for total 
           internal reflection, calculate partial transmittance of the ray.
           Note: The sign of ray.u[2] was inverted before, and it is now
           pointing downwards, so ray.u[2] > muc is the correct test (the 
           cosine will be larger if the angle is smaller).
           */
           if ( ray.u[2] > muc )
           {
             // Calculate refraction angle and Fresnel reflection:
             // ra_mu = snell_u(ray.u[2], n_rat_wa);
             SNELL_U(ray.u[2], u_refrac[2], n_rat_wa, s_sin);
             FRESNEL_MC_U(ray.u[2], u_refrac[2], iop_na, iop_nw, Fmat);

             /*
             Update backward diffuse accumulator for the unshadowed condition:
             Note that the current azimuthal rotation of ray.s[1] + sun_rra[cs]
             is done since sun azimuth is arbitrarily defined while sky rad is
             fixed.
             */
             rttaz = ray.s[1] + sun_rra[0];
             rttaz = (rttaz > K_2PI) ? rttaz - K_2PI : rttaz;
             s_refrac[0] = acos(u_refrac[2]);
             s_refrac[1] = rttaz;
             ACCM_B_ADD(accm_df_f, skr, ray_vrt.b, s_refrac, ray.stks, 
               Fmat.T[0], 0, 0);

             #ifdef SHADOWING
             /*
             If ray under shadowing condition is still active, update the
             diffuse component accumulator for the shadowed condition if no
             intersection occurs in air.
             */
             if ( ray.nsdw_f ) 
             {
               u_refrac[0] = sin(M_PI - s_refrac[0]) * cos(ray.s[1]); 
               u_refrac[1] = sin(M_PI - s_refrac[0]) * sin(ray.s[1]);
               u_refrac[2] = cos(M_PI - s_refrac[0]);
               if ( !chk_intrs (ray_vrt.b, ray_vrt.b, u_refrac, 1000.0,  
                      str_ncl, cylns, str_ncn, cones, str_ncb, cubds) )
               {
                 ACCM_B_ADD(accm_df_s, skr, ray_vrt.b, s_refrac, ray.stks, 
                   Fmat.T[0], 0, 0);
               }
             }
             #endif // SHADOWING

             // Update the stokes intensities of reflected ray:
             for (cw = 0; cw < iop_nw0; cw++)
             {
               for (cb = 0; cb < btt_nbr; cb++)
               {
                 ray.stks[cw][cb][0] *= Fmat.R[0]; //fr;
               } // Loop: btt_nbr
             } // Loop: iop_nw0
           } // Logical: Does ray refract?
         } // Logical: Is ray in air?

         /*
         Interaction of rays with the bottom:

         The interaction with the bottom require that at each reflection the
         potential contribution of the ray in the direction of the Sun be 
         accounted for, and this take full account of possible multiple 
         reflections that can happen between bottom and surface when in clear 
         and shallow waters.
         
         Note that for simplicity, we enter only the loop if depth is higher 
         than bottom, not exactly equal. The reason is that a sampled reflected
         direction from the bottom BRDF can be 90 degrees, such that ray will
         stay at the bottom surface and the loop will not break (additionally
         some division per 0 might occur as Z-axis Cartesian direction will be 
         0.

         The pseudo-code for bottom interaction is as follows:
         ( 1) Calculate path back to depth equal to bottom depth;
         ( 2) Move ray terminal point backwards to bottom depth;
         ( 3) SHADOWING: Did ray intersect a structure in its way downwards?
         ( 4) Get DHR for incident angle on bottom;
         ( 5) Loop over all Sun angles:
             (5.1) Get BDR for incident angle and reflection angle to Sun;
             (5.2) Propagate virtual ray from bottom to surface;
             (5.3) Update Stokes vector of virtual ray;
             (5.4) Update accumulator for FREE condition;
             (5.5) SHADOWING: If no intersection occurs in path upwards and in
                  air, update accumulator for SHADOWING condition;
         ( 6) Sample reflected direction from BRDF;
         ( 7) Update Stokes vector of the ray;
         ( 8) Update ray initial point to be at the bottom;
         ( 9) Move ray terminal point by its remainder path;
         (10) AIRTIGHT_ALEX: Test if ray reaches the interface at the opening of
              the source in its way upwards;

         It is not necessary to test the AIRTIGHT_ALEX condition for the virtual
         ray from bottom to surface, as AIRTIGHT_ALEX produces only SHADOWING
         results and if virtual ray would reach the interface at the opening of
         the source, it would also intersect the source structure and its
         contribution will not be accounted for.
         */
         if ( ray.b[2] > btt->depth )
         {
           /*
           Calculate path back to bottom and set current position to be at the 
           bottom. Since in the 'if' above we restrict the loop if positions are 
           higher than bottom depth, the current code below is sufficient. 
           Otherwise, some of the rays will reflect at 90 degrees so will have 
           their depth still equal to the bottom and the u[2] will be zero, 
           resulting in NaN.
           */
           l_vrt = (ray.b[2] - btt->depth) / ray.u[2];
           ray.b[0] += (-ray.u[0] * l_vrt);
           ray.b[1] += (-ray.u[1] * l_vrt);
           ray.b[2]  = btt->depth;

           #ifdef SHADOWING
           /*
           Check for intersection with shadowing structures in the path to 
           bottom:
           */
           if (ray.nsdw_f)
           {
             ray.nsdw_f = !chk_intrs (ray.a, ray.b, ray.u, l-l_vrt, str_ncl, 
               cylns, str_ncn, cones, str_ncb, cubds);
           }
           #endif //  SHADOWING

           /*
           Update accumulators for the direct component:
           
           Note that this has to be done after the update of the ray Stokes 
           vector for interaction with the bottom since the sun weighting does 
           not include the bottom bi-hemisphercal reflectance - note however 
           that this might be an issue in the future if polarized BRDFs are 
           required, as the Stokes vector will require the final direction 
           before being updated.
           */
           s_scat_i[0] = ray.u[2];
           s_scat_i[1] = ray.s[1];
           s_scat_i[2] = 1.0;
           REFL_DHR(&btt->refl, s_scat_i, btt_dhr);

           for (cs = 0; cs < sim_ns; cs++)
           {
             s_scat_r[0] = sun_uw[cs][2];
             s_scat_r[1] = sun_sw[cs][1];
             s_scat_r[2] = 1.0;
             REFL_BRDF(&btt->refl, s_scat_i, s_scat_r, btt_bdr);

             /*
             Propagate virtual ray to surface (Z == 0):
             To avoid rounding errors in Z, it is directly set to 0.0.
             */
             ray_vrt.b[0] = ray.b[0] + (sun_uw[cs][0] * sun_btt_l[cs]);
             ray_vrt.b[1] = ray.b[1] + (sun_uw[cs][1] * sun_btt_l[cs]);
             ray_vrt.b[2] = 0.0;

             /*
             Calculate Stokes intensity of virtual ray accounting for the 
             probability of exiting to Sun. Note that in this code, the BRDF
             is not normalized, that is, it includes both the probability of
             scattering to a given angle and of ray survival.
             */
             for (cw = 0; cw < iop_nw0; cw++)
             {
               for (cb = 0; cb < btt_nbr; cb++)
               {
                 ray_vrt.stks[cw][cb][0] = ray.stks[cw][cb][0] * 
                   btt_bdr[cb] * sun_btt_scl[cs];
               }
             }

             // Update accumulator for unshadowed condition:
             s_refrac[0] = sim_sza[cs];
             s_refrac[1] = sim_saa[cs];
             ACCM_B_ADD(accm_dr_f, skr, ray_vrt.b, s_refrac, ray_vrt.stks, 
               1.0, 1, cs);

             #ifdef SHADOWING
             /*
             If ray under shadowing condition is still active, update its 
             accumulator for the direct component.
             */
             if ( ray.nsdw_f )
             {
               /*
               If no intersection occurs in its path to surface, and no 
               intersection occurs in air, add its contribution:
               */
               if ( !chk_intrs (ray.b, ray_vrt.b, sun_uw[cs], sun_btt_l[cs], 
                      str_ncl, cylns, str_ncn, cones, str_ncb, cubds) )
               {
                 if ( !chk_intrs (ray_vrt.b, ray_vrt.b, sun_ua[cs], 1000.0, 
                        str_ncl, cylns, str_ncn, cones, str_ncb, cubds) )
                 {
                   ACCM_B_ADD(accm_dr_s, skr, ray_vrt.b, s_refrac, 
                     ray_vrt.stks, 1.0, 1, cs);
                 } // Logical: chk_intrs (air)
               }  // Logical: chk_intrs (water)
             } // Logical: Is ray active for shadow?
             #endif // SHADOWING
           } // Loop: sim_ns

           // Update stokes intensity after bottom interaction:
           for (cw = 0; cw < iop_nw0; cw++)
           {
             for (cb = 0; cb < btt_nbr; cb++)
             {
               ray.stks[cw][cb][0] *= btt_dhr[cb];
             }
           }

           // Get reflected direction after bottom interaction:
           prob[0] = gsl_rng_uniform(random);
           prob[1] = gsl_rng_uniform(random);
           s_scat_i[0] = ray.u[2];
           s_scat_i[1] = ray.s[1];
           REFL_QTL(&btt->refl, s_scat_i, s_scat_r, prob);

           l = l_vrt;
           ray.s[0] = acos(-s_scat_r[0]);
           ray.s[1] = s_scat_r[1];
           ray.s[2] = 1.0;
           SPH_TO_COS_UNIT(ray.u, ray.s, sin_t);

           // Update ray's position:
           RAY_MOV(l, ray);

           #ifdef AIRTIGHT_ALEX
           /*
           For Alex's airtight system, check if in the path up it reached 
           air-water interface on the tip.
           */
           check_interface_up_airtight_alex(&ray, INFINITY, cylns[0], src, muc, 
             n_rat_wa, iop_na, iop_nw, Fmat, btt_nbr, iop_nw0, ref_o_1,
             ref_o_2, (double const **) M_open, random);
           #endif // AIRTIGHT_ALEX
         } // Logical: Is ray below the bottom?
       } // Loop: Is ray in air or below the bottom?

       #ifdef SHADOWING
       /*
       If the ray terminal point is not above the water or below the bottom,
       the surface-bottom loop was not entered or has just finished and is 
       necessary to check for intersection with shadowing structures.
       */
       if ( ray.nsdw_f )
       {
         ray.nsdw_f = !chk_intrs (ray.a, ray.b, ray.u, l, str_ncl, cylns, 
           str_ncn, cones, str_ncb, cubds);
       }
       #endif // SHADOWING

       /*
       Update accumulators for the direct component:

       Here we calculate the probability of scattering to the necessary 
       direction, and to be transmitted to and across the surface. This is 
       similar to what is done for the case of reflection at the bottom, 
       however, the angular probability here is given by the phase function 
       instead of the BRDF.

       Looping over all Sun angles, the pseudo-code is:
       (1) Get scattering angle from current direction to the Sun direction;
       (2) Get PF value for the scattering angle;
       (3) Get pathlength to the surface;
       (4) Propagate virtual ray to the surface;
       (5) Update the Stokes vector of the virtual ray;
       (6) Update accumulator for FREE condition;
       (7) SHADOWING: If the virtual ray does not intersect structures in its 
           path upwards and in air, update accumulator for SHADOWING condition;

       Since here we are dealing with a "virtual" ray generated out of the real
       ray for the local estimate method, the final fate of the ray is not
       relevant. Therefore it is not necessary to include the tests for Alex's 
       airtight system. This is because is a ray would intersect the interface
       at the opening of the cylinder, it would necesserally intersect the
       cylinder and therefore not contribute to the direct component
       accumulator. No specular reflection is necessary as this is a virtual
       and is no longer tracked.

       Because in the multiple interaction between bottom and surface the ray
       Stokes intensity can get lower than CRTCW, the surface-bottom loop will
       break with the ray in an arbitrary position below the bottom or above the
       water. The ray should be extinguished but it is not possible to skip to
       next ray at this point, so we only proceed here if the Stokes intensity
       is still higher than CRTCW and the ray will be extinguished after
       scattering, with no consequence for accumulators.
       */
       if ( ray.stks[iop_nw0 - 1][btt_nbr - 1][0] > CRTCW )
       {
         for (cs = 0; cs < sim_ns; cs++)
         {
           /*
           Get scattering angle from current direction to Sun direction, the 
           scattering phase function value for that scattering angle and the
           pathlength to the surface.
           */
           s_scat[0] = acos(ray.u[2] * sun_uw[cs][2] + 
                       sqrt(1.0 - ray.u[2] * ray.u[2]) * 
                       sqrt(1.0 - sun_uw[cs][2] * sun_uw[cs][2]) * 
                       cos(ray.s[1]));
           s_scat[2] = 1.0;
           l_vrt = ray.b[2] / -sun_uw[cs][2];
           SCAT_PF(scat, s_scat, pfv);

           /*
           Propagate virtual ray to surface (Z == 0). To avoid rounding errors 
           in Z, it is directly set to 0.0.
           */
           ray_vrt.b[0] = ray.b[0] + (sun_uw[cs][0] * l_vrt);
           ray_vrt.b[1] = ray.b[1] + (sun_uw[cs][1] * l_vrt);
           ray_vrt.b[2] = 0.0;

           // Update accumulator for unshadowed condition:
           sclwdr = SUN_OMG_W * exp(-l_vrt * iop_c) * pfv[0] * 
             sun_Fmat[cs].T[0];
           for (cw = 0; cw < iop_nw0; cw++)
           {
             for (cb = 0; cb < btt_nbr; cb++)
             {
               ray_vrt.stks[cw][cb][0] = ray.stks[cw][cb][0] * iop_w0[cw] * 
                 sclwdr;
             }
           }

           s_refrac[0] = sim_sza[cs];
           s_refrac[1] = sim_saa[cs];
           ACCM_B_ADD(accm_dr_f, skr, ray_vrt.b, s_refrac, ray_vrt.stks, 1.0,
             1, cs);

           #ifdef SHADOWING
           /*
           If no intersection occurs in its path to surface, and no intersection
           occurs in air, add its contribution.
           */
           if ( ray.nsdw_f )
           {
             if ( !chk_intrs (ray.b, ray_vrt.b, sun_uw[cs], l_vrt, str_ncl, 
                    cylns, str_ncn, cones, str_ncb, cubds) )
             {
               if ( !chk_intrs (ray_vrt.b, ray_vrt.b, sun_ua[cs], 1000.0, 
                      str_ncl, cylns, str_ncn, cones, str_ncb, cubds) )
               {
                 ACCM_B_ADD(accm_dr_s, skr, ray_vrt.b, s_refrac, ray_vrt.stks,
                   1.0, 1, cs);
               }
             }
           }
           #endif // SHADOWING
         }
       }

       // Update stokes intensity of the photon package after interaction with
       // matter and sample new scattered direction:
       for (cw = 0; cw < iop_nw0; cw++)
       {
         for (cb = 0; cb < btt_nbr; cb++)
         {
           ray.stks[cw][cb][0] *= iop_w0[cw];
         }
       }
       prob[0] = gsl_rng_uniform(random);
       prob[1] = gsl_rng_uniform(random);

       SCAT_QTL(scat, s_scat, prob);
       ROT_MAT_ZY0(M_scat, ray.u, sin_t);
       s_scat[2] = 1.0;
       SPH_TO_COS_UNIT(u_scat, s_scat, sin_t);
       RAY_SCAT(s_scat, u_scat, M_scat, sin_t, ray);
     }
   }

   // Normalize the results: ***************************************************
   /*
   The simulation is completed and we set the normalization of the accumulators.

   The normalization is set here for underwater detectors acting as sources. If 
   source is above the water (negative Z), the iop_nw^2/iop_na^2 factor is
   removed by multiplying the normalization factors by iop_na^2 / iop_nw^2.
   */

   // Set normalization for radiance sensors acting as sources:
   if ( strcmp (src->tp, "rad") == 0 )
   {
     for (cs = 0; cs < (sim_ns + 2); cs++)
     {
       dif_n[cs] = n_rat_wa[1] / sim_nr;
       if ( cs < sim_ns )
       {
         dir_n[cs] = 1.0 / (sim_nr * -sun_uw[cs][2] * SUN_OMG_W);
       }
     }
   }

   // Set normalization for scalar irradiance sensors acting as sources:
   if ( strcmp (src->tp, "sir") == 0 )
   {
     dif_n[0] = K_2PI * n_rat_wa[1] / sim_nr;
     dir_n[0] = K_2PI / (sim_nr * -sun_uw[cs][2] * SUN_OMG_W);
     for (cs = 1; cs < (sim_ns + 2); cs++)
     {
       dif_n[cs] = dif_n[0];
       if ( cs < sim_ns )
         dir_n[cs] = dir_n[0];
     }
   }

   // Set normalization for plane irradiance sensors acting as sources:
   if ( strcmp (src->tp, "pir") == 0 )
   {
     for (cs = 0; cs < (sim_ns + 2); cs++) 
     {
       dif_n[cs] = M_PI * n_rat_wa[1] / sim_nr;
       if (cs < sim_ns)
       {
         dir_n[cs] = M_PI / (sim_nr * -sun_uw[cs][2] * SUN_OMG_W);
       }
     }
   }

   // If necessary, remove the iop_nw^2/iop_na^2 factor:
   if ( (src->o[2] < 0) || ((1.0 / src->o[2]) == -INFINITY) )
   {
     for(cs = 0; cs < (sim_ns + 2); cs++)
     {
       dif_n[cs] *= n_rat_aw[1];
       if ( cs < sim_ns )
       {
         dir_n[cs] *= n_rat_aw[1];
       }
     }
   }

   // Normalize results:
   accm_b_norm(accm_df_f, sim_f0, dif_n);
   accm_b_norm(accm_dr_f, sim_f0, dir_n);
   #ifdef SHADOWING
   accm_b_norm(accm_df_s, sim_f0, dif_n);
   accm_b_norm(accm_dr_s, sim_f0, dir_n);
   #endif // SHADOWING

   /*
   For the direct sensor, it is necessary to add the direct transmission from
   the sensor to the Sun if the Sun is in the FOV of the sensor.
   
   For this purpose, we first find if the direction to the Sun is in the FOV of
   the sensor. If so, we find the probability of transmission to and across the
   surface. Note that the relevant direction is just that towards the Sun, such
   that the transmittance provides all normalization necessary. This value is
   finally scaled by the refractive index change if necessary.
   
   Finally, it is possible that arbitrary structures will block that direction
   and so if structures are present a test is first made to evaluate if the path
   to Sun is free for the accumulator under shadowing.
   */

   for (cs = 0; cs < sim_ns; cs++)
   {
     s_scat[0] = acos(src->u[0] * sun_uw[cs][0] +
                 src->u[1] * sun_uw[cs][1] +
                 src->u[2] * sun_uw[cs][2]);
     s_scat[2] = 1.0;
     if ( s_scat[0] < src->hfov )
     {
       l    = src->o[2] / -sun_uw[cs][2];
       sclw = exp(-l * iop_c) * sun_Fmat[cs].T[0];
       if ( inv_src_zpos )
       {
         sclw *= n_rat_aw[1];
       }
       RAY_INI_U(src->o, sun_uw[cs], src->stks, ray, random);
       RAY_MOV(l, ray);

       FINDBIN(ray.b[1], accm_dr_f->miny, accm_dr_f->ky_inv, accm_dr_f->ny, 
         rid);
       FINDBIN(ray.b[0], accm_dr_f->minx, accm_dr_f->kx_inv, accm_dr_f->nx,
         cid);
       accm_b_add_k(accm_dr_f, sclw, cs, rid, cid);

       #ifdef SHADOWING
       if ( !chk_intrs (ray.a, ray.b, sun_uw[cs], l, str_ncl, cylns, str_ncn, 
              cones, str_ncb, cubds) )
       {
         if ( !chk_intrs (ray.b, ray.b, sun_ua[cs], 1000.0, str_ncl, cylns, 
                str_ncn, cones, str_ncb, cubds) )
         {
           accm_b_add_k(accm_dr_s, sclw, cs, rid, cid);
         }
       }
       #endif // SHADOWING
     }
   }

   // Deallocate memory: *******************************************************
   gsl_rng_free(random);
   free_2d(sim_ns, &sun_sw);
   free_2d(sim_ns, &sun_sa);
   free_2d(sim_ns, &sun_uw);
   free_2d(sim_ns, &sun_ua);
   free(sun_Fmat);
   free_1d(&sun_rra);
   free_1d(&sun_btt_l);
   free_1d(&sun_btt_scl);
   free_2d(3, &M_scat);
   free_1d(&btt_dhr);
   free_1d(&btt_bdr);
   free_1d(&dif_n);
   free_1d(&dir_n);
   ray_free(&ray);
   ray_free(&ray_vrt);
   #ifdef AIRTIGHT_ALEX
   free_2d(3, &M_open);
   #endif // AIRTIGHT_ALEX

   return 0;
 }

