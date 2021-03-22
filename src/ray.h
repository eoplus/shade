
 #ifndef RAY
 #define RAY

 #include "config.h"
 #include "constants.h"
 #include "aux.h"
 #include "geometry.h"
 #include "rotation.h"

/*******************************************************************************
 light_ray

 This structure defines the light ray's properties. Several informations are 
 stored to facilitate calculations: the initial and terminal points of the ray 
 path, its spherical and Cartesian directions.

 The light properties are stored as a Stokes vector {I,Q,U,V}. Polarization is 
 not tracked in the current version, but will be added in the future. Since all 
 optical properties that do not affect directionallity can be easily included in
 a simulation, the Stokes parameters are stores as a three-dimensional array 
 with all combinations of number of single scattering albedos (w0) and bottom 
 bi-hemispherical reflectances (BHR).

 Two flags indicate properties of the ray: (1) If its active for the detectors
 under shadowed condition; and (2) If it is a virtual ray, such as used for the
 local estimate method. Virtual ray carries information only of directions and 
 positions but not on the Stokes parameters.

 COMPONENTS:
 a      - An array containing the initial point of the ray, meters, [0] X,  
          [1] Y, [2] Z;
 b      - An array containing the terminal point of the ray, meters, [0] X, 
          [1] Y, [2] Z;
 u      - An array containing the Cartesian directions of the ray, unitless, 
          [0] X, [1] Y, [2] Z;
 s      - An array containing the spherical directions of the ray, radians, 
          [0] Theta, [1] Phi, [2] Radius (always 1.0);
 nsdw_f - Flag to indicate if ray is active for shadowing detectors. 
          0 = inactive, 1 = active;
 real_f - Flag to indicate if ray is a virtual ray. 0 = virtual, 1 = real;
 nw0    - Number of single scattering albedos;
 nbr    - Number of bi-hemispherical bottom reflectances;
 stks   - Pointer to an array of pointer arrays containing the Stokes parameters
          for all w0 and BHR combinations, dimensions are: [nw0][nbr][4], and 
          the Stokes components in the last dimension are: [0] I, [1] Q, [2] U, 
          [3] V;

*******************************************************************************/

 struct light_ray
 {
   double a[3];
   double b[3];
   double u[3];
   double s[3];
   int    nsdw_f;
   int    real_f;
   size_t nw0;
   size_t nbr;
   double ***stks;
 };

/* Initiate ray: ***************************************************************

 ray_ini_s, ray_ini_u

 Initiate a light_ray based on initial position and spherical direction (s) or
 Cartesian direction (u).

 INPUT:
 p   - Pointer to array containing the position of emission, meters, X Y Z;
 s   - Pointer to array containing the spherical directions, radians, Theta Phi;
 ray - Pointer to a light_ray.

 OUTPUT:
 None. Update values of the light_ray components.

*******************************************************************************/

 static inline __attribute__((always_inline)) void
 ray_ini_s
 (
   double const * o,
   double const * s,
   double const * stks,
   struct light_ray * ray
 )
 {
   for (size_t i = 0; i < 3; i++)
   {
     ray->a[i] = o[i];
     ray->b[i] = o[i];
     ray->s[i] = s[i];
   }
   sph_to_cos_unit(ray->u, ray->s);

   if ( ray->real_f )
   {
     for (size_t cw = 0; cw < ray->nw0; cw++)
     {
       for (size_t cb = 0; cb < ray->nbr; cb++)
       {
         for (size_t ci = 0; ci < STKS_N; ci++)
         {
           ray->stks[cw][cb][ci] = stks[ci];
         }
       }
     }
   }
   ray->nsdw_f = 1;
 }

 #define RAY_INI_S(I_o, I_s, I_skts, I_ray, I_sin_t) \
         for (size_t i = 0; i < 3; i++) \
         { \
           (I_ray).a[i] = (I_o)[i]; \
           (I_ray).b[i] = (I_o)[i]; \
           (I_ray).s[i] = (I_s)[i]; \
         } \
         SPH_TO_COS_UNIT( (I_ray).u, (I_ray).s, (I_sin_t) ); \
         if ( (I_ray).real_f ) \
         { \
           for (size_t cw = 0; cw < (I_ray).nw0; cw++) \
           { \
             for (size_t cb = 0; cb < (I_ray).nbr; cb++) \
             { \
               for (size_t ci = 0; ci < STKS_N; ci++) \
               { \
                 (I_ray).stks[cw][cb][ci] = (I_skts)[ci]; \
               } \
             } \
           } \
         } \
         (I_ray).nsdw_f = 1;

 static inline __attribute__((always_inline)) void
 ray_ini_u
 (
   double const * o,
   double const * u,
   double const * stks,
   gsl_rng * random,
   struct light_ray * ray
 )
 {
   for (size_t i = 0; i < 3; i++)
   {
     ray->a[i] = o[i];
     ray->b[i] = o[i];
     ray->u[i] = u[i];
   }
   cos_to_sph_unit(ray->s, ray->u, random);

   if ( ray->real_f )
   {
     for (size_t cw = 0; cw < ray->nw0; cw++)
     {
       for (size_t cb = 0; cb < ray->nbr; cb++)
       {
         for (size_t ci = 0; ci < STKS_N; ci++)
         {
           ray->stks[cw][cb][ci] = stks[ci];
         }
       }
     }
   }
   ray->nsdw_f = 1;
 }

 #define RAY_INI_U(I_o, I_u, I_skts, I_ray, I_random) \
         for (size_t i = 0; i < 3; i++) \
         { \
           (I_ray).a[i] = (I_o)[i]; \
           (I_ray).b[i] = (I_o)[i]; \
           (I_ray).u[i] = (I_u)[i]; \
         } \
         COS_TO_SPH_UNIT( (I_ray).s, (I_ray).u, (I_random) ); \
         if ( (I_ray).real_f ) \
         { \
           for (size_t cw = 0; cw < (I_ray).nw0; cw++) \
           { \
             for (size_t cb = 0; cb < (I_ray).nbr; cb++) \
             { \
               for (size_t ci = 0; ci < STKS_N; ci++) \
               { \
                 (I_ray).stks[cw][cb][ci] = (I_skts)[ci]; \
               } \
             } \
           } \
         } \
         (I_ray).nsdw_f = 1;

/* Move ray: *******************************************************************

 ray_mov

 Move a ray along its direction.

 INPUT:
 l   - Pathlength to travel, meters.
 ray - Pointer to a light_ray.

 OUTPUT:
 None. Update values of the light_ray components.

*******************************************************************************/

 static inline __attribute__((always_inline)) void
 ray_mov 
 (
   double l, 
   struct light_ray * ray
 )
 {
   for (size_t i = 0; i < 3; i++)
   {
     ray->a[i]  = ray->b[i];
     ray->b[i] += ray->u[i] * l;
   }
 }

 #define RAY_MOV(I_l, I_ray) \
         for (size_t i = 0; i < 3; i++) \
         { \
           (I_ray).a[i]  = (I_ray).b[i]; \
           (I_ray).b[i] += (I_ray).u[i] * (I_l); \
         }

/* Update ray direction after scattering: **************************************

 ray_scat

 Update the ray's direction after a scattering event. In the current version
 (1.6) polariziation is not fully implemented and scattering uses only the 
 first element of the Muller matrix.

 INPUT:
 s_scat - Spehrical directions of scatter in the scattering reference frame.
 ray    - Pointer to a light_ray.
 random - Pointer to the GSL random number generator (potentially necessary for 
          the conversion from Cartesian to spherical directions). 

 OUTPUT:
 None. Update values of the the u and s members of the light_ray.

*******************************************************************************/

 static inline __attribute__((always_inline)) void
 ray_scat 
 (
   double const * s_scat,
   double * u_scat,
   double ** M_scat, 
   gsl_rng * random,
   struct light_ray *ray
 )
 {
   sph_to_cos_unit(u_scat, s_scat);

   if( (1.0 - ABS(ray->u[2])) < TOLERANCE )
   {
     /*
     If |cos (theta)| is ~ 1, the ray is traveling in a direction that is 
     parallel to the Z-axis, either the same or opposite direction. In this
     case we keep the scattering directions (reflected if necessary, in the 
     case that the ray was travelling -Z-axis). 
     */
     for (size_t i = 0; i < 3; i++)
     {
       ray->u[i] = SIGN(ray->u[2]) * u_scat[i];     
     }
   } 
   else
   {
     /*
     Otherwise, we compute the rotation matrix to bring the scattering direction
     from the scattering reference frame to the system's reference frame.
     */
     #ifdef VECTOR_RTE
     printf("\nERROR: VECTOR SCATTERING NOT IMPLEMENTED");
     exit(-1);
     #else
     rot_mat_ZY0(M_scat, ray->u);
     mat_vec_prod(ray->u, u_scat, (double const **) M_scat, 3, 3);
     #endif // VECTOR_RTE
   }

   /* Update the spherical directions */
   cos_to_sph_unit(ray->s, ray->u, random);   
 }

 #define RAY_SCAT(I_s_scat, I_u_scat, I_M_scat, I_sin_t, I_ray) \
         SPH_TO_COS_UNIT( (I_u_scat), (I_s_scat), (I_sin_t) ); \
         if( (1.0 - ABS( (I_ray).u[2] )) < TOLERANCE ) \
         { \
           for (size_t i = 0; i < 3; i++) \
           { \
             (I_ray).u[i] = SIGN( (I_ray).u[2] ) * (I_u_scat)[i]; \
           } \
         } \
         else \
         { \
           rot_mat_ZY0( (I_M_scat), (I_ray).u ); \
           MAT_VEC_PROD( (I_ray).u, (I_u_scat), (I_M_scat), 3, 3 ); \
        } \
        COS_TO_SPH_UNIT( (I_ray).s, (I_ray).u, (I_random));  

/* Function prototypes: *******************************************************/

 void
 ray_alloc
 (
   int const real_f,
   int const iop_nw0,
   int const btt_nbr,
   struct light_ray * ray
 );

 void
 ray_free
 ( struct light_ray * ray );

 void
 ray_fprintf
 (
   FILE *odv,
   struct light_ray *ray,
   int const indent 
 );

 #endif // RAY

