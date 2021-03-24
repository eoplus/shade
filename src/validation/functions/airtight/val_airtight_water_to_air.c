
/*******************************************************************************
 val_airtight_water_to_air.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 This is a driver program to test the interface interaction at the opening of an
 airtight tube.

 Compile instruction:
 gcc val_airtight_water_to_air.c ../../../geometry.c ../../../structures.c \
  ../../../structures_cylinder.c ../../../rotation.c ../../../intersect.c \
  ../../../intersect_cylinder.c ../../../intersect_cone.c \
  ../../../intersect_cuboid.c ../../../accumulators.c ../../../ray.c \
  ../../../memory.c ../../../aux.c ../../../skyrad.c ../../../sources.c \
  ../../../structures_cone.c ../../../tr_interface.c ../../../airtight_alex.c \
  -DSHADOWING -DAIRTIGHT_ALEX -lm -lgsl -O3 -o val_airtight_water_to_air.o

 ./val_airtight_water_to_air.o >> res_airtight_water_to_air.txt

*******************************************************************************/

 #include <stdio.h>
 #include <math.h>
 #include <gsl/gsl_rng.h>
 #include <sys/time.h>

 #include "../../../config.h"
 #include "../../../aux.h"
 #include "../../../memory.h"
 #include "../../../constants.h"
 #include "../../../geometry.h"
 #include "../../../structures.h"
 #include "../../../rotation.h"
 #include "../../../intersect.h"
 #include "../../../skyrad.h"
 #include "../../../accumulators.h"
 #include "../../../ray.h"
 #include "../../../sources.h"
 #include "../../../tr_interface.h"
 #include "../../../airtight_alex.h"

 void
 main
 ( void )
 {
   int ci;
   FILE * fi;
   long int fpos;
   double origin[3];
   double axis[2];
   double alpha;
   double radius;
   double height;
   double base_s[3];
   double top_s[3];
   int closed[2];

   int iop_nw0 = 1;
   int btt_nbr = 1;
   double iop_nw = 1.34;
   double iop_na = 1.0;

   // Variables for random sampling:
   double prob[2];
   struct timeval tv;
   gettimeofday (&tv, 0);
   unsigned long int seed = tv.tv_sec + tv.tv_usec;
   gsl_rng * random = gsl_rng_alloc (gsl_rng_ranlxs0);
   gsl_rng_set (random, seed);

   // Variables for light ray:
   struct light_ray ray;
   ray_alloc(1, iop_nw0, btt_nbr, &ray);
   double l = 0.0;

   // Variables for light source:
   double inv_src_zpos = 1.0;
   double s_em[3] = {1.0};
   double u_rot[3];

   // Variables for interaction with the interface:
   struct Fm Fmat;
   double sin_t;
   double CMPLX_T ra_mu = 0.0;
   double CMPLX_T s_sin;
   double s_refrac[3];
   double u_refrac[3];
   double CMPLX_T n_rat_wa[3];
   double CMPLX_T n_rat_aw[3];
   n_rat_wa[0] = iop_nw / iop_na;
   n_rat_wa[1] = n_rat_wa[0] * n_rat_wa[0];
   n_rat_wa[2] = n_rat_wa[1] * n_rat_wa[0];
   n_rat_aw[0] = 1.0 / n_rat_wa[0];
   n_rat_aw[1] = 1.0 / n_rat_wa[1];
   n_rat_aw[2] = 1.0 / n_rat_wa[2];
   double muc = CMPLX_F(creal)(CMPLX_P(sqrt)(1.0 - n_rat_aw[1]));
   FRESNEL_MC_VAR

   // Variables necessary to propagate from source to exit of airtight tube:
   double ** M_open = calloc_2d(3, 3, "M_open (val_airtight.c)");
   double ref_o_1[3];
   double ref_o_2[3];

   /* Setup cylinders: *********************************************************

    Two right and two oblique cylinders are used. For each, one is aligned with
    the Z-axis and one is inclined 30 degrees.

    The cylinders are based on the PONDER skylight-blocking structure, with a
    radius of 0.0125 m, a height of 0.11 m. The specific position depends on the 
    inclination, such as to result in 0.02 m below the water level.

   ****************************************************************************/

   int str_ncl = 4;
   struct str_cyln ** cylns;
   cylns = (struct str_cyln **) calloc(str_ncl, sizeof(struct str_cyln *));

   fi = fopen("input_cylinders.txt", "r");
   if(fi == NULL)
   {
     printf("\nERROR: Failed to open input file input_cylinders.txt\n\n");
     exit(-1);
   }
   fpos = ftell(fi);

   for (size_t i = 0; i < str_ncl; i++)
   {
     cylns[i] = str_cyln_alloc();
     str_cyln_read(fi, &fpos, origin, axis, &alpha, &radius, &height, base_s, 
       top_s, closed);
     str_cyln_setup(cylns[i], origin, axis, alpha, radius, height, base_s, 
       top_s, closed);   
   }

   /* Setup sources: ***********************************************************

    All sources are setup to be at 0.085 m from the center of the top opening
    (0.195 m of the base opening, the origin of the cylinder).

   ****************************************************************************/

   char src_tp[] = {"rad"};
   double src_ref_o[3];
   double src_rel_o[3] = {0.0, 0.0, -0.195};
   double src_s[3];
   double src_fov = 7.5 * RAD;
   double src_stks[STKS_N] = {1.0};

   struct source ** srcs;
   srcs = (struct source **) calloc(str_ncl, sizeof(struct source *));
   struct source ** srcs_2;
   srcs_2 = (struct source **) calloc(str_ncl, sizeof(struct source *));

   for (size_t i = 0; i < str_ncl; i++)
   {
     srcs[i] = src_alloc();
     srcs_2[i] = src_alloc();
     for (size_t j = 0; j < 3; j++)
     {
       src_ref_o[j] = cylns[i]->o[j];
       src_s[j] = cylns[i]->s[j];
     }
     src_setup(srcs_2[i], src_fov, src_s, src_ref_o, src_rel_o, src_stks, src_tp);
     src_s[0] = M_PI - src_s[0];
     src_s[1] = M_PI - src_s[1];
     src_setup(srcs[i], src_fov, src_s, src_ref_o, src_rel_o, src_stks, src_tp);
   }

   l = 0.09;
 
   /* Print cylinders and sources: *********************************************

   ****************************************************************************/

   printf("\n");
   str_cylns_fprintf(stdout, (struct str_cyln const **) cylns, str_ncl, 1);
   printf("\n");

   /* From source to the opening of the tube: **********************************

    For each cylinder, a source is setup with the same orientation axis as the 
    cylinder and rays are emitted. The test tracks the directions and positions
    in the cylinder reference frame and rotated to the system's reference frame.

   ****************************************************************************/

   printf("\n\n%s", ANSI_BOLD);
   printf("  Right-circular cylinder (aligned) =============================\n\n%s",
     ANSI_RESET);
   ci = 0;

   str_cyln_fprintf(stdout, cylns[ci], 1);
   printf("    Base opening:\n");
   str_ellp_fprintf(stdout, (str_ellp const *) cylns[ci]->base, 3);
   printf("    Top opening:\n");
   str_ellp_fprintf(stdout, (str_ellp const *) cylns[ci]->top, 3);
   printf("\n");
   src_fprintf(stdout, srcs[ci], 1);

   printf("\n Emission at 0.0º, 0.0º --------------------------------------\n");

   mat_transpose(M_open, (double const **) cylns[ci]->top->M, 3, 3);
   for (size_t i = 0; i < 3; i++)
   {
     ref_o_2[i] = -cylns[ci]->h * cylns[ci]->top->M[i][2];			
     ref_o_1[i] = -cylns[ci]->o[2] * cylns[ci]->M[i][2];			
   }
     // Right-circular cylinder:
     s_em[0] = 0.0;
     s_em[1] = 0.0;
     RAY_INI_S(srcs[ci]->o, s_em, srcs[ci]->stks, ray, sin_t)
     if ( srcs[ci]->rotate_f )
     {
       ROT_VEC_UNIT(u_rot, ray.u, srcs[ci]->M);
       for (size_t i = 0; i < 3; i++)
         ray.u[i] = u_rot[i];
       COS_TO_SPH_UNIT(ray.s, ray.u, random);
     }

     printf("\n  As emitted:\n");
     ray_fprintf(stdout, &ray, 2);
     printf("\n");
     check_interface_up_airtight_alex(&ray, l, cylns[ci], srcs_2[ci],
       muc, n_rat_aw, iop_na, iop_nw, Fmat, btt_nbr, iop_nw0, ref_o_1, ref_o_2,
       (double const **) M_open, random);
     printf("\n  Reflected:\n");
     ray_fprintf(stdout, &ray, 2);
     printf("\n");

   printf("\n Emission at random ------------------------------------------\n");

   mat_transpose(M_open, (double const **) cylns[ci]->top->M, 3, 3);
   for (size_t i = 0; i < 3; i++)
   {
     ref_o_2[i] = -cylns[ci]->h * cylns[ci]->top->M[i][2];			
     ref_o_1[i] = -cylns[ci]->o[2] * cylns[ci]->M[i][2];			
   }
     // Right-circular cylinder:
     prob[0] = gsl_rng_uniform (random);
     prob[1] = gsl_rng_uniform (random);
     SRC_QTL(srcs[ci], prob, s_em);
     RAY_INI_S(srcs[ci]->o, s_em, srcs[ci]->stks, ray, sin_t)
     if ( srcs[ci]->rotate_f )
     {
       ROT_VEC_UNIT(u_rot, ray.u, srcs[ci]->M);
       for (size_t i = 0; i < 3; i++)
         ray.u[i] = u_rot[i];
       COS_TO_SPH_UNIT(ray.s, ray.u, random);
     }

     printf("\n  As emitted:\n");
     ray_fprintf(stdout, &ray, 2);
     printf("\n");
     check_interface_up_airtight_alex(&ray, l, cylns[ci], srcs_2[ci],
       muc, n_rat_aw, iop_na, iop_nw, Fmat, btt_nbr, iop_nw0, ref_o_1, ref_o_2,
       (double const **) M_open, random);
     printf("\n  Reflected:\n");
     ray_fprintf(stdout, &ray, 2);
     printf("\n");






   printf("\n\n%s", ANSI_BOLD);
   printf("  Right-circular cylinder (inclined) ============================\n\n%s",
     ANSI_RESET);
   ci = 1;

   str_cyln_fprintf(stdout, cylns[ci], 1);
   printf("    Base opening:\n");
   str_ellp_fprintf(stdout, (str_ellp const *) cylns[ci]->base, 3);
   printf("    Top opening:\n");
   str_ellp_fprintf(stdout, (str_ellp const *) cylns[ci]->top, 3);
   printf("\n");
   src_fprintf(stdout, srcs[ci], 1);

   printf("\n Emission at 0.0º, 0.0º --------------------------------------\n");

   mat_transpose(M_open, (double const **) cylns[ci]->top->M, 3, 3);
   for (size_t i = 0; i < 3; i++)
   {
     ref_o_2[i] = -cylns[ci]->h * cylns[ci]->top->M[i][2];			
     ref_o_1[i] = -cylns[ci]->o[2] * cylns[ci]->M[i][2];			
   }
     // Right-circular cylinder:
     s_em[0] = 0.0;
     s_em[1] = 0.0;
     RAY_INI_S(srcs[ci]->o, s_em, srcs[ci]->stks, ray, sin_t)
     if ( srcs[ci]->rotate_f )
     {
       ROT_VEC_UNIT(u_rot, ray.u, srcs[ci]->M);
       for (size_t i = 0; i < 3; i++)
         ray.u[i] = u_rot[i];
       COS_TO_SPH_UNIT(ray.s, ray.u, random);
     }

     printf("\n  As emitted:\n");
     ray_fprintf(stdout, &ray, 2);
     printf("\n");
     check_interface_up_airtight_alex(&ray, l, cylns[ci], srcs_2[ci],
       muc, n_rat_aw, iop_na, iop_nw, Fmat, btt_nbr, iop_nw0, ref_o_1, ref_o_2,
       (double const **) M_open, random);
     printf("\n  Reflected:\n");
     ray_fprintf(stdout, &ray, 2);
     printf("\n");

   printf("\n Emission at random ------------------------------------------\n");

   mat_transpose(M_open, (double const **) cylns[ci]->top->M, 3, 3);
   for (size_t i = 0; i < 3; i++)
   {
     ref_o_2[i] = -cylns[ci]->h * cylns[ci]->top->M[i][2];			
     ref_o_1[i] = -cylns[ci]->o[2] * cylns[ci]->M[i][2];			
   }
     // Right-circular cylinder:
     prob[0] = gsl_rng_uniform (random);
     prob[1] = gsl_rng_uniform (random);
     SRC_QTL(srcs[ci], prob, s_em);
     RAY_INI_S(srcs[ci]->o, s_em, srcs[ci]->stks, ray, sin_t)

     if ( srcs[ci]->rotate_f )
     {
       ROT_VEC_UNIT(u_rot, ray.u, srcs[ci]->M);
       for (size_t i = 0; i < 3; i++)
         ray.u[i] = u_rot[i];
       COS_TO_SPH_UNIT(ray.s, ray.u, random);
     }

     printf("\n  As emitted:\n");
     ray_fprintf(stdout, &ray, 2);
     printf("\n");
     check_interface_up_airtight_alex(&ray, l, cylns[ci], srcs_2[ci],
       muc, n_rat_aw, iop_na, iop_nw, Fmat, btt_nbr, iop_nw0, ref_o_1, ref_o_2,
       (double const **) M_open, random);
     printf("\n  Reflected:\n");
     ray_fprintf(stdout, &ray, 2);
     printf("\n");





   printf("\n\n%s", ANSI_BOLD);
   printf("  Oblique-circular cylinder (aligned) ===========================\n\n%s",
     ANSI_RESET);
   ci = 2;

   str_cyln_fprintf(stdout, cylns[ci], 1);
   printf("    Base opening:\n");
   str_ellp_fprintf(stdout, (str_ellp const *) cylns[ci]->base, 3);
   printf("    Top opening:\n");
   str_ellp_fprintf(stdout, (str_ellp const *) cylns[ci]->top, 3);
   printf("\n");
   src_fprintf(stdout, srcs[ci], 1);

   printf("\n Emission at 0.0º, 0.0º --------------------------------------\n");

   mat_transpose(M_open, (double const **) cylns[ci]->top->M, 3, 3);
   for (size_t i = 0; i < 3; i++)
   {
     ref_o_2[i] = -cylns[ci]->h * cylns[ci]->top->M[i][2];			
     ref_o_1[i] = -cylns[ci]->o[2] * cylns[ci]->M[i][2];			
   }
     // Right-circular cylinder:
     s_em[0] = 0.0;
     s_em[1] = 0.0;
     RAY_INI_S(srcs[ci]->o, s_em, srcs[ci]->stks, ray, sin_t)
     if ( srcs[ci]->rotate_f )
     {
       ROT_VEC_UNIT(u_rot, ray.u, srcs[ci]->M);
       for (size_t i = 0; i < 3; i++)
         ray.u[i] = u_rot[i];
       COS_TO_SPH_UNIT(ray.s, ray.u, random);
     }

     printf("\n  As emitted:\n");
     ray_fprintf(stdout, &ray, 2);
     printf("\n");
     check_interface_up_airtight_alex(&ray, l, cylns[ci], srcs_2[ci],
       muc, n_rat_aw, iop_na, iop_nw, Fmat, btt_nbr, iop_nw0, ref_o_1, ref_o_2,
       (double const **) M_open, random);
     printf("\n  Reflected:\n");
     ray_fprintf(stdout, &ray, 2);
     printf("\n");

   printf("\n Emission at random ------------------------------------------\n");

   mat_transpose(M_open, (double const **) cylns[ci]->top->M, 3, 3);
   for (size_t i = 0; i < 3; i++)
   {
     ref_o_2[i] = -cylns[ci]->h * cylns[ci]->top->M[i][2];			
     ref_o_1[i] = -cylns[ci]->o[2] * cylns[ci]->M[i][2];			
   }
     // Right-circular cylinder:
     prob[0] = gsl_rng_uniform (random);
     prob[1] = gsl_rng_uniform (random);
     SRC_QTL(srcs[ci], prob, s_em);
     RAY_INI_S(srcs[ci]->o, s_em, srcs[ci]->stks, ray, sin_t)
     if ( srcs[ci]->rotate_f )
     {
       ROT_VEC_UNIT(u_rot, ray.u, srcs[ci]->M);
       for (size_t i = 0; i < 3; i++)
         ray.u[i] = u_rot[i];
       COS_TO_SPH_UNIT(ray.s, ray.u, random);
     }

     printf("\n  As emitted:\n");
     ray_fprintf(stdout, &ray, 2);
     printf("\n");
     check_interface_up_airtight_alex(&ray, l, cylns[ci], srcs_2[ci],
       muc, n_rat_aw, iop_na, iop_nw, Fmat, btt_nbr, iop_nw0, ref_o_1, ref_o_2,
       (double const **) M_open, random);
     printf("\n  Reflected:\n");
     ray_fprintf(stdout, &ray, 2);
     printf("\n");





   printf("\n\n%s", ANSI_BOLD);
   printf("  Oblique-circular cylinder (inclined) ==========================\n\n%s",
     ANSI_RESET);
   ci = 3;

   str_cyln_fprintf(stdout, cylns[ci], 1);
   printf("    Base opening:\n");
   str_ellp_fprintf(stdout, (str_ellp const *) cylns[ci]->base, 3);
   printf("    Top opening:\n");
   str_ellp_fprintf(stdout, (str_ellp const *) cylns[ci]->top, 3);
   printf("\n");
   src_fprintf(stdout, srcs[ci], 1);

   printf("\n Emission at 0.0º, 0.0º --------------------------------------\n");

   mat_transpose(M_open, (double const **) cylns[ci]->top->M, 3, 3);
   for (size_t i = 0; i < 3; i++)
   {
     ref_o_2[i] = -cylns[ci]->h * cylns[ci]->top->M[i][2];			
     ref_o_1[i] = -cylns[ci]->o[2] * cylns[ci]->M[i][2];			
   }
     // Right-circular cylinder:
     s_em[0] = 0.0;
     s_em[1] = 0.0;
     RAY_INI_S(srcs[ci]->o, s_em, srcs[ci]->stks, ray, sin_t)
     if ( srcs[ci]->rotate_f )
     {
       ROT_VEC_UNIT(u_rot, ray.u, srcs[ci]->M);
       for (size_t i = 0; i < 3; i++)
         ray.u[i] = u_rot[i];
       COS_TO_SPH_UNIT(ray.s, ray.u, random);
     }

     printf("\n  As emitted:\n");
     ray_fprintf(stdout, &ray, 2);
     printf("\n");
     check_interface_up_airtight_alex(&ray, l, cylns[ci], srcs_2[ci],
       muc, n_rat_aw, iop_na, iop_nw, Fmat, btt_nbr, iop_nw0, ref_o_1, ref_o_2,
       (double const **) M_open, random);
     printf("\n  Reflected:\n");
     ray_fprintf(stdout, &ray, 2);
     printf("\n");

   printf("\n Emission at random ------------------------------------------\n");

   mat_transpose(M_open, (double const **) cylns[ci]->top->M, 3, 3);
   for (size_t i = 0; i < 3; i++)
   {
     ref_o_2[i] = -cylns[ci]->h * cylns[ci]->top->M[i][2];			
     ref_o_1[i] = -cylns[ci]->o[2] * cylns[ci]->M[i][2];			
   }
     // Right-circular cylinder:
     prob[0] = gsl_rng_uniform (random);
     prob[1] = gsl_rng_uniform (random);
     SRC_QTL(srcs[ci], prob, s_em);
     RAY_INI_S(srcs[ci]->o, s_em, srcs[ci]->stks, ray, sin_t)
     if ( srcs[ci]->rotate_f )
     {
       ROT_VEC_UNIT(u_rot, ray.u, srcs[ci]->M);
       for (size_t i = 0; i < 3; i++)
         ray.u[i] = u_rot[i];
       COS_TO_SPH_UNIT(ray.s, ray.u, random);
     }

     printf("\n  As emitted:\n");
     ray_fprintf(stdout, &ray, 2);
     printf("\n");
     check_interface_up_airtight_alex(&ray, l, cylns[ci], srcs_2[ci],
       muc, n_rat_aw, iop_na, iop_nw, Fmat, btt_nbr, iop_nw0, ref_o_1, ref_o_2,
       (double const **) M_open, random);
     printf("\n  Reflected:\n");
     ray_fprintf(stdout, &ray, 2);
     printf("\n");

 }
