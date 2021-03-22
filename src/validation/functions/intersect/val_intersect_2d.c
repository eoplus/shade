
/*******************************************************************************
 val_intersect_2d.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 This is a driver program to test if the intersection against bi-dimensional
 structures work as expected.

 Compile instruction:
 gcc val_intersect_2d.c ../../../src/geometry.c ../../../src/structures.c ../../../src/rotation.c ../../../src/intersect.c ../../../src/accumulators.c ../../../src/ray.c ../../../src/memory.c ../../../src/aux.c -lm -lgsl -O3 -o val_intersect_2d.o

 ./val_intersect_2d.o
*******************************************************************************/

 #include <stdio.h>
 #include <math.h>
 #include <gsl/gsl_blas.h>
 #include <gsl/gsl_rng.h>
 #include <sys/time.h>

 #include "../../../src/config.h"
 #include "../../../src/aux.h"
 #include "../../../src/memory.h"
 #include "../../../src/constants.h"
 #include "../../../src/geometry.h"
 #include "../../../src/structures.h"
 #include "../../../src/rotation.h"
 #include "../../../src/intersect.h"
 #include "../../../src/skyrad.h"
 #include "../../../src/accumulators.h"
 #include "../../../src/ray.h"

 void main ( void )
 {
   str_rect *rect;
   str_rect **rects;

   FILE      *fi;
   long int fpos;
   double origin[3];
   double axis[2];
   double lengths[2];
   double alpha;

   double *v_rot = calloc_1d(3, "v_rot in main (val_intersect_2d)");
   double *v     = calloc_1d(3, "v in main (val_intersect_2d)");
   double *u_rot = calloc_1d(3, "u_rot in main (val_intersect_2d)");
   double *u     = calloc_1d(3, "u in main (val_intersect_2d)");
   double *s_rot = calloc_1d(3, "s_rot in main (val_intersect_2d)");
   double *s     = calloc_1d(3, "s in main (val_intersect_2d)");
   double l;

   int res;
   int ci;

   double ext  = 1.0;
   double resy = 0.01;
   double resx = 0.01; 
   accm_bmc *accm = accm_b_grid_alloc(1, 1, 1, ext, resy, resx);

   double sclw;
   skyrad *skr;
   double btt_nbr = 1.0;
   double iop_w0  = 1.0;
   double iop_nw0 = 1.0;

   light_ray* ray = ray_alloc(1, iop_nw0, btt_nbr);

   double src_pos[3];

  /* Tests for rectangles *****************************************************/
   int n_set01 = 5;
   int *res_set01 = (int *) calloc(n_set01, sizeof(int));

   printf("\n%s", ANSI_BOLD);
   printf("Rectangles ====================================================\n%s",
     ANSI_RESET);

   int str_nrt = 2;
   rects = (str_rect**) calloc(str_nrt, sizeof(str_rect*));

   fi = fopen("input_rectangle.txt", "r");
   if(fi == NULL)
   {
     printf("\nERROR: Failed to open input file input_rectangle.txt\n\n");
     exit(-1);
   }
   fpos = ftell(fi);

   for (size_t i = 0; i < str_nrt; i++)
   {
     rects[i] = str_rect_alloc();
     str_rect_read(fi, &fpos, origin, axis, lengths, &alpha);   
     str_rect_setup(rects[i], origin, axis, lengths, alpha);   
   }
   printf("\n");
   str_rects_printf(rects, str_nrt);
   printf("\n");
   str_rect_printf(rects[0]);
   printf("\n");
   str_rect_printf(rects[1]);
   printf("\n");

  /* Set 01 */
   l = 2.0;
   v[0] =  0.0;
   v[1] =  0.0;
   v[2] = -1.0;
   s[0] =  0.0 * RAD;
   s[1] =  0.0 * RAD;
   s[2] =  1.0;
   sph_to_cos(u, s);

   ci = 0;
   rot_vec(v_rot, v, rects[0]->o, (double const **) rects[0]->M);
   rot_vec_unit(u_rot, u, (double const **) rects[0]->M);
   res = intrs_rect(v_rot, u_rot, l, rects[0]);
   res_set01[ci] = (res == 1) ? 1 : 0;
   printf("  Set 01: ------------------------------------------------------\n");
   printf("  Ray straight down. Should intersect.                   %s%s%s\n", 
     res_set01[ci] ? ANSI_BLUE : ANSI_RED, 
     res_set01[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

  /* Set 02 */
   ci = 1;
   rot_vec(v_rot, v, rects[1]->o, (double const **) rects[1]->M);
   rot_vec_unit(u_rot, u, (double const **) rects[1]->M);
   res = intrs_rect(v_rot, u_rot, l, rects[1]);
   res_set01[ci] = (res == 0) ? 1 : 0;
   printf("  Set 02: ------------------------------------------------------\n");
   printf("  Ray straight down, to the side. Should not intersect.  %s%s%s\n", 
     res_set01[ci] ? ANSI_BLUE : ANSI_RED, 
     res_set01[ci] ? "PASSED!" : "FAILED!", ANSI_RESET);

  /* Set 03 test 01*/
   l = 2.0;
   sclw = 1.0;
   s[0] = 0.0 * RAD;
   s[1] = 0.0 * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[1] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, s, ray);
       ray_mov(l, ray);
       u_rot[0] = ray->u[0];
       u_rot[1] = ray->u[1];
       u_rot[2] = ray->u[2];
       v_rot[0] = ray->a[0];
       v_rot[1] = ray->a[1];
       v_rot[2] = ray->a[2];
       if ( rects[0]->rotate_f )
       {
         rot_vec_unit(u_rot, u, (double const **) rects[0]->M);
         rot_vec(v_rot, ray->a, rects[0]->o, (double const **) rects[0]->M);
       }
       res = !intrs_rect(v_rot, u_rot, l, rects[0]);
       if ( res )
       {
         accm_b_grid(ray->a, s, ray->stks, accm, skr, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/rect_set_03_test_01", "_out_spr", &iop_w0, 
   &btt_nbr, &s[0]);
   printf ("Result: %sCheck image in intersect_cylnd_res after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

  /* Set 03 test 02*/
  l = 2.0;
   sclw = 1.0;
   s[0] = 45.0 * RAD;
   s[1] = 0.0 * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[1] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       src_pos[2] = -0.2;
       ray_ini(src_pos, s, ray);
       ray_mov(l, ray);
       u_rot[0] = ray->u[0];
       u_rot[1] = ray->u[1];
       u_rot[2] = ray->u[2];
       v_rot[0] = ray->a[0];
       v_rot[1] = ray->a[1];
       v_rot[2] = ray->a[2];
       if ( rects[1]->rotate_f )
       {
         rot_vec_unit(u_rot, ray->u, (double const **) rects[1]->M);
         rot_vec(v_rot, ray->a, rects[1]->o, (double const **) rects[1]->M);
       }
       res = !intrs_rect(v_rot, u_rot, l, rects[1]);
       if ( res )
       {
         accm_b_grid(ray->a, s, ray->stks, accm, skr, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/rect_set_03_test_02", "_out_spr", &iop_w0, 
   &btt_nbr, &s[0]);
   printf ("Result: %sCheck image in intersect_cylnd_res after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);


  /* Tests for ellipse ********************************************************/
   int n_set02 = 5;
   int *res_set02 = (int *) calloc(n_set02, sizeof(int));

   printf("\n%s", ANSI_BOLD);
   printf("Ellipses ======================================================\n%s",
     ANSI_RESET);

   int str_nel = 5;
   str_ellp **ellps = (str_ellp**) calloc(str_nel, sizeof(str_ellp*));

   fi = fopen("input_ellipse.txt", "r");
   if(fi == NULL)
   {
     printf("\nERROR: Failed to open input file input_ellipse.txt\n\n");
     exit(-1);
   }
   fpos = ftell(fi);

   for (size_t i = 0; i < str_nel; i++)
   {
     ellps[i] = str_ellp_alloc();
     str_ellp_read(fi, &fpos, origin, axis, lengths, &alpha);   
     str_ellp_setup(ellps[i], origin, axis, lengths, alpha);   
   }
   printf("\n");
   str_ellps_printf(ellps, str_nel);
   printf("\n");
   str_ellp_printf(ellps[0]);
   printf("\n");
   str_ellp_printf(ellps[1]);
   printf("\n");

  /* Set 01 test 01 */
   l = 2.0;
   sclw = 1.0;
   s[0] = 0.0 * RAD;
   s[1] = 0.0 * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[1] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, s, ray);
       ray_mov(l, ray);
       u_rot[0] = ray->u[0];
       u_rot[1] = ray->u[1];
       u_rot[2] = ray->u[2];
       v_rot[0] = ray->a[0];
       v_rot[1] = ray->a[1];
       v_rot[2] = ray->a[2];
       if ( ellps[0]->rotate_f )
       {
         rot_vec_unit(u_rot, u, (double const **) ellps[0]->M);
         rot_vec(v_rot, ray->a, ellps[0]->o, (double const **) ellps[0]->M);
       }
       res = !intrs_ellp(v_rot, u_rot, l, ellps[0]);
       if ( res )
       {
         accm_b_grid(ray->a, s, ray->stks, accm, skr, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/ellp_set_01_test_01", "_out_spr", &iop_w0, 
   &btt_nbr, &s[0]);
   printf ("Result: %sCheck image in intersect_cylnd_res after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

  /* Set 01 test 02 */
   l = 2.0;
   sclw = 1.0;
   s[0] = 0.0 * RAD;
   s[1] = 0.0 * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[1] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, s, ray);
       ray_mov(l, ray);
       u_rot[0] = ray->u[0];
       u_rot[1] = ray->u[1];
       u_rot[2] = ray->u[2];
       v_rot[0] = ray->a[0];
       v_rot[1] = ray->a[1];
       v_rot[2] = ray->a[2];
       if ( ellps[1]->rotate_f )
       {
         rot_vec_unit(u_rot, ray->u, (double const **) ellps[1]->M);
         rot_vec(v_rot, ray->a, ellps[1]->o, (double const **) ellps[1]->M);
       }
       res = !intrs_ellp(v_rot, u_rot, l, ellps[1]);
       if ( res )
       {
         accm_b_grid(ray->a, s, ray->stks, accm, skr, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/ellp_set_01_test_02", "_out_spr", &iop_w0, 
   &btt_nbr, &s[0]);
   printf ("Result: %sCheck image in intersect_cylnd_res after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

  /* Set 01 test 03 */
   l = 2.0;
   sclw = 1.0;
   s[0] = 0.0 * RAD;
   s[1] = 0.0 * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[1] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, s, ray);
       ray_mov(l, ray);
       u_rot[0] = ray->u[0];
       u_rot[1] = ray->u[1];
       u_rot[2] = ray->u[2];
       v_rot[0] = ray->a[0];
       v_rot[1] = ray->a[1];
       v_rot[2] = ray->a[2];
       if ( ellps[1]->rotate_f )
       {
         rot_vec_unit(u_rot, ray->u, (double const **) ellps[2]->M);
         rot_vec(v_rot, ray->a, ellps[2]->o, (double const **) ellps[2]->M);
       }
       res = !intrs_ellp(v_rot, u_rot, l, ellps[2]);
       if ( res )
       {
         accm_b_grid(ray->a, s, ray->stks, accm, skr, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/ellp_set_01_test_03", "_out_spr", &iop_w0, 
   &btt_nbr, &s[0]);
   printf ("Result: %sCheck image in intersect_cylnd_res after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

  /* Set 01 test 04 */
   l = 2.0;
   sclw = 1.0;
   s[0] = 0.0 * RAD;
   s[1] = 0.0 * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[1] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, s, ray);
       ray_mov(l, ray);
       u_rot[0] = ray->u[0];
       u_rot[1] = ray->u[1];
       u_rot[2] = ray->u[2];
       v_rot[0] = ray->a[0];
       v_rot[1] = ray->a[1];
       v_rot[2] = ray->a[2];
       if ( ellps[1]->rotate_f )
       {
         rot_vec_unit(u_rot, ray->u, (double const **) ellps[3]->M);
         rot_vec(v_rot, ray->a, ellps[3]->o, (double const **) ellps[3]->M);
       }
       res = !intrs_ellp(v_rot, u_rot, l, ellps[3]);
       if ( res )
       {
         accm_b_grid(ray->a, s, ray->stks, accm, skr, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/ellp_set_01_test_04", "_out_spr", &iop_w0, 
   &btt_nbr, &s[0]);
   printf ("Result: %sCheck image in intersect_cylnd_res after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

  /* Set 01 test 05 */
   l = 2.0;
   sclw = 1.0;
   s[0] = 0.0 * RAD;
   s[1] = 0.0 * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[1] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, s, ray);
       ray_mov(l, ray);
       u_rot[0] = ray->u[0];
       u_rot[1] = ray->u[1];
       u_rot[2] = ray->u[2];
       v_rot[0] = ray->a[0];
       v_rot[1] = ray->a[1];
       v_rot[2] = ray->a[2];
       if ( ellps[1]->rotate_f )
       {
         rot_vec_unit(u_rot, ray->u, (double const **) ellps[4]->M);
         rot_vec(v_rot, ray->a, ellps[4]->o, (double const **) ellps[4]->M);
       }
       res = !intrs_ellp(v_rot, u_rot, l, ellps[4]);
       if ( res )
       {
         accm_b_grid(ray->a, s, ray->stks, accm, skr, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/ellp_set_01_test_05", "_out_spr", &iop_w0, 
   &btt_nbr, &s[0]);
   printf ("Result: %sCheck image in intersect_cylnd_res after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

  printf("\nRunning R script val_intersect_2d.R:\n");
   system("Rscript val_intersect_2d.R");
   printf("Done!\n");

 }


