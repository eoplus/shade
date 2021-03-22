
/*******************************************************************************
 val_intersect_3d.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 This is a driver program to test if the intersection against three-dimensional
 structures work as expected.

 Compile instruction:
 gcc val_intersect_3d.c ../../../geometry.c ../../../structures.c \
 ../../../rotation.c ../../../intersect.c ../../../intersect_cylinder.c \
 ../../../intersect_cuboid.c ../../../intersect_cone.c ../../../accumulators.c \
 ../../../ray.c ../../../memory.c ../../../aux.c ../../../skyrad.c \
 ../../../sources.c -DSHADOWING -DSPATIALLY_RESOLVED -lm -lgsl -O3 \
 -o val_intersect_3d.o

 ./val_intersect_3d.o
*******************************************************************************/

 #include <stdio.h>
 #include <math.h>
 #include <gsl/gsl_blas.h>
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

 void main ( void )
 {
   FILE      *fi;
   long int fpos;
   double origin[3];
   double axis[2];
   double lengths[3];
   double alpha;
   double base_s[3];
   double top_s[3];
   int    closed[2];
   double height;
   double radius;

   double *v_rot = calloc_1d(3, "v_rot in main (val_intersect_2d)");
   double *v     = calloc_1d(3, "v in main (val_intersect_2d)");
   double *u_rot = calloc_1d(3, "u_rot in main (val_intersect_2d)");
   double *u     = calloc_1d(3, "u in main (val_intersect_2d)");
   double *s_rot = calloc_1d(3, "s_rot in main (val_intersect_2d)");
   double *s     = calloc_1d(3, "s in main (val_intersect_2d)");
   double l;

   int res;
   int ci;

   double sclw;
   struct skyradiance *skr = skr_alloc();
   double iop_w0  = 1.0;
   double btt_bh = 1.0;
   int sim_ns = 1;
   int btt_nbr = 1;
   int iop_nw0 = 1;

   double ext  = 1.0;
   double resy = 0.01;
   double resx = 0.01; 
   struct accumulator_bmc * accm = accm_b_alloc();
   accm_b_setup_grid(accm, sim_ns, iop_nw0, btt_nbr, ext, resy, resx);

   struct light_ray ray;
   ray_alloc(1, iop_nw0, btt_nbr, &ray);

   double src_pos[3];

   double stks[STKS_N] = {0.0};
   stks[0] = 1.0;

  /* Tests for cuboids ********************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("Cuboids =======================================================\n%s",
     ANSI_RESET);

   str_cubd **cubds;
   int str_ncb = 4;
   cubds = (str_cubd**) calloc(str_ncb, sizeof(str_cubd*));

   fi = fopen("input_cuboid.txt", "r");
   if(fi == NULL)
   {
     printf("\nERROR: Failed to open input file input_cuboid.txt\n\n");
     exit(-1);
   }
   fpos = ftell(fi);

   for (size_t i = 0; i < str_ncb; i++)
   {
     cubds[i] = str_cubd_alloc();
     str_cubd_read(fi, &fpos, origin, axis, &alpha, lengths, base_s, top_s, 
       closed);
     str_cubd_setup(cubds[i], origin, axis, alpha, lengths, base_s, top_s, 
       closed);   
   }
   printf("\n");
   str_cubds_printf((str_cubd const **) cubds, str_ncb);
   printf("\n");
   for (size_t i = 0; i < str_ncb; i++)
   {
     str_cubd_fprintf(stdout, cubds[i], 1);
     printf("\n");
   }

  /* Cuboid 01 */
   ci   = 0;
   l    = 2.0;
   sclw = 1.0;
   s[0] = 0.0 * RAD;
   s[1] = 0.0 * RAD;
   s[2] = 1.0;


   printf("%d %d %lf %lf\n", accm->nx, accm->ny, accm->xbrks[1], accm->ybrks[1]);

   accm_b_reset(accm);
//printf("GOT HERE 01\n");
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
//printf("GOT HERE 02\n");
       src_pos[0] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[1] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       src_pos[2] = -1.0;
//printf("GOT HERE 03\n");
       ray_ini_s(src_pos, s, stks, &ray);
//printf("GOT HERE 04\n");
       ray_mov(l, &ray);
//printf("GOT HERE 05\n");
       res = !intrs_cubd(ray.a, ray.b, ray.u, cubds[ci]);
//printf("GOT HERE 06\n");
       if ( res )
       {
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
//printf("GOT HERE 07\n");
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cuboid_01_0top", "_out_spr", &iop_w0, 
     &btt_bh, &s[0]);
//printf("GOT HERE 08\n");
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   ci = 0;
   l = 2.0;
   sclw = 1.0;
   s[0] = 90.0 * RAD;
   s[1] = 0.0  * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = -1.0;
       src_pos[1] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[2] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cubd(ray.a, ray.b, ray.u, cubds[ci]);
       if ( res )
       {
         ray.a[0] = ray.a[1];
         ray.a[1] = ray.a[2];
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cuboid_01_1side", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

  /* Cuboid 02 */
   ci = 1;
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
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cubd(ray.a, ray.b, ray.u, cubds[ci]);
       if ( res )
       {
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cuboid_02_0top", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   ci = 1;
   l = 2.0;
   sclw = 1.0;
   s[0] = 90.0 * RAD;
   s[1] = 0.0  * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = -1.0;
       src_pos[1] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[2] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cubd(ray.a, ray.b, ray.u, cubds[ci]);
       if ( res )
       {
         ray.a[0] = ray.a[1];
         ray.a[1] = ray.a[2];
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cuboid_02_1side", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

  /* Cuboid 03 */
   ci = 2;
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
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cubd(ray.a, ray.b, ray.u, cubds[ci]);
       if ( res )
       {
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cuboid_03_0top", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   ci = 2;
   l = 2.0;
   sclw = 1.0;
   s[0] = 90.0 * RAD;
   s[1] = 0.0  * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = -1.0;
       src_pos[1] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[2] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cubd(ray.a, ray.b, ray.u, cubds[ci]);
       if ( res )
       {
         ray.a[0] = ray.a[1];
         ray.a[1] = ray.a[2];
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cuboid_03_1side", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

  /* Cuboid 04 */
   ci = 3;
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
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cubd(ray.a, ray.b, ray.u, cubds[ci]);
       if ( res )
       {
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cuboid_04_0top", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   ci = 3;
   l = 2.0;
   sclw = 1.0;
   s[0] = 90.0 * RAD;
   s[1] = 0.0  * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = -1.0;
       src_pos[1] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[2] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cubd(ray.a, ray.b, ray.u, cubds[ci]);
       if ( res )
       {
         ray.a[0] = ray.a[1];
         ray.a[1] = ray.a[2];
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cuboid_04_1side", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);


  /* Tests for right-cylinders*************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("Right-Cylinders ===============================================\n%s",
     ANSI_RESET);

   str_cyln **cylns;
   int str_ncl = 4;
   cylns = (str_cyln**) calloc(str_ncl, sizeof(str_cyln*));

   fi = fopen("input_cylinder_r.txt", "r");
   if(fi == NULL)
   {
     printf("\nERROR: Failed to open input file input_cylinder_r.txt\n\n");
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

   printf("\n");
   str_cylns_printf( (str_cyln const **) cylns, str_ncl);
   printf("\n");
   for (size_t i = 0; i < str_ncl; i++)
   {
     str_cyln_fprintf(stdout, cylns[i], 1);
     printf("\n");
   }

  /* Right-cylinder 01 */
   ci = 0;
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
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cyln(ray.a, ray.b, ray.u, l, cylns[ci]);
       if ( res )
       {
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cylinder_r_01_0top", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   ci = 0;
   l = 2.0;
   sclw = 1.0;
   s[0] = 90.0 * RAD;
   s[1] = 0.0  * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = -1.0;
       src_pos[1] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[2] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cyln(ray.a, ray.b, ray.u, l, cylns[ci]);
       if ( res )
       {
         ray.a[0] = ray.a[1];
         ray.a[1] = ray.a[2];
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cylinder_r_01_1side", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

  /* Right-cylinder 02 */
   ci = 1;
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
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cyln(ray.a, ray.b, ray.u, l, cylns[ci]);
       if ( res )
       {
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cylinder_r_02_0top", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   ci = 1;
   l = 2.0;
   sclw = 1.0;
   s[0] = 90.0 * RAD;
   s[1] = 0.0  * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = -1.0;
       src_pos[1] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[2] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cyln(ray.a, ray.b, ray.u, l, cylns[ci]);
       if ( res )
       {
         ray.a[0] = ray.a[1];
         ray.a[1] = ray.a[2];
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cylinder_r_02_1side", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

  /* Right-cylinder 03 */
   ci = 2;
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
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cyln(ray.a, ray.b, ray.u, l, cylns[ci]);
       if ( res )
       {
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cylinder_r_03_0top", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   ci = 2;
   l = 2.0;
   sclw = 1.0;
   s[0] = 90.0 * RAD;
   s[1] = 0.0  * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = -1.0;
       src_pos[1] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[2] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cyln(ray.a, ray.b, ray.u, l, cylns[ci]);
       if ( res )
       {
         ray.a[0] = ray.a[1];
         ray.a[1] = ray.a[2];
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cylinder_r_03_1side", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

  /* Right-cylinder 04 */
   ci = 3;
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
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cyln(ray.a, ray.b, ray.u, l, cylns[ci]);
       if ( res )
       {
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cylinder_r_04_0top", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   ci = 3;
   l = 2.0;
   sclw = 1.0;
   s[0] = 90.0 * RAD;
   s[1] = 0.0  * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = -1.0;
       src_pos[1] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[2] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cyln(ray.a, ray.b, ray.u, l, cylns[ci]);
       if ( res )
       {
         ray.a[0] = ray.a[1];
         ray.a[1] = ray.a[2];
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cylinder_r_04_1side", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

  /* Tests for oblique-cylinders **********************************************/
   printf("\n%s", ANSI_BOLD);
   printf("Oblique-Cylinders ==============================================\n%s",
     ANSI_RESET);

//   str_cyln **cylns;
//   int str_ncl = 4;
//   cylns = (str_cyln**) calloc(str_ncl, sizeof(str_cyln*));

   fi = fopen("input_cylinder_o.txt", "r");
   if(fi == NULL)
   {
     printf("\nERROR: Failed to open input file input_cylinder_o.txt\n\n");
     exit(-1);
   }
   fpos = ftell(fi);

   for (size_t i = 0; i < str_ncl; i++)
   {
//     cylns[i] = str_cyln_alloc();
     str_cyln_read(fi, &fpos, origin, axis, &alpha, &radius, &height, base_s, 
       top_s, closed);
     str_cyln_setup(cylns[i], origin, axis, alpha, radius, height, base_s, 
       top_s, closed);   
   }
   printf("\n");
   str_cylns_printf( (str_cyln const **) cylns, str_ncl);
   printf("\n");
   for (size_t i = 0; i < str_ncl; i++)
   {
     str_cyln_fprintf(stdout, cylns[i], 1);
     printf("\n");
   }

  /* Oblique-cylinder 01 */
   ci = 0;
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
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cyln(ray.a, ray.b, ray.u, l, cylns[ci]);
       if ( res )
       {
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cylinder_o_01_0top", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   ci = 0;
   l = 2.0;
   sclw = 1.0;
   s[0] = 90.0 * RAD;
   s[1] = 0.0  * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = -1.0;
       src_pos[1] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[2] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cyln(ray.a, ray.b, ray.u, l, cylns[ci]);
       if ( res )
       {
         ray.a[0] = ray.a[1];
         ray.a[1] = ray.a[2];
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cylinder_o_01_1side", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

  /* Oblique-cylinder 02 */
   ci = 1;
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
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cyln(ray.a, ray.b, ray.u, l, cylns[ci]);
       if ( res )
       {
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cylinder_o_02_0top", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   ci = 1;
   l = 2.0;
   sclw = 1.0;
   s[0] = 90.0 * RAD;
   s[1] = 0.0  * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = -1.0;
       src_pos[1] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[2] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cyln(ray.a, ray.b, ray.u, l, cylns[ci]);
       if ( res )
       {
         ray.a[0] = ray.a[1];
         ray.a[1] = ray.a[2];
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cylinder_o_02_1side", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

  /* Oblique-cylinder 03 */
   ci = 2;
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
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cyln(ray.a, ray.b, ray.u, l, cylns[ci]);
       if ( res )
       {
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cylinder_o_03_0top", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   ci = 2;
   l = 2.0;
   sclw = 1.0;
   s[0] = 90.0 * RAD;
   s[1] = 0.0  * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = -1.0;
       src_pos[1] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[2] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cyln(ray.a, ray.b, ray.u, l, cylns[ci]);
       if ( res )
       {
         ray.a[0] = ray.a[1];
         ray.a[1] = ray.a[2];
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cylinder_o_03_1side", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

  /* Oblique-cylinder 04 */
   ci = 3;
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
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cyln(ray.a, ray.b, ray.u, l, cylns[ci]);
       if ( res )
       {
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cylinder_o_04_0top", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   ci = 3;
   l = 2.0;
   sclw = 1.0;
   s[0] = 90.0 * RAD;
   s[1] = 0.0  * RAD;
   s[2] = 1.0;

   accm_b_reset(accm);
   for (size_t cx = 1; cx < (accm->nx - 1); cx++)
   {
     for (size_t cy = 1; cy < (accm->ny - 1); cy++)
     {
       src_pos[0] = -1.0;
       src_pos[1] = (accm->xbrks[cx] + accm->xbrks[cx+1]) / 2.0;
       src_pos[2] = (accm->ybrks[cy] + accm->ybrks[cy+1]) / 2.0;
       ray_ini_s(src_pos, s, stks, &ray);
       ray_mov(l, &ray);
       res = !intrs_cyln(ray.a, ray.b, ray.u, l, cylns[ci]);
       if ( res )
       {
         ray.a[0] = ray.a[1];
         ray.a[1] = ray.a[2];
         accm_b_add(accm, skr, ray.a, s, (double const ***) ray.stks, sclw, 1, 0);
       }
     }
   }

   accm_b_write_grid (accm, 1, "plots/cylinder_o_04_1side", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

  /* Cylinders from PONDER *************************************************/
   // Light source is isotropic, inside the cylinder.

   // Variables for random sampling:
   double prob[2];
   struct timeval tv;
   gettimeofday (&tv, 0);
   unsigned long int seed = tv.tv_sec + tv.tv_usec;
   gsl_rng *random = gsl_rng_alloc (gsl_rng_ranlxs0);
   gsl_rng_set (random, seed);

   printf("\n%s", ANSI_BOLD);
   printf("PONDER ===============================================\n%s",
     ANSI_RESET);

   origin[0] = 0.0;
   origin[1] = 0.0;
   origin[2] = -0.09;
   axis[0] = 0.0;
   axis[1] = 0.0;
   alpha = 0.0;
   radius = 0.0125;
   height = 0.11;
   base_s[0] = 0.0;
   base_s[1] = 0.0;
   top_s[0] = 0.0;
   top_s[1] = 0.0; 
   closed[0] = 1;
   closed[1] = 0;

   struct source *src = src_alloc();
   double src_stks[4] = {1.0};
   double src_fov = M_PI;
   double src_s[3] = {0.0, 0.0, 1.0};
   double src_ref_o[3] = {0.0, 0.0, -0.065};
   double src_rel_o[3] = {0.0};
   char src_tp[STRMXLEN] = "sir";
   double s_em[3];

   src_setup(src, src_fov, src_s, src_ref_o, src_rel_o, src_stks, src_tp);

   str_cyln_setup(cylns[0], origin, axis, alpha, radius, height, base_s, 
       top_s, closed);   

   str_cyln_fprintf(stdout, cylns[0], 1);
   printf("\n");
   src_fprintf(stdout, src, 1);
   printf("\n");

  /* Starting inside cylinder */
   sclw = 1.0;

   accm_b_reset(accm);
   for (size_t i = 0; i < 10000000; i++)
   {
     prob[0] = gsl_rng_uniform (random);
     prob[1] = gsl_rng_uniform (random);
     SRC_QTL(src, prob, s_em);
     ray_ini_s(src->o, s_em, stks, &ray);
     if ( src->rotate_f )
     {
       ROT_VEC_UNIT(u_rot, ray.u, src->M);
       for (ci = 0; ci < 3; ci++)
         ray.u[ci] = u_rot[ci];
       COS_TO_SPH_UNIT(ray.s, ray.u, random);
     }
     l = 2.0 / ray.u[2];
     ray_mov(l, &ray);
     ray.nsdw_f = !intrs_cyln(ray.a, ray.b, ray.u, l, cylns[0]);
//     if (gsl_rng_uniform (random) > 0.99999)
//       ray_fprintf(stdout, &ray, 1);
     if ( ray.nsdw_f )
     {
       accm_b_add(accm, skr, ray.b, s_em, (double const ***) ray.stks, sclw, 1, 0);
     }
   }

   s[0] = 0.0;
   accm_b_write_grid (accm, 1, "plots/PONDER_01", "_out_spr", &iop_w0, 
   &btt_bh, &s[0]);
   printf ("Result: %sCheck image in plots after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   printf("\nRunning R script val_intersect_3d.R:\n");
   system("Rscript val_intersect_3d.R");
   printf("Done!\n");


 }


