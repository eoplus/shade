
 #include <stdio.h>
 #include <stdlib.h>
 #include <math.h>

 #include "../../../src/aux.h"
 #include "../../../src/geometry.h"
 #include "../../../src/constants.h"
 #include "../../../src/structures.h"
 #include "../../../src/ray.h"
 #include "../../../src/intersect.h"
 #include "../../../src/acm_geom.h"
 #include "../../../src/config.h"
 #include "../../../src/output.h"

 #ifndef SPRES
   #error SPRES must be defined in config.h
 #endif

 void main (void)
 {
   FILE   *p_fi;
   long int fpos;

   int    i, ci, cj;
   int    str_nel;
   struct str_ellp *p_ellps;

   double origin[3];
   double axis[2];
   double radius[2];

   double s;
   double psi;
   double phi;
   double src_pos[3] = {0.0};
   int    sim_nbr = 1;
   int    sim_nw0 = 1;
   struct str_ray ray, ray_vrt;
   ray.ppp_stks = calloc_3d(1, 1, 4, "ray.ppp_stks in val_intersect_cylinder");
   ray_vrt.ppp_stks = calloc_3d(1, 1, 4, "ray_vrt.ppp_stks in val_intersect_cylinder");

   int sdwf;

   double ext  = 1.0;
   double resy = 0.01;
   double resx = 0.01; 
   struct str_acm dirw;
   acm_grid_alloc(1, 1, 1, ext, resy, resx, &dirw);
   double sclw;
   struct str_skr skr;
   double btthdr = 1.0;
   double iop_w0 = 1.0;

   // Set 1: -------------------------------------------------------------------
   int n_set01 = 16;
   int *res_set01 = (int *) calloc(n_set01, sizeof(int));

   printf("\n%s", ANSI_BOLD);
   printf("Set 01 =========================================================\n");
   printf("Z-axis oriented circle. \n\n%s", ANSI_RESET);

   str_nel = 1;
   p_ellps = (struct str_ellp*) calloc(str_nel, sizeof(struct str_ellp));

   origin[0] = 0.0;
   origin[1] = 0.0;
   origin[2] = 0.0;
   axis[0]   = 0.0 * RAD;
   axis[1]   = 0.0 * RAD;
   radius[0] = 1.0;
   radius[1] = 1.0;

   setup_ellp (&p_ellps[0], origin, axis, radius);
   print_ellp (&p_ellps[0]);

   // Set 1, test 1 ------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 01, Test 01 ------------------------------------------------\n");
   printf("Ray segment starting from origin, straight downwards. Should not "
     "intersect. \n\n%s", ANSI_RESET);

   i = 0;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = 0.0;
   psi = 0.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_circle (ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set01[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set01[i] ? ANSI_BLUE : ANSI_RED, 
     res_set01[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 1, test 2 ------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 01, Test 02 ------------------------------------------------\n");
   printf("Ray segment starting from origin, at an angle of 30º. Should not "
     "intersect. \n\n%s", ANSI_RESET);

   i = 1;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = 0.0;
   psi = 30.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_circle (ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set01[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set01[i] ? ANSI_BLUE : ANSI_RED, 
     res_set01[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 1, test 3 ------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 01, Test 03 ------------------------------------------------\n");
   printf("Ray segment starting after origin stright down. Should not "
     "intersect. \n\n%s", ANSI_RESET);

   i = 2;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = 0.5;
   psi = 0.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_circle (ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set01[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set01[i] ? ANSI_BLUE : ANSI_RED, 
     res_set01[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 1, test 4 ------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 01, Test 04 ------------------------------------------------\n");
   printf("Ray segment starting after origin, at an angle of 30º. Should not "
     "intersect. \n\n%s", ANSI_RESET);

   i = 3;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = 0.5;
   psi = 30.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_circle(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set01[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set01[i] ? ANSI_BLUE : ANSI_RED, 
     res_set01[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 1, test 5 ------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 01, Test 05 ------------------------------------------------\n");
   printf("Ray segment starting before origin stright up. Should not "
     "intersect. \n\n%s", ANSI_RESET);

   i = 4;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = -0.5;
   psi = 180.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_circle (ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set01[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set01[i] ? ANSI_BLUE : ANSI_RED, 
     res_set01[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 1, test 6 ------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 01, Test 06 ------------------------------------------------\n");
   printf("Ray segment starting before origin, at an angle of 150º. Should not "
     "intersect. \n\n%s", ANSI_RESET);

   i = 5;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = -0.5;
   psi = 150.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_circle(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set01[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set01[i] ? ANSI_BLUE : ANSI_RED, 
     res_set01[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 1, test 7 ------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 01, Test 07 ------------------------------------------------\n");
   printf("Ray segment starting after origin, moving sideways. Should not "
     "intersect. \n\n%s", ANSI_RESET);

   i = 6;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = -0.5;
   psi = 90.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_circle(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set01[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set01[i] ? ANSI_BLUE : ANSI_RED, 
     res_set01[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 1, test 8 ------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 01, Test 08 ------------------------------------------------\n");
   printf("Ray segment starting before origin, moving sideways. Should not "
     "intersect. \n\n%s", ANSI_RESET);

   i = 7;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = -0.5;
   psi = 90.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_circle(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set01[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set01[i] ? ANSI_BLUE : ANSI_RED, 
     res_set01[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 1, test 9 ------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 01, Test 09 ------------------------------------------------\n");
   printf("Ray segment starting before origin, straight down. Should "
     "intersect. \n\n%s", ANSI_RESET);

   i = 8;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = -0.5;
   psi = 0.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_circle(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set01[i] = (sdwf == 1) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set01[i] ? ANSI_BLUE : ANSI_RED, 
     res_set01[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 1, test 10 -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 01, Test 10 ------------------------------------------------\n");
   printf("Ray segment starting after origin, straight up. Should "
     "intersect. \n\n%s", ANSI_RESET);

   i = 9;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = 0.5;
   psi = 180.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_circle(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set01[i] = (sdwf == 1) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set01[i] ? ANSI_BLUE : ANSI_RED, 
     res_set01[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 1, test 11 -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 01, Test 11 ------------------------------------------------\n");
   printf("Ray segment starting before origin, at an angle of 30º. Should "
     "intersect. \n\n%s", ANSI_RESET);

   i = 10;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = -0.5;
   psi = 30.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_circle(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set01[i] = (sdwf == 1) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set01[i] ? ANSI_BLUE : ANSI_RED, 
     res_set01[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 1, test 12 -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 01, Test 12 ------------------------------------------------\n");
   printf("Ray segment starting after origin, at an angle of 125º. Should "
     "intersect. \n\n%s", ANSI_RESET);

   i = 11;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = 0.5;
   psi = 125.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_circle(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set01[i] = (sdwf == 1) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set01[i] ? ANSI_BLUE : ANSI_RED, 
     res_set01[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 1, test 13 -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 01, Test 13 ------------------------------------------------\n");
   printf("Ray segment starting after origin, straight up tangent to border. Should "
     "intersect. \n\n%s", ANSI_RESET);

   i = 12;
   src_pos[0] = 1.0;
   src_pos[1] = 0.0;
   src_pos[2] = 1.0;
   psi = 180.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_circle(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set01[i] = (sdwf == 1) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set01[i] ? ANSI_BLUE : ANSI_RED, 
     res_set01[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 1, test 14 -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 01, Test 14 ------------------------------------------------\n");
   printf("Ray segment starting after origin to its side, straight up. Should "
     " not intersect. \n\n%s", ANSI_RESET);

   i = 13;
   src_pos[0] = 1.5;
   src_pos[1] = 0.0;
   src_pos[2] = 0.5;
   psi = 180.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_circle(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set01[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set01[i] ? ANSI_BLUE : ANSI_RED, 
     res_set01[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 1, test 15 -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 01, Test 15 ------------------------------------------------\n");
   printf("Ray segment starting after origin, tangent path to border. Should "
     "intersect. \n\n%s", ANSI_RESET);

   i = 14;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = 1.0;
   psi = 135.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0 / -cos(psi);

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_circle(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set01[i] = (sdwf == 1) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set01[i] ? ANSI_BLUE : ANSI_RED, 
     res_set01[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 1, test 16 -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 01, Test 16 ------------------------------------------------\n");
   printf("Ray segment starting after origin, at an angle of 135º. Should "
     "not intersect. \n\n%s", ANSI_RESET);

   i = 15;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = 1.5;
   psi = 135.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.5 / -cos(psi);

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_circle(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set01[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set01[i] ? ANSI_BLUE : ANSI_RED, 
     res_set01[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 2: -------------------------------------------------------------------
   int n_set02 = 18;
   int *res_set02 = (int *) calloc(n_set02, sizeof(int));

   printf("\n%s", ANSI_BOLD);
   printf("Set 02 =========================================================\n");
   printf("Z-axis oriented ellipse. \n\n%s", ANSI_RESET);

   str_nel = 1;
   p_ellps = (struct str_ellp*) calloc(str_nel, sizeof(struct str_ellp));

   origin[0] = 0.0;
   origin[1] = 0.0;
   origin[2] = 0.0;
   axis[0]   = 0.0 * RAD;
   axis[1]   = 0.0 * RAD;
   radius[0] = 1.0;
   radius[1] = 2.0;

   setup_ellp (&p_ellps[0], origin, axis, radius);
   print_ellp (&p_ellps[0]);

   // Set 2, test 1 ------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 02, Test 01 ------------------------------------------------\n");
   printf("Ray segment starting from origin, straight downwards. Should not "
     "intersect. \n\n%s", ANSI_RESET);

   i = 0;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = 0.0;
   psi = 0.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_ellipse (ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set02[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set02[i] ? ANSI_BLUE : ANSI_RED, 
     res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 2, test 2 ------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 02, Test 02 ------------------------------------------------\n");
   printf("Ray segment starting from origin, at an angle of 30º. Should not "
     "intersect. \n\n%s", ANSI_RESET);

   i = 1;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = 0.0;
   psi = 30.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_ellipse (ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set02[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set02[i] ? ANSI_BLUE : ANSI_RED, 
     res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 2, test 3 ------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 02, Test 03 ------------------------------------------------\n");
   printf("Ray segment starting after origin stright down. Should not "
     "intersect. \n\n%s", ANSI_RESET);

   i = 2;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = 0.5;
   psi = 0.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_ellipse (ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set02[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set02[i] ? ANSI_BLUE : ANSI_RED, 
     res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 2, test 4 ------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 02, Test 04 ------------------------------------------------\n");
   printf("Ray segment starting after origin, at an angle of 30º. Should not "
     "intersect. \n\n%s", ANSI_RESET);

   i = 3;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = 0.5;
   psi = 30.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set02[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set02[i] ? ANSI_BLUE : ANSI_RED, 
     res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 2, test 5 ------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 02, Test 05 ------------------------------------------------\n");
   printf("Ray segment starting before origin stright up. Should not "
     "intersect. \n\n%s", ANSI_RESET);

   i = 4;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = -0.5;
   psi = 180.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set02[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set02[i] ? ANSI_BLUE : ANSI_RED, 
     res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 2, test 6 ------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 02, Test 06 ------------------------------------------------\n");
   printf("Ray segment starting before origin, at an angle of 150º. Should not "
     "intersect. \n\n%s", ANSI_RESET);

   i = 5;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = -0.5;
   psi = 150.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set02[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set02[i] ? ANSI_BLUE : ANSI_RED, 
     res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 2, test 7 ------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 02, Test 07 ------------------------------------------------\n");
   printf("Ray segment starting after origin, moving sideways. Should not "
     "intersect. \n\n%s", ANSI_RESET);

   i = 6;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = -0.5;
   psi = 90.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set02[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set02[i] ? ANSI_BLUE : ANSI_RED, 
     res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 2, test 8 ------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 02, Test 08 ------------------------------------------------\n");
   printf("Ray segment starting before origin, moving sideways. Should not "
     "intersect. \n\n%s", ANSI_RESET);

   i = 7;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = -0.5;
   psi = 90.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set02[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set02[i] ? ANSI_BLUE : ANSI_RED, 
     res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 2, test 9 ------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 02, Test 09 ------------------------------------------------\n");
   printf("Ray segment starting before origin, straight down. Should "
     "intersect. \n\n%s", ANSI_RESET);

   i = 8;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = -0.5;
   psi = 0.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set02[i] = (sdwf == 1) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set02[i] ? ANSI_BLUE : ANSI_RED, 
     res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 2, test 10 -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 02, Test 10 ------------------------------------------------\n");
   printf("Ray segment starting after origin, straight up. Should "
     "intersect. \n\n%s", ANSI_RESET);

   i = 9;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = 0.5;
   psi = 180.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set02[i] = (sdwf == 1) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set02[i] ? ANSI_BLUE : ANSI_RED, 
     res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 2, test 11 -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 02, Test 11 ------------------------------------------------\n");
   printf("Ray segment starting before origin, at an angle of 30º. Should "
     "intersect. \n\n%s", ANSI_RESET);

   i = 10;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = -0.5;
   psi = 30.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set02[i] = (sdwf == 1) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set02[i] ? ANSI_BLUE : ANSI_RED, 
     res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 2, test 12 -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 02, Test 12 ------------------------------------------------\n");
   printf("Ray segment starting after origin, at an angle of 125º. Should "
     "intersect. \n\n%s", ANSI_RESET);

   i = 11;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = 0.5;
   psi = 125.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set02[i] = (sdwf == 1) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set02[i] ? ANSI_BLUE : ANSI_RED, 
     res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 2, test 13 -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 02, Test 13 ------------------------------------------------\n");
   printf("Ray segment starting after origin, straight up tangent to border. Should "
     "intersect. \n\n%s", ANSI_RESET);

   i = 12;
   src_pos[0] = 1.0;
   src_pos[1] = 0.0;
   src_pos[2] = 1.0;
   psi = 180.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set02[i] = (sdwf == 1) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set02[i] ? ANSI_BLUE : ANSI_RED, 
     res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 2, test 14 -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 02, Test 14 ------------------------------------------------\n");
   printf("Ray segment starting after origin to its side, straight up. Should "
     " not intersect. \n\n%s", ANSI_RESET);

   i = 13;
   src_pos[0] = 1.5;
   src_pos[1] = 0.0;
   src_pos[2] = 0.5;
   psi = 180.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set02[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set02[i] ? ANSI_BLUE : ANSI_RED, 
     res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 2, test 15 -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 02, Test 15 ------------------------------------------------\n");
   printf("Ray segment starting after origin, tangent path to border. Should "
     "intersect. \n\n%s", ANSI_RESET);

   i = 14;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = 1.0;
   psi = 135.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0 / -cos(psi);

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set02[i] = (sdwf == 1) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set02[i] ? ANSI_BLUE : ANSI_RED, 
     res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 2, test 16 -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 02, Test 16 ------------------------------------------------\n");
   printf("Ray segment starting after origin, at an angle of 135º. Should "
     "not intersect. \n\n%s", ANSI_RESET);

   i = 15;
   src_pos[0] = 0.0;
   src_pos[1] = 0.0;
   src_pos[2] = 1.5;
   psi = 135.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.5 / -cos(psi);

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set02[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set02[i] ? ANSI_BLUE : ANSI_RED, 
     res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 2, test 17 -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 02, Test 17 ------------------------------------------------\n");
   printf("Ray segment starting after origin, tangent path to border in Y. Should "
     "intersect. \n\n%s", ANSI_RESET);

   i = 16;
   src_pos[0] = 0.0;
   src_pos[1] = 2.0;
   src_pos[2] = 1.0;
   psi = 180.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set02[i] = (sdwf == 1) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set02[i] ? ANSI_BLUE : ANSI_RED, 
     res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 2, test 18 -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 02, Test 18 ------------------------------------------------\n");
   printf("Ray segment starting after origin to the side. Should not "
     "intersect. \n\n%s", ANSI_RESET);

   i = 17;
   src_pos[0] = 0.0;
   src_pos[1] = 2.1;
   src_pos[2] = 1.0;
   psi = 180.0 * RAD;
   phi = 0.0 * RAD;
   s = 1.0;

   ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
   ray_mov(s, &ray);
   sdwf = intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);

   print_ray(&ray);
   printf ("\nIntersection ellipse flag: %d\n\n", sdwf);
   res_set02[i] = (sdwf == 0) ? 1 : 0;
   printf ("Result: %s %s %s\n\n", res_set02[i] ? ANSI_BLUE : ANSI_RED, 
     res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);

   // Set 3: -------------------------------------------------------------------
   int n_set03 = 4;
   int *res_set03 = (int *) calloc(n_set03, sizeof(int));

   printf("\n%s", ANSI_BOLD);
   printf("Set 03 =========================================================\n");
   printf("Z-axis oriented circle, imaging approach. \n\n%s", ANSI_RESET);

   str_nel = 1;
   p_ellps = (struct str_ellp*) calloc(str_nel, sizeof(struct str_ellp));

   origin[0] = 0.0;
   origin[1] = 0.0;
   origin[2] = 0.0;
   axis[0]   = 0.0 * RAD;
   axis[1]   = 0.0 * RAD;
   radius[0] = 0.5;
   radius[1] = 0.5;

   setup_ellp (&p_ellps[0], origin, axis, radius);
   print_ellp (&p_ellps[0]);

   // Set 3, test 1: -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 03, Test 01 ------------------------------------------------\n");
   printf("Rays from 0 theta and 0 phi. \n\n%s", ANSI_RESET);

   acm_reset(&dirw);
   psi  = 0.0 * RAD;
   phi  = 0.0 * RAD;
   s    = 10.0;
   sclw = 1.0;

   for (ci = 1; ci < (dirw.nx - 1); ci++)
   {
     for (cj = 1; cj < (dirw.ny - 1); cj++)
     {
       src_pos[0] = (dirw.p_xbrks[ci] + dirw.p_xbrks[ci+1]) / 2.0;
       src_pos[1] = (dirw.p_ybrks[cj] + dirw.p_ybrks[cj+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
       ray_mov(s, &ray);
       sdwf = !intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);
       if ( sdwf )
       {
         acm_grid(ray.ppos, psi, phi, ray.ppp_stks, &dirw, &skr, sclw, 1, 0);
       }
     }
   }

   write_grid (&dirw, 1, "plots/ellipse/set_03_test_01", "_out_spr",  
     &iop_w0, &btthdr, &psi);
   printf ("Result: %sCheck image in plots/ellipse after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   // Set 3, test 2: -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 03, Test 02 ------------------------------------------------\n");
   printf("Rays from 15 theta and 0 phi. \n\n%s", ANSI_RESET);

   acm_reset(&dirw);
   psi  = 15.0 * RAD;
   phi  = 0.0 * RAD;
   s    = 10.0;
   sclw = 1.0;

   for (ci = 1; ci < (dirw.nx - 1); ci++)
   {
     for (cj = 1; cj < (dirw.ny - 1); cj++)
     {
       src_pos[0] = (dirw.p_xbrks[ci] + dirw.p_xbrks[ci+1]) / 2.0;
       src_pos[1] = (dirw.p_ybrks[cj] + dirw.p_ybrks[cj+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
       ray_mov(s, &ray);
       sdwf = !intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);
       if ( sdwf )
       {
         acm_grid(ray.ppos, psi, phi, ray.ppp_stks, &dirw, &skr, sclw, 1, 0);
       }
     }
   }

   write_grid (&dirw, 1, "plots/ellipse/set_03_test_02", "_out_spr",  
     &iop_w0, &btthdr, &psi);
   printf ("Result: %sCheck image in plots/ellipse after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   // Set 3, test 3: -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 03, Test 03 ------------------------------------------------\n");
   printf("Rays from 30 theta and 0 phi. \n\n%s", ANSI_RESET);

   acm_reset(&dirw);
   psi  = 30.0 * RAD;
   phi  = 0.0 * RAD;
   s    = 10.0;
   sclw = 1.0;

   for (ci = 1; ci < (dirw.nx - 1); ci++)
   {
     for (cj = 1; cj < (dirw.ny - 1); cj++)
     {
       src_pos[0] = (dirw.p_xbrks[ci] + dirw.p_xbrks[ci+1]) / 2.0;
       src_pos[1] = (dirw.p_ybrks[cj] + dirw.p_ybrks[cj+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
       ray_mov(s, &ray);
       sdwf = !intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);
       if ( sdwf )
       {
         acm_grid(ray.ppos, psi, phi, ray.ppp_stks, &dirw, &skr, sclw, 1, 0);
       }
     }
   }

   write_grid (&dirw, 1, "plots/ellipse/set_03_test_03", "_out_spr",  
     &iop_w0, &btthdr, &psi);
   printf ("Result: %sCheck image in plots/ellipse after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   // Set 3, test 4: -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 03, Test 04 ------------------------------------------------\n");
   printf("Rays from 30 theta and 90 phi. \n\n%s", ANSI_RESET);

   acm_reset(&dirw);
   psi  = 30.0 * RAD;
   phi  = 90.0 * RAD;
   s    = 10.0;
   sclw = 1.0;

   for (ci = 1; ci < (dirw.nx - 1); ci++)
   {
     for (cj = 1; cj < (dirw.ny - 1); cj++)
     {
       src_pos[0] = (dirw.p_xbrks[ci] + dirw.p_xbrks[ci+1]) / 2.0;
       src_pos[1] = (dirw.p_ybrks[cj] + dirw.p_ybrks[cj+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
       ray_mov(s, &ray);
       sdwf = !intrs_ellipse(ray.ppos, ray.cdir, s, &p_ellps[0]);
       if ( sdwf )
       {
         acm_grid(ray.ppos, psi, phi, ray.ppp_stks, &dirw, &skr, sclw, 1, 0);
       }
     }
   }

   write_grid (&dirw, 1, "plots/ellipse/set_03_test_04", "_out_spr",  
     &iop_w0, &btthdr, &psi);
   printf ("Result: %sCheck image in plots/ellipse after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   // Set 4: -------------------------------------------------------------------
   int n_set04 = 6;
   int *res_set04 = (int *) calloc(n_set04, sizeof(int));

   printf("\n%s", ANSI_BOLD);
   printf("Set 04 =========================================================\n");
   printf("Off-axis circle, 30º 0º imaging approach. \n\n%s", ANSI_RESET);

   str_nel = 1;
   p_ellps = (struct str_ellp*) calloc(str_nel, sizeof(struct str_ellp));

   origin[0] = 0.0;
   origin[1] = 0.0;
   origin[2] = 0.0;
   axis[0]   = 45.0 * RAD;
   axis[1]   = 0.0 * RAD;
   radius[0] = 0.5;
   radius[1] = 0.5;

   setup_ellp (&p_ellps[0], origin, axis, radius);
   print_ellp (&p_ellps[0]);

   // Set 4, test 1: -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 04, Test 01 ------------------------------------------------\n");
   printf("Rays from 0 theta and 0 phi. \n\n%s", ANSI_RESET);

   acm_reset(&dirw);
   psi  = 0.0 * RAD;
   phi  = 0.0 * RAD;
   s    = 10.0;
   sclw = 1.0;

   for (ci = 1; ci < (dirw.nx - 1); ci++)
   {
     for (cj = 1; cj < (dirw.ny - 1); cj++)
     {
       src_pos[0] = (dirw.p_xbrks[ci] + dirw.p_xbrks[ci+1]) / 2.0;
       src_pos[1] = (dirw.p_ybrks[cj] + dirw.p_ybrks[cj+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
       ray_mov(s, &ray);
       copy_vec(ray_vrt.ppos, ray.ppos, 3);
       copy_vec(ray_vrt.cdir, ray.cdir, 3);
       center_vec(ray_vrt.ppos, p_ellps[0].p_origin, 3);
       rotate_vec(ray_vrt.ppos, p_ellps[0].p_caxs, 0);
       rotate_vec(ray_vrt.cdir, p_ellps[0].p_caxs, 0);
       sdwf = !intrs_ellipse(ray_vrt.ppos, ray_vrt.cdir, s, &p_ellps[0]);
       if ( sdwf )
       {
         acm_grid(ray.ppos, psi, phi, ray.ppp_stks, &dirw, &skr, sclw, 1, 0);
       }
     }
   }

   write_grid (&dirw, 1, "plots/ellipse/set_04_test_01", "_out_spr",  
     &iop_w0, &btthdr, &psi);
   printf ("Result: %sCheck image in plots/ellipse after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   // Set 4, test 2: -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 04, Test 02 ------------------------------------------------\n");
   printf("Rays from 15 theta and 0 phi. \n\n%s", ANSI_RESET);

   acm_reset(&dirw);
   psi  = 15.0 * RAD;
   phi  =  0.0 * RAD;
   s    = 10.0;
   sclw = 1.0;

   for (ci = 1; ci < (dirw.nx - 1); ci++)
   {
     for (cj = 1; cj < (dirw.ny - 1); cj++)
     {
       src_pos[0] = (dirw.p_xbrks[ci] + dirw.p_xbrks[ci+1]) / 2.0;
       src_pos[1] = (dirw.p_ybrks[cj] + dirw.p_ybrks[cj+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
       ray_mov(s, &ray);
       copy_vec(ray_vrt.ppos, ray.ppos, 3);
       copy_vec(ray_vrt.cdir, ray.cdir, 3);
       center_vec(ray_vrt.ppos, p_ellps[0].p_origin, 3);
       rotate_vec(ray_vrt.ppos, p_ellps[0].p_caxs, 0);
       rotate_vec(ray_vrt.cdir, p_ellps[0].p_caxs, 0);
       sdwf = !intrs_ellipse(ray_vrt.ppos, ray_vrt.cdir, s, &p_ellps[0]);
       if ( sdwf )
       {
         acm_grid(ray.ppos, psi, phi, ray.ppp_stks, &dirw, &skr, sclw, 1, 0);
       }
     }
   }

   write_grid (&dirw, 1, "plots/ellipse/set_04_test_02", "_out_spr",  
     &iop_w0, &btthdr, &psi);
   printf ("Result: %sCheck image in plots/ellipse after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   // Set 4, test 3: -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 04, Test 03 ------------------------------------------------\n");
   printf("Rays from 30 theta and 0 phi. \n\n%s", ANSI_RESET);

   acm_reset(&dirw);
   psi  = 30.0 * RAD;
   phi  = 0.0 * RAD;
   s    = 10.0;
   sclw = 1.0;

   for (ci = 1; ci < (dirw.nx - 1); ci++)
   {
     for (cj = 1; cj < (dirw.ny - 1); cj++)
     {
       src_pos[0] = (dirw.p_xbrks[ci] + dirw.p_xbrks[ci+1]) / 2.0;
       src_pos[1] = (dirw.p_ybrks[cj] + dirw.p_ybrks[cj+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
       ray_mov(s, &ray);
       copy_vec(ray_vrt.ppos, ray.ppos, 3);
       copy_vec(ray_vrt.cdir, ray.cdir, 3);
       center_vec(ray_vrt.ppos, p_ellps[0].p_origin, 3);
       rotate_vec(ray_vrt.ppos, p_ellps[0].p_caxs, 0);
       rotate_vec(ray_vrt.cdir, p_ellps[0].p_caxs, 0);
       sdwf = !intrs_ellipse(ray_vrt.ppos, ray_vrt.cdir, s, &p_ellps[0]);
       if ( sdwf )
       {
         acm_grid(ray.ppos, psi, phi, ray.ppp_stks, &dirw, &skr, sclw, 1, 0);
       }
     }
   }

   write_grid (&dirw, 1, "plots/ellipse/set_04_test_03", "_out_spr",  
     &iop_w0, &btthdr, &psi);
   printf ("Result: %sCheck image in plots/ellipse after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   // Set 4, test 4: -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 04, Test 04 ------------------------------------------------\n");
   printf("Rays from 30 theta and 90 phi. \n\n%s", ANSI_RESET);

   acm_reset(&dirw);
   psi  = 30.0 * RAD;
   phi  = 90.0 * RAD;
   s    = 10.0;
   sclw = 1.0;

   for (ci = 1; ci < (dirw.nx - 1); ci++)
   {
     for (cj = 1; cj < (dirw.ny - 1); cj++)
     {
       src_pos[0] = (dirw.p_xbrks[ci] + dirw.p_xbrks[ci+1]) / 2.0;
       src_pos[1] = (dirw.p_ybrks[cj] + dirw.p_ybrks[cj+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
       ray_mov(s, &ray);
       copy_vec(ray_vrt.ppos, ray.ppos, 3);
       copy_vec(ray_vrt.cdir, ray.cdir, 3);
       center_vec(ray_vrt.ppos, p_ellps[0].p_origin, 3);
       rotate_vec(ray_vrt.ppos, p_ellps[0].p_caxs, 0);
       rotate_vec(ray_vrt.cdir, p_ellps[0].p_caxs, 0);
       sdwf = !intrs_ellipse(ray_vrt.ppos, ray_vrt.cdir, s, &p_ellps[0]);
       if ( sdwf )
       {
         acm_grid(ray.ppos, psi, phi, ray.ppp_stks, &dirw, &skr, sclw, 1, 0);
       }
     }
   }

   write_grid (&dirw, 1, "plots/ellipse/set_04_test_04", "_out_spr",  
     &iop_w0, &btthdr, &psi);
   printf ("Result: %sCheck image in plots/ellipse after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   // Set 4, test 5: -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 04, Test 05 ------------------------------------------------\n");
   printf("Rays from 45 theta and 0 phi. \n\n%s", ANSI_RESET);

   acm_reset(&dirw);
   psi  = 45.0 * RAD;
   phi  = 0.0 * RAD;
   s    = 10.0;
   sclw = 1.0;

   for (ci = 1; ci < (dirw.nx - 1); ci++)
   {
     for (cj = 1; cj < (dirw.ny - 1); cj++)
     {
       src_pos[0] = (dirw.p_xbrks[ci] + dirw.p_xbrks[ci+1]) / 2.0;
       src_pos[1] = (dirw.p_ybrks[cj] + dirw.p_ybrks[cj+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
       ray_mov(s, &ray);
       copy_vec(ray_vrt.ppos, ray.ppos, 3);
       copy_vec(ray_vrt.cdir, ray.cdir, 3);
       center_vec(ray_vrt.ppos, p_ellps[0].p_origin, 3);
       rotate_vec(ray_vrt.ppos, p_ellps[0].p_caxs, 0);
       rotate_vec(ray_vrt.cdir, p_ellps[0].p_caxs, 0);
       sdwf = !intrs_ellipse(ray_vrt.ppos, ray_vrt.cdir, s, &p_ellps[0]);
       if ( sdwf )
       {
         acm_grid(ray.ppos, psi, phi, ray.ppp_stks, &dirw, &skr, sclw, 1, 0);
       }
     }
   }

   write_grid (&dirw, 1, "plots/ellipse/set_04_test_05", "_out_spr",  
     &iop_w0, &btthdr, &psi);
   printf ("Result: %sCheck image in plots/ellipse after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   // Set 4, test 6: -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 04, Test 06 ------------------------------------------------\n");
   printf("Rays from 30 theta and 180 phi. \n\n%s", ANSI_RESET);

   acm_reset(&dirw);
   psi  = 30.0 * RAD;
   phi  = 180.0 * RAD;
   s    = 10.0;
   sclw = 1.0;

   for (ci = 1; ci < (dirw.nx - 1); ci++)
   {
     for (cj = 1; cj < (dirw.ny - 1); cj++)
     {
       src_pos[0] = (dirw.p_xbrks[ci] + dirw.p_xbrks[ci+1]) / 2.0;
       src_pos[1] = (dirw.p_ybrks[cj] + dirw.p_ybrks[cj+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
       ray_mov(s, &ray);
       copy_vec(ray_vrt.ppos, ray.ppos, 3);
       copy_vec(ray_vrt.cdir, ray.cdir, 3);
       center_vec(ray_vrt.ppos, p_ellps[0].p_origin, 3);
       rotate_vec(ray_vrt.ppos, p_ellps[0].p_caxs, 0);
       rotate_vec(ray_vrt.cdir, p_ellps[0].p_caxs, 0);
       sdwf = !intrs_ellipse(ray_vrt.ppos, ray_vrt.cdir, s, &p_ellps[0]);
       if ( sdwf )
       {
         acm_grid(ray.ppos, psi, phi, ray.ppp_stks, &dirw, &skr, sclw, 1, 0);
       }
     }
   }

   write_grid (&dirw, 1, "plots/ellipse/set_04_test_06", "_out_spr",  
     &iop_w0, &btthdr, &psi);
   printf ("Result: %sCheck image in plots/ellipse after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);












   // Set 5: -------------------------------------------------------------------
   int n_set05 = 6;
   int *res_set05 = (int *) calloc(n_set05, sizeof(int));

   printf("\n%s", ANSI_BOLD);
   printf("Set 05 =========================================================\n");
   printf("Off-axis ellipse, 45º 0º imaging approach. \n\n%s", ANSI_RESET);

   str_nel = 1;
   p_ellps = (struct str_ellp*) calloc(str_nel, sizeof(struct str_ellp));

   origin[0] = 0.0;
   origin[1] = 0.0;
   origin[2] = 0.0;
   axis[0]   = 45.0 * RAD;
   axis[1]   = 0.0 * RAD;
   radius[0] = 0.25;
   radius[1] = 0.75;

   setup_ellp (&p_ellps[0], origin, axis, radius);
   print_ellp (&p_ellps[0]);

   // Set 5, test 1: -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 05, Test 01 ------------------------------------------------\n");
   printf("Rays from 0 theta and 0 phi. \n\n%s", ANSI_RESET);

   acm_reset(&dirw);
   psi  = 0.0 * RAD;
   phi  = 0.0 * RAD;
   s    = 10.0;
   sclw = 1.0;

   for (ci = 1; ci < (dirw.nx - 1); ci++)
   {
     for (cj = 1; cj < (dirw.ny - 1); cj++)
     {
       src_pos[0] = (dirw.p_xbrks[ci] + dirw.p_xbrks[ci+1]) / 2.0;
       src_pos[1] = (dirw.p_ybrks[cj] + dirw.p_ybrks[cj+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
       ray_mov(s, &ray);
       copy_vec(ray_vrt.ppos, ray.ppos, 3);
       copy_vec(ray_vrt.cdir, ray.cdir, 3);
       center_vec(ray_vrt.ppos, p_ellps[0].p_origin, 3);
       rotate_vec(ray_vrt.ppos, p_ellps[0].p_caxs, 0);
       rotate_vec(ray_vrt.cdir, p_ellps[0].p_caxs, 0);
       sdwf = !intrs_ellipse(ray_vrt.ppos, ray_vrt.cdir, s, &p_ellps[0]);
       if ( sdwf )
       {
         acm_grid(ray.ppos, psi, phi, ray.ppp_stks, &dirw, &skr, sclw, 1, 0);
       }
     }
   }

   write_grid (&dirw, 1, "plots/ellipse/set_05_test_01", "_out_spr",  
     &iop_w0, &btthdr, &psi);
   printf ("Result: %sCheck image in plots/ellipse after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   // Set 5, test 2: -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 05, Test 02 ------------------------------------------------\n");
   printf("Rays from 15 theta and 0 phi. \n\n%s", ANSI_RESET);

   acm_reset(&dirw);
   psi  = 15.0 * RAD;
   phi  =  0.0 * RAD;
   s    = 10.0;
   sclw = 1.0;

   for (ci = 1; ci < (dirw.nx - 1); ci++)
   {
     for (cj = 1; cj < (dirw.ny - 1); cj++)
     {
       src_pos[0] = (dirw.p_xbrks[ci] + dirw.p_xbrks[ci+1]) / 2.0;
       src_pos[1] = (dirw.p_ybrks[cj] + dirw.p_ybrks[cj+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
       ray_mov(s, &ray);
       copy_vec(ray_vrt.ppos, ray.ppos, 3);
       copy_vec(ray_vrt.cdir, ray.cdir, 3);
       center_vec(ray_vrt.ppos, p_ellps[0].p_origin, 3);
       rotate_vec(ray_vrt.ppos, p_ellps[0].p_caxs, 0);
       rotate_vec(ray_vrt.cdir, p_ellps[0].p_caxs, 0);
       sdwf = !intrs_ellipse(ray_vrt.ppos, ray_vrt.cdir, s, &p_ellps[0]);
       if ( sdwf )
       {
         acm_grid(ray.ppos, psi, phi, ray.ppp_stks, &dirw, &skr, sclw, 1, 0);
       }
     }
   }

   write_grid (&dirw, 1, "plots/ellipse/set_05_test_02", "_out_spr",  
     &iop_w0, &btthdr, &psi);
   printf ("Result: %sCheck image in plots/ellipse after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   // Set 5, test 3: -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 05, Test 03 ------------------------------------------------\n");
   printf("Rays from 30 theta and 0 phi. \n\n%s", ANSI_RESET);

   acm_reset(&dirw);
   psi  = 30.0 * RAD;
   phi  = 0.0 * RAD;
   s    = 10.0;
   sclw = 1.0;

   for (ci = 1; ci < (dirw.nx - 1); ci++)
   {
     for (cj = 1; cj < (dirw.ny - 1); cj++)
     {
       src_pos[0] = (dirw.p_xbrks[ci] + dirw.p_xbrks[ci+1]) / 2.0;
       src_pos[1] = (dirw.p_ybrks[cj] + dirw.p_ybrks[cj+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
       ray_mov(s, &ray);
       copy_vec(ray_vrt.ppos, ray.ppos, 3);
       copy_vec(ray_vrt.cdir, ray.cdir, 3);
       center_vec(ray_vrt.ppos, p_ellps[0].p_origin, 3);
       rotate_vec(ray_vrt.ppos, p_ellps[0].p_caxs, 0);
       rotate_vec(ray_vrt.cdir, p_ellps[0].p_caxs, 0);
       sdwf = !intrs_ellipse(ray_vrt.ppos, ray_vrt.cdir, s, &p_ellps[0]);
       if ( sdwf )
       {
         acm_grid(ray.ppos, psi, phi, ray.ppp_stks, &dirw, &skr, sclw, 1, 0);
       }
     }
   }

   write_grid (&dirw, 1, "plots/ellipse/set_05_test_03", "_out_spr",  
     &iop_w0, &btthdr, &psi);
   printf ("Result: %sCheck image in plots/ellipse after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   // Set 5, test 4: -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 05, Test 04 ------------------------------------------------\n");
   printf("Rays from 30 theta and 90 phi. \n\n%s", ANSI_RESET);

   acm_reset(&dirw);
   psi  = 30.0 * RAD;
   phi  = 90.0 * RAD;
   s    = 10.0;
   sclw = 1.0;

   for (ci = 1; ci < (dirw.nx - 1); ci++)
   {
     for (cj = 1; cj < (dirw.ny - 1); cj++)
     {
       src_pos[0] = (dirw.p_xbrks[ci] + dirw.p_xbrks[ci+1]) / 2.0;
       src_pos[1] = (dirw.p_ybrks[cj] + dirw.p_ybrks[cj+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
       ray_mov(s, &ray);
       copy_vec(ray_vrt.ppos, ray.ppos, 3);
       copy_vec(ray_vrt.cdir, ray.cdir, 3);
       center_vec(ray_vrt.ppos, p_ellps[0].p_origin, 3);
       rotate_vec(ray_vrt.ppos, p_ellps[0].p_caxs, 0);
       rotate_vec(ray_vrt.cdir, p_ellps[0].p_caxs, 0);
       sdwf = !intrs_ellipse(ray_vrt.ppos, ray_vrt.cdir, s, &p_ellps[0]);
       if ( sdwf )
       {
         acm_grid(ray.ppos, psi, phi, ray.ppp_stks, &dirw, &skr, sclw, 1, 0);
       }
     }
   }

   write_grid (&dirw, 1, "plots/ellipse/set_05_test_04", "_out_spr",  
     &iop_w0, &btthdr, &psi);
   printf ("Result: %sCheck image in plots/ellipse after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   // Set 5, test 5: -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 05, Test 05 ------------------------------------------------\n");
   printf("Rays from 30 theta and 180 phi. \n\n%s", ANSI_RESET);

   acm_reset(&dirw);
   psi  = 30.0 * RAD;
   phi  = 180.0 * RAD;
   s    = 10.0;
   sclw = 1.0;

   for (ci = 1; ci < (dirw.nx - 1); ci++)
   {
     for (cj = 1; cj < (dirw.ny - 1); cj++)
     {
       src_pos[0] = (dirw.p_xbrks[ci] + dirw.p_xbrks[ci+1]) / 2.0;
       src_pos[1] = (dirw.p_ybrks[cj] + dirw.p_ybrks[cj+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
       ray_mov(s, &ray);
       copy_vec(ray_vrt.ppos, ray.ppos, 3);
       copy_vec(ray_vrt.cdir, ray.cdir, 3);
       center_vec(ray_vrt.ppos, p_ellps[0].p_origin, 3);
       rotate_vec(ray_vrt.ppos, p_ellps[0].p_caxs, 0);
       rotate_vec(ray_vrt.cdir, p_ellps[0].p_caxs, 0);
       sdwf = !intrs_ellipse(ray_vrt.ppos, ray_vrt.cdir, s, &p_ellps[0]);
       if ( sdwf )
       {
         acm_grid(ray.ppos, psi, phi, ray.ppp_stks, &dirw, &skr, sclw, 1, 0);
       }
     }
   }

   write_grid (&dirw, 1, "plots/ellipse/set_05_test_05", "_out_spr",  
     &iop_w0, &btthdr, &psi);
   printf ("Result: %sCheck image in plots/ellipse after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);

   // Set 5, test 6: -----------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Set 05, Test 06 ------------------------------------------------\n");
   printf("Rays from 45 theta and 0 phi. \n\n%s", ANSI_RESET);

   acm_reset(&dirw);
   psi  = 45.0 * RAD;
   phi  = 0.0 * RAD;
   s    = 10.0;
   sclw = 1.0;

   for (ci = 1; ci < (dirw.nx - 1); ci++)
   {
     for (cj = 1; cj < (dirw.ny - 1); cj++)
     {
       src_pos[0] = (dirw.p_xbrks[ci] + dirw.p_xbrks[ci+1]) / 2.0;
       src_pos[1] = (dirw.p_ybrks[cj] + dirw.p_ybrks[cj+1]) / 2.0;
       src_pos[2] = -1.0;
       ray_ini(src_pos, psi, phi, sim_nw0, sim_nbr, &ray);
       ray_mov(s, &ray);
       copy_vec(ray_vrt.ppos, ray.ppos, 3);
       copy_vec(ray_vrt.cdir, ray.cdir, 3);
       center_vec(ray_vrt.ppos, p_ellps[0].p_origin, 3);
       rotate_vec(ray_vrt.ppos, p_ellps[0].p_caxs, 0);
       rotate_vec(ray_vrt.cdir, p_ellps[0].p_caxs, 0);
       sdwf = !intrs_ellipse(ray_vrt.ppos, ray_vrt.cdir, s, &p_ellps[0]);
       if ( sdwf )
       {
         acm_grid(ray.ppos, psi, phi, ray.ppp_stks, &dirw, &skr, sclw, 1, 0);
       }
     }
   }

   write_grid (&dirw, 1, "plots/ellipse/set_05_test_06", "_out_spr",  
     &iop_w0, &btthdr, &psi);
   printf ("Result: %sCheck image in plots/ellipse after R script%s\n\n", 
     ANSI_YELLOW, ANSI_RESET);


   // ==========================================================================

   printf("Summary ========================================================\n");
   printf("  Set 01: ------------------------------------------------------\n");
   for (i = 0; i < n_set01; i++)
   {
     printf("    Test %02d: %s %s %s\n", i+1, 
       res_set01[i] ? ANSI_BLUE : ANSI_RED, 
       res_set01[i] ? "PASSED!" : "FAILED!", ANSI_RESET);
   }
   printf("  Set 02: ------------------------------------------------------\n");
   for (i = 0; i < n_set02; i++)
   {
     printf("    Test %02d: %s %s %s\n", i+1, 
       res_set02[i] ? ANSI_BLUE : ANSI_RED, 
       res_set02[i] ? "PASSED!" : "FAILED!", ANSI_RESET);
   }
   printf("  Set 03: ------------------------------------------------------\n");
   for (ci = 0; ci < n_set03; ci++)
   {
     printf("    Test %02d: %s %s %s\n", ci+1, ANSI_YELLOW, 
       "Check image in plots/ellipse after R script", ANSI_RESET);
   }
   printf("  Set 04: ------------------------------------------------------\n");
   for (ci = 0; ci < n_set04; ci++)
   {
     printf("    Test %02d: %s %s %s\n", ci+1, ANSI_YELLOW, 
       "Check image in plots/ellipse after R script", ANSI_RESET);
   }
   printf("  Set 05: ------------------------------------------------------\n");
   for (ci = 0; ci < n_set05; ci++)
   {
     printf("    Test %02d: %s %s %s\n", ci+1, ANSI_YELLOW, 
       "Check image in plots/ellipse after R script", ANSI_RESET);
   }
   printf("\nRunning R script test_ellipse.R:\n");
   system("Rscript test_ellipse.R");
   printf("Done!\n");
 }
