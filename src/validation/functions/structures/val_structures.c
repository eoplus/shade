
/*******************************************************************************
 val_structures.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 This is a driver program to test if the structures function work as expected.

 Compile instruction:
 gcc val_structures.c ../../../src/geometry.c ../../../src/structures.c ../../../src/rotation.c ../../../src/memory.c -lm -lgsl -O3 -o val_structures.o

 ./val_structures.o
*******************************************************************************/

 #include <stdio.h>
 #include <math.h>
 #include <gsl/gsl_blas.h>
 #include <gsl/gsl_rng.h>
 #include <sys/time.h>

 #include "../../../src/config.h"
 #include "../../../src/aux.h"
 #include "../../../src/constants.h"
 #include "../../../src/memory.h"
 #include "../../../src/geometry.h"
 #include "../../../src/structures.h"

 void main ( void )
 {
   str_rect *rect;
   str_rect **rects;

   FILE      *fi;
   long int fpos;
   double origin[3];
   double axis[2];
   double lengths[3];
   double alpha;

  /* Tests for rectangles *****************************************************/
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
     str_rect_read(fi, &fpos, origin, axis, &alpha, lengths);   
     str_rect_setup(rects[i], origin, axis, alpha, lengths);   
   }
   printf("\n");
   str_rects_printf(rects, str_nrt);
   printf("\n");
   str_rect_printf(rects[0]);
   printf("\n");
   str_rect_printf(rects[1]);
   printf("\n");

  /* Tests for ellipses *******************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("Ellipses ======================================================\n%s",
     ANSI_RESET);

   int str_nel = 2;
   str_ellp **ellps;
   ellps = (str_ellp**) calloc(str_nel, sizeof(str_ellp*));
   double radius[2];

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
     str_ellp_read(fi, &fpos, origin, axis, &alpha, radius);
     str_ellp_setup(ellps[i], origin, axis, alpha, radius);   
   }
   printf("\n");
   str_ellps_printf(ellps, str_nel);
   printf("\n");
   str_ellp_printf(ellps[0]);
   printf("\n");
   str_ellp_printf(ellps[1]);
   printf("\n");

  /* Tests for cuboids ********************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("Cuboids =======================================================\n%s",
     ANSI_RESET);

   double base_s[3];
   double top_s[3];
   int    closed[2];
   int str_ncb = 2;
   str_cubd **cubds;
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
   str_cubds_printf(cubds, str_ncb);
   printf("\n");
   str_cubd_printf(cubds[0]);
   printf("\n");
   str_cubd_printf(cubds[1]);
   printf("\n");

  /* Tests for cylinders ******************************************************/
   printf("\n%s", ANSI_BOLD);
   printf("Cylinders =====================================================\n%s",
     ANSI_RESET);

   double height;
   int str_ncl = 2;
   str_cyln **cylns;
   cylns = (str_cyln**) calloc(str_ncl, sizeof(str_cyln*));

   fi = fopen("input_cylinder.txt", "r");
   if(fi == NULL)
   {
     printf("\nERROR: Failed to open input file input_cylinder.txt\n\n");
     exit(-1);
   }
   fpos = ftell(fi);

   for (size_t i = 0; i < str_ncl; i++)
   {
     cylns[i] = str_cyln_alloc();
     str_cyln_read(fi, &fpos, origin, axis, &alpha, &radius[0], &height, base_s,
       top_s, closed);
     str_cyln_setup(cylns[i], origin, axis, alpha, radius[0], height, base_s,
       top_s, closed);   
   }
   printf("\n");
   str_cylns_printf(cylns, str_ncb);
   printf("\n");
   str_cyln_printf(cylns[0]);
   printf("\n");
   str_cyln_printf(cylns[1]);
   printf("\n");
 }

