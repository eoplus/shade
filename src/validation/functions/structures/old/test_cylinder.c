
 #include <stdio.h>
 #include <stdlib.h>

 #include "../../../src/constants.h"
 #include "../../../src/structures.h"

 void main (void)
 {
   FILE   *p_fi;
   long int fpos;
   int    str_ncl;
   struct str_cylnd *p_cylnd;
   double base[3];
   double axis[2];
   double axis0[2];
   double axis1[2];
   double radius;
   double height;
   int    closed[2];

   str_ncl = 1;
   p_cylnd = (struct str_cylnd*) calloc(str_ncl, sizeof(struct str_cylnd));

   // Test 1: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 01 ========================================================\n");
   printf("Cylinder as in our PONDER study. \n\n%s", ANSI_RESET);

   base[0]   = 0.0;
   base[1]   = 0.0;
   base[2]   = -0.09;
   axis[0]   = 0.0 * RAD;
   axis[1]   = 0.0 * RAD;
   axis0[0]  = 0.0 * RAD;
   axis0[1]  = 0.0 * RAD;
   axis1[0]  = 0.0 * RAD;
   axis1[1]  = 0.0 * RAD;
   radius    = 0.0125;
   height    = 0.11;
   closed[0] = 1;
   closed[1] = 0;

   setup_cylnd (&p_cylnd[0], base, axis, axis0, axis1, closed, radius, height);
   print_cylnd (&p_cylnd[0]);

   // Test 2: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 02 ========================================================\n");
   printf("Cylinder inclined by 10º. \n\n%s", ANSI_RESET);

   base[0]   = 0.0;
   base[1]   = 0.0;
   base[2]   = -0.09;
   axis[0]   = 10.0 * RAD;
   axis[1]   = 0.0 * RAD;
   axis0[0]  = 0.0 * RAD;
   axis0[1]  = 0.0 * RAD;
   axis1[0]  = 0.0 * RAD;
   axis1[1]  = 0.0 * RAD;
   radius    = 0.0125;
   height    = 0.11;
   closed[0] = 1;
   closed[1] = 0;

   setup_cylnd (&p_cylnd[0], base, axis, axis0, axis1, closed, radius, height);
   print_cylnd (&p_cylnd[0]);

   // Test 3: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 03 ========================================================\n");
   printf("Cylinder inclined by 30º. \n\n%s", ANSI_RESET);

   base[0]   = 0.0;
   base[1]   = 0.0;
   base[2]   = -0.09;
   axis[0]   = 30.0 * RAD;
   axis[1]   = 0.0 * RAD;
   axis0[0]  = 0.0 * RAD;
   axis0[1]  = 0.0 * RAD;
   axis1[0]  = 0.0 * RAD;
   axis1[1]  = 0.0 * RAD;
   radius    = 0.0125;
   height    = 0.11;
   closed[0] = 1;
   closed[1] = 0;

   setup_cylnd (&p_cylnd[0], base, axis, axis0, axis1, closed, radius, height);
   print_cylnd (&p_cylnd[0]);

   // Test 4: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 04 ========================================================\n");
   printf("Cylinder inclined by 60º. \n\n%s", ANSI_RESET);

   base[0]   = 0.0;
   base[1]   = 0.0;
   base[2]   = -0.09;
   axis[0]   = 60.0 * RAD;
   axis[1]   = 0.0 * RAD;
   axis0[0]  = 0.0 * RAD;
   axis0[1]  = 0.0 * RAD;
   axis1[0]  = 0.0 * RAD;
   axis1[1]  = 0.0 * RAD;
   radius    = 0.0125;
   height    = 0.11;
   closed[0] = 1;
   closed[1] = 0;

   setup_cylnd (&p_cylnd[0], base, axis, axis0, axis1, closed, radius, height);
   print_cylnd (&p_cylnd[0]);

   // Test 5: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 05 ========================================================\n");
   printf("Oblique cylinder inclined by 10º. \n\n%s", ANSI_RESET);

   base[0]   = 0.0;
   base[1]   = 0.0;
   base[2]   = -0.09;
   axis[0]   = 10.0 * RAD;
   axis[1]   = 0.0 * RAD;
   axis0[0]  = 10.0 * RAD;
   axis0[1]  = 180.0 * RAD;
   axis1[0]  = 10.0 * RAD;
   axis1[1]  = 180.0 * RAD;
   radius    = 0.0125;
   height    = 0.11;
   closed[0] = 1;
   closed[1] = 0;

   setup_cylnd (&p_cylnd[0], base, axis, axis0, axis1, closed, radius, height);
   print_cylnd (&p_cylnd[0]);

   // Test 6: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 06 ========================================================\n");
   printf("Oblique cylinder inclined by 30º. \n\n%s", ANSI_RESET);

   base[0]   = 0.0;
   base[1]   = 0.0;
   base[2]   = -0.09;
   axis[0]   = 30.0 * RAD;
   axis[1]   = 0.0 * RAD;
   axis0[0]  = 30.0 * RAD;
   axis0[1]  = 180.0 * RAD;
   axis1[0]  = 30.0 * RAD;
   axis1[1]  = 180.0 * RAD;
   radius    = 0.0125;
   height    = 0.11;
   closed[0] = 1;
   closed[1] = 0;

   setup_cylnd (&p_cylnd[0], base, axis, axis0, axis1, closed, radius, height);
   print_cylnd (&p_cylnd[0]);

   // Test 7: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 07 ========================================================\n");
   printf("Oblique cylinder inclined by 60º. \n\n%s", ANSI_RESET);

   base[0]   = 0.0;
   base[1]   = 0.0;
   base[2]   = -0.09;
   axis[0]   = 60.0 * RAD;
   axis[1]   = 0.0 * RAD;
   axis0[0]  = 60.0 * RAD;
   axis0[1]  = 180.0 * RAD;
   axis1[0]  = 60.0 * RAD;
   axis1[1]  = 180.0 * RAD;
   radius    = 0.0125;
   height    = 0.11;
   closed[0] = 1;
   closed[1] = 0;

   setup_cylnd (&p_cylnd[0], base, axis, axis0, axis1, closed, radius, height);
   print_cylnd (&p_cylnd[0]);

   // Test 8: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 08 ========================================================\n");
   printf("Oblique at the top only, and inclined by 60º. \n\n%s", ANSI_RESET);

   base[0]   = 0.0;
   base[1]   = 0.0;
   base[2]   = -0.09;
   axis[0]   = 60.0 * RAD;
   axis[1]   = 0.0 * RAD;
   axis0[0]  = 0.0 * RAD;
   axis0[1]  = 0.0 * RAD;
   axis1[0]  = 60.0 * RAD;
   axis1[1]  = 180.0 * RAD;
   radius    = 0.0125;
   height    = 0.11;
   closed[0] = 1;
   closed[1] = 0;

   setup_cylnd (&p_cylnd[0], base, axis, axis0, axis1, closed, radius, height);
   print_cylnd (&p_cylnd[0]);

   // Test 9: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 09 ========================================================\n");
   printf("The same as test 08, but print summary. \n\n%s", ANSI_RESET);
   print_cylnds(p_cylnd, str_ncl);

   // Test 10: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 10 ========================================================\n");
   printf("Read cylinder from disk. \n\n%s", ANSI_RESET);

   str_ncl = 2;
   p_cylnd = (struct str_cylnd*) calloc(str_ncl, sizeof(struct str_cylnd));

   p_fi = fopen("test_cylinder.txt", "r");
   if(p_fi == NULL)
   {
     printf("\nERROR: Failed to open input file test_cylinder.txt\n\n");
     exit(-1);
   }
   fpos = ftell(p_fi);
   read_cylnd (p_fi, &fpos, base, axis, axis0, axis1, closed, &radius, &height);   
   setup_cylnd (&p_cylnd[0], base, axis, axis0, axis1, closed, radius, height);
   read_cylnd (p_fi, &fpos, base, axis, axis0, axis1, closed, &radius, &height);   
   setup_cylnd (&p_cylnd[1], base, axis, axis0, axis1, closed, radius, height);
   fclose(p_fi);

   print_cylnds(p_cylnd, str_ncl);
   printf("\n");
   print_cylnd (&p_cylnd[0]);
   printf("\n");
   print_cylnd (&p_cylnd[1]);
   printf("\n");
 }

