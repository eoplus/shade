
 #include <stdio.h>
 #include <stdlib.h>

 #include "../../../src/constants.h"
 #include "../../../src/structures.h"

 void main (void)
 {
   FILE   *p_fi;
   long int fpos;

   int    i;
   int    str_nel;
   struct str_ellp *p_ellps;

   double origin[3];
   double axis[2];
   double radius[2];

   str_nel = 1;
   p_ellps = (struct str_ellp*) calloc(str_nel, sizeof(struct str_ellp));

   // Test 1: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 01 ========================================================\n");
   printf("Z-axis oriented circle. \n\n%s", ANSI_RESET);

   origin[0] = 0.0;
   origin[1] = 0.0;
   origin[2] = 0.0;
   axis[0]   = 0.0 * RAD;
   axis[1]   = 0.0 * RAD;
   radius[0] = 1.0;
   radius[1] = 1.0;

   setup_ellp (&p_ellps[0], origin, axis, radius);
   print_ellp (&p_ellps[0]);

   // Test 2: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 02 ========================================================\n");
   printf("Circle rotated by 30º. \n\n%s", ANSI_RESET);

   origin[0] = 0.0;
   origin[1] = 0.0;
   origin[2] = 0.0;
   axis[0]   = 30.0 * RAD;
   axis[1]   = 0.0 * RAD;
   radius[0] = 1.0;
   radius[1] = 1.0;

   setup_ellp (&p_ellps[0], origin, axis, radius);
   print_ellp (&p_ellps[0]);

   // Test 3: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 03 ========================================================\n");
   printf("Circle rotated by 60º. \n\n%s", ANSI_RESET);

   origin[0] = 0.0;
   origin[1] = 0.0;
   origin[2] = 0.0;
   axis[0]   = 60.0 * RAD;
   axis[1]   = 0.0 * RAD;
   radius[0] = 1.0;
   radius[1] = 1.0;

   setup_ellp (&p_ellps[0], origin, axis, radius);
   print_ellp (&p_ellps[0]);

   // Test 4: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 04 ========================================================\n");
   printf("Circle rotated by 120º. \n\n%s", ANSI_RESET);

   origin[0] = 0.0;
   origin[1] = 0.0;
   origin[2] = 0.0;
   axis[0]   = 120.0 * RAD;
   axis[1]   = 0.0 * RAD;
   radius[0] = 1.0;
   radius[1] = 1.0;

   setup_ellp (&p_ellps[0], origin, axis, radius);
   print_ellp (&p_ellps[0]);

   // Test 5: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 05 ========================================================\n");
   printf("Circle rotated by 120º polar and 30º azimuth. \n\n%s", ANSI_RESET);

   origin[0] = 0.0;
   origin[1] = 0.0;
   origin[2] = 0.0;
   axis[0]   = 120.0 * RAD;
   axis[1]   = 30.0 * RAD;
   radius[0] = 1.0;
   radius[1] = 1.0;

   setup_ellp (&p_ellps[0], origin, axis, radius);
   print_ellp (&p_ellps[0]);

   // Test 6: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 06 ========================================================\n");
   printf("Z-axis oriented ellipse. \n\n%s", ANSI_RESET);

   origin[0] = 0.0;
   origin[1] = 0.0;
   origin[2] = 0.0;
   axis[0]   = 0.0 * RAD;
   axis[1]   = 0.0 * RAD;
   radius[0] = 1.0;
   radius[1] = 2.0;

   setup_ellp (&p_ellps[0], origin, axis, radius);
   print_ellp (&p_ellps[0]);

   // Test 7: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 07 ========================================================\n");
   printf("Ellipse rotated by 120º polar and 30º azimuth. \n\n%s", ANSI_RESET);

   origin[0] = 0.0;
   origin[1] = 0.0;
   origin[2] = 0.0;
   axis[0]   = 120.0 * RAD;
   axis[1]   = 30.0 * RAD;
   radius[0] = 1.0;
   radius[1] = 2.0;

   setup_ellp (&p_ellps[0], origin, axis, radius);
   print_ellp (&p_ellps[0]);

   // Test 8: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 08 ========================================================\n");
   printf("Same as test 07, but a summary. \n\n%s", ANSI_RESET);

   origin[0] = 0.0;
   origin[1] = 0.0;
   origin[2] = 0.0;
   axis[0]   = 120.0 * RAD;
   axis[1]   = 30.0 * RAD;
   radius[0] = 1.0;
   radius[1] = 2.0;

   setup_ellp (&p_ellps[0], origin, axis, radius);
   print_ellps (p_ellps, str_nel);

   // Test 9: ------------------------------------------------------------------
   printf("\n%s", ANSI_BOLD);
   printf("Test 09 ========================================================\n");
   printf("Read ellipses from disk. \n\n%s", ANSI_RESET);

   str_nel = 3;
   p_ellps = (struct str_ellp*) calloc(str_nel, sizeof(struct str_ellp));

   p_fi = fopen("test_ellipse.txt", "r");
   if(p_fi == NULL)
   {
     printf("\nERROR: Failed to open input file test_ellipse.txt\n\n");
     exit(-1);
   }
   fpos = ftell(p_fi);

     read_ellp (p_fi, &fpos, origin, axis, radius);   
     setup_ellp (&p_ellps[0], origin, axis, radius);
     read_ellp (p_fi, &fpos, origin, axis, radius);   
     setup_ellp (&p_ellps[1], origin, axis, radius);
     read_ellp (p_fi, &fpos, origin, axis, radius);   
     setup_ellp (&p_ellps[2], origin, axis, radius);
/*
   for (i = 0; i < str_nel; i++)
   {
     read_ellp (p_fi, &fpos, origin, axis, radius);   
     setup_ellp (&p_ellps[i], origin, axis, radius);
   }
*/   fclose(p_fi);

   print_ellps(p_ellps, str_nel);
   printf("\n");

   for (i = 0; i < str_nel; i++)
   {
     print_ellp (&p_ellps[i]);
     printf("\n");
   }
}

