
/*******************************************************************************
 structures_cylinder.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 Functions to read, write and setup the internal representation of a positive
 finite circular cylinder. The definition of the cone struct is detailed in the
 structures_cone.h header file.

 The functions are:
 str_cone_alloc   - Allocate structure;
 str_cone_free    - Release memory of structure;
 str_cone_read    - Read parameters from input file;
 str_cone_setup   - Initialize structure with parameter values;
 str_cone_printf  - Formatted detailed print of a structure;
 str_cones_fprintf - Formatted summary print of an array of structures;

*******************************************************************************/

 #include <stdio.h>			// FILE
 #include <stddef.h>			// size_t
 #include <stdlib.h>			// exit, malloc, free
 #include <string.h>			// strcat, strncpy
 #include <math.h>			// M_PI, M_PI_2, cos, sin

 #include "config.h"			// STRMXLEN, TOLERANCE
 #include "aux.h"			// ABS, SIGN
 #include "constants.h"			// RAD, DEG, K_2PI
 #include "memory.h"			// calloc_2d, free_2d
 #include "geometry.h"			// sph_to_cos_unit, mat_transpose, mat_fprintf
 #include "rotation.h"			// rot_mat_ZYZ
 #include "structures.h"		// str_ellp
 #include "structures_cylinder.h"	// str_cylinder

/* Allocate cylinder struct: ***************************************************

 str_cyln_alloc

 Allocate memory for a cylinder structure (and its components).

 INPUT:
 None.

 OUTPUT:
 Pointer to a cylinder structure.

*******************************************************************************/

 struct str_cyln * 
 str_cyln_alloc
 ( void )
 {
   struct str_cyln * cyln = 
     (struct str_cyln *) malloc( sizeof(struct str_cyln) );
   cyln->base = str_ellp_alloc();
   cyln->top = str_ellp_alloc();
   cyln->M = calloc_2d(3, 3, "cyln->M in str_cyln_alloc (structures.h)");
   cyln->alp = 0.0;
   cyln->closed_f[0] = 1;
   cyln->closed_f[1] = 1;
   cyln->rotate_f[0] = 0;
   cyln->rotate_f[1] = 0;
   cyln->rotate_f[2] = 0;
   cyln->rotate_f[3] = 0;

   return cyln;
 }

/* Deallocate cylinder struct: *************************************************

 str_cyln_free

 Release memory of an cylinder structure.

 INPUT:
 cyln - Pointer to an cylinder structure pointer.

 OUTPUT:
 None. Releases memory block and set the cylinder pointer to NULL.

*******************************************************************************/

 void 
 str_cyln_free
 ( struct str_cyln ** cyln )
 {
   str_ellp_free( &(*cyln)->base );
   str_ellp_free( &(*cyln)->top );
   free_2d( 3, &(*cyln)->M );
   free( *cyln );
   *cyln = NULL;
 }

/* Read cylinder parameters from input file: ***********************************

 str_cyln_read

 Reads an cylinder from a formatted input file.

 INPUT:
 fi      - File stream pointer;
 fpos    - File stream position pointer;
 origin  - Pointer to array to receive center of the cylinder's base;
 axis    - Pointer to array to receive orientation of the cylinder's normal;
 alpha   - Pointer to azimuthal angle around the cylinder's normal;
 radius  - Pointer to double to receive the radius;
 height  - Pointer to double to receive the height;
 base_s  - Pointer to array to receive the spherical coordinates of the base
           openings relative to cylinder's axis;
 top_s   - Pointer to array to receive the spherical coordinates of the top
           openings relative to cylinder's axis;
 closed  - Pointer to array to receive flag if base and top are open or
           closed;

 OUTPUT:
 None. Updates the values of origin, axis, alpha, radius, height, base_s, top_s,
 and closed.

*******************************************************************************/

 void
 str_cyln_read
 (
   FILE * fi, 
   long int * fpos,
   double * origin,
   double * axis,
   double * alpha,
   double * radius,
   double * height,
   double * s_base,
   double * s_top,
   int * closed
 )
 {
   fseek(fi, *fpos, SEEK_SET);
   fscanf(fi, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d", 
     &origin[0], 
     &origin[1],
     &origin[2], 
     &axis[0],
     &axis[1],
     alpha,
     radius,
     height,
     &s_base[0],
     &s_base[1],
     &s_top[0],
     &s_top[1],
     &closed[0],
     &closed[1]);
   *fpos = ftell(fi);

   axis[0] *= RAD;
   axis[1] *= RAD;
   *alpha *= RAD;
   s_base[0] *= RAD;
   s_base[1] *= RAD;
   s_top[0] *= RAD;
   s_top[1] *= RAD;

   if ( (axis[0] < 0.0) || (axis[0] > M_PI) || 
        (axis[1] < 0.0) || (axis[1] > K_2PI) )
   {
     printf("\nERROR: The orientation axis of the cylinder must be between "
       "0º and 180º for the polar angle and 0º and 360º for the azimuthal "
       "angle.\n\n");
     exit(-1);
   }
   if ( ((*alpha) < 0.0) || ((*alpha) > K_2PI) )
   {
     printf("\nERROR: The azimuthal rotation angle around the normal must be "
       "between 0º and 360º.\n\n");
     exit(-1);
   }
   if ( (s_base[0] < 0.0) || (s_base[0] > (60.0 * RAD)) || 
        ( (s_base[1] != 0.0) && (s_base[1] != K_2PI) ) )
   {
     printf("\nERROR: The orientation axis of the cylinder's base must be "
       "between 0º and 60º for the polar angle and 0º or 180º for the "
       "azimuthal angle.\n\n");
     exit(-1);
   }
   if ( (s_top[0] < 0.0) || (s_top[0] > (60.0 * RAD)) || 
        ( (s_top[1] > TOLERANCE) && (ABS(s_top[1] - M_PI) > TOLERANCE) ) )
   {
      printf("\nERROR: The orientation axis of the cylinder's top must be "
       "between 0º and 60º for the polar angle and 0º or 180º for the "
       "azimuthal angle.\n\n");
     exit(-1);
   }
   if ( ((*radius) < TOLERANCE) || 
        ((*height) < TOLERANCE) )
   {
     printf("\nERROR: The radius and height must be greater than 0.\n\n");
     exit(-1);
   }
 }

/* Setup cylinder: *************************************************************

 str_cyln_setup

 Setup the cylinder parameter structure. The details of the cylinder structure
 are presented in the structures.h header file.

 INPUT:
 cyln    - Pointer to the cylinder structure;
 origin  - Pointer to array containing the center of the cylinder's base;
 axis    - Pointer to array containing the orientation of the cylinder's normal;
 alpha   - Pointer to azimuthal angle around the cylinder's normal;
 radius  - Pointer to array containing the radius of the cylinder;
 height  - Pointer to array containing the height of the cylinder along its
           axis;
 base_s  - Pointer to array containing the orientation of the cylinder's base
           normal relative to its axis;
 top_s   - Pointer to array containing the orientation of the cylinder's top
           normal relative to its axis;
 closed  - Pointer to array containing the flag if openings are open or close.

 OUTPUT:
 None. Updates the cylinder structure pointed by cyln.

*******************************************************************************/

 void 
 str_cyln_setup
 (
   struct str_cyln * cyln,
   double * origin,
   double const * axis,
   double const alpha,
   double const radius,
   double const height,
   double const * s_base,
   double const * s_top,
   int const * closed
 )
 {
   for (size_t i = 0; i < 2; i++)
   {
     cyln->o[i]  = origin[i];
     cyln->s[i]  = axis[i];
     cyln->closed_f[i] = closed[i];
     cyln->s_base[i] = s_base[i];
     cyln->s_top[i] = s_top[i];
   }
   cyln->o[2] = origin[2];
   cyln->s[2] = 1.0;
   cyln->r   = radius;
   cyln->rsq = radius * radius;
   cyln->h   = height;
   cyln->alp = alpha;

   sph_to_cos_unit(cyln->u, cyln->s);
   if ( (1.0 - ABS(cyln->u[2]) < TOLERANCE) && 
        (axis[1] < TOLERANCE) && (alpha < TOLERANCE) )
   {
     cyln->u[0] = 0.0;
     cyln->u[1] = 0.0;
     cyln->u[2] = 1.0;
     cyln->M[0][0] = SIGN(cyln->u[2]);
     cyln->M[1][1] = SIGN(cyln->u[2]);
     cyln->M[2][2] = SIGN(cyln->u[2]);
     cyln->rotate_f[0] = 0;
   }
   else
   {
     double **R = calloc_2d(3, 3, "R in str_cyln_setup (structures.c)");
     rot_mat_ZYZ (R, cyln->u, cyln->s, alpha);
     mat_transpose(cyln->M, (double const **) R, 3, 3);
     free_2d(3, &R);
     cyln->rotate_f[0] = 1;
   }

   // Set openings:
   double open_r[2];
   open_r[0] = radius / cos(s_base[0]);
   open_r[1] = radius;
   cyln->dhmx[0] = tan(s_base[0]); //-sin(s_base[0]);
   str_ellp_setup(cyln->base, origin, s_base, 0.0, open_r);
   cyln->rotate_f[1] = (s_base[0] < TOLERANCE)? 0 : 1;
   for (size_t i = 0; i < 3; i++)
   {
     origin[i] += cyln->u[i] * height;
   }
   open_r[0] = radius / cos(s_top[0]);
   cyln->dhmx[1] = tan(s_top[0]); // -sin(s_top[0]);
   str_ellp_setup(cyln->top, origin, s_top, 0.0, open_r);
   cyln->rotate_f[2] = (s_top[0] < TOLERANCE)? 0 : 1;
   cyln->rotate_f[3] = (ABS(s_base[0] - s_top[0]) < TOLERANCE)? 0 : 1;
 }

/* Print summary or detailed information of a cylinder(s): *********************

 str_cyln_fprintf, str_cylns_fprintf

 Print detailed information of a cylinder or a summary of multiple cylinders.

 INPUT:
 cyln - Pointer to the cylinder structure;

 OUTPUT:
 None. Prints the cylinder parameters.

*******************************************************************************/

 void
 str_cyln_fprintf 
 (
   FILE * odv,
   struct str_cyln const * cyln,
   int const indent
 )
 {
   char pre_0[STRMXLEN] = "";
   char pre_1[STRMXLEN] = "";
   for (size_t i = 0; i < indent; i++)
     strcat(pre_0, "  ");
   strncpy(pre_1, pre_0, STRMXLEN);
   strcat(pre_1, "  ");

   fprintf(odv, "%sCylinder parameters:\n", pre_0);
   fprintf(odv, "%sOrigin:   % .2e, % .2e, % .2e (m)\n", pre_1, 
     cyln->o[0], 
     cyln->o[1], 
     cyln->o[2]);
   fprintf(odv, "%sRadius:  % .2e (m)\n", pre_1,
     cyln->r);
   fprintf(odv, "%sHeight:  % .2e (m)\n", pre_1,
     cyln->h);
   fprintf(odv, "%sAxis:      %6.2lfº, %6.2lfº (% .2e, % .2e, % .2e)\n", pre_1, 
     cyln->s[0] * DEG, 
     cyln->s[1] * DEG, 
     cyln->u[0], 
     cyln->u[1],
     cyln->u[2]);
   fprintf(odv, "%sBase axis: %6.2lfº, %6.2lfº (% .2e, % .2e, % .2e)\n", pre_1, 
     cyln->base->s[0] * DEG,
     cyln->base->s[1] * DEG, 
     cyln->base->u[0],
     cyln->base->u[1],
     cyln->base->u[2]);
   fprintf(odv, "%sTop axis:  %6.2lfº, %6.2lfº (% .2e, % .2e, % .2e)\n", pre_1, 
     cyln->top->s[0] * DEG,
     cyln->top->s[1] * DEG, 
     cyln->top->u[0],
     cyln->top->u[1],
     cyln->top->u[2]);
   fprintf(odv, "%sAlpha:     %6.2lfº\n", pre_1, cyln->alp * DEG);
   fprintf(odv, "%sRotate:    %d, %d, %d, %d\n", pre_1,
     cyln->rotate_f[0],
     cyln->rotate_f[1], 
     cyln->rotate_f[2],
     cyln->rotate_f[3]);
   fprintf(odv, "%sClosed:    %d, %d\n", pre_1, cyln->closed_f[0], 
     cyln->closed_f[1]);
   mat_fprintf(odv, (double const **) cyln->M, 3, 3, 1.0, "Rotation matrix", 
     indent + 1);
 }

 void
 str_cylns_fprintf
 (
   FILE * odv,
   struct str_cyln const ** cylns,
   size_t n,
   int const indent
 )
 {
   char pre_0[STRMXLEN] = "";
   for (size_t i = 0; i < indent; i++)
     strcat(pre_0, "  ");

   fprintf(odv, "%sNumber of cylinders: %d\n", pre_0, n);
   fprintf(odv, "%s                type                     origin (m)       "
     "axis (º)  alpha (º)  Radius (m)  Height (m)  Closed      Rotate\n", 
     pre_0);
   for(size_t i = 0; i < n; i++)
   {
     fprintf(odv, "%s%02d  %s  % .2e,% .2e,% .2e  %6.2lf,%6.2lf     %6.2lf   "
       "% .2e   % .2e    %d, %d  %d, %d, %d, %d\n",
       pre_0,
       i + 1,
       (cylns[i]->rotate_f[1] || cylns[i]->rotate_f[2]) ? 
         "Oblique-Circular" : "  Right-Circular",
       cylns[i]->o[0], 
       cylns[i]->o[1], 
       cylns[i]->o[2], 
       cylns[i]->s[0] * DEG, 
       cylns[i]->s[1] * DEG, 
       cylns[i]->alp  * DEG,
       cylns[i]->r,
       cylns[i]->h,
       cylns[i]->closed_f[0],
       cylns[i]->closed_f[1],
       cylns[i]->rotate_f[0],
       cylns[i]->rotate_f[1],
       cylns[i]->rotate_f[2],
       cylns[i]->rotate_f[3]);
   }
 }

