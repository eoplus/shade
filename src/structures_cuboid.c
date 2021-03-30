
/*******************************************************************************
 structures_cuboid.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 Version: 1.6
 Date: 2021-03-25
 License: GPL-3.0

 Functions to read, write and setup the internal representation of a cuboid. 
 The definition of the cone struct is detailed in the structures_cuboid.h header
 file.

 The functions are:
 str_cubd_alloc    - Allocate structure;
 str_cubd_free     - Release memory of structure;
 str_cubd_read     - Read parameters from input file;
 str_cubd_setup    - Initialize structure with parameter values;
 str_cubd_fprintf  - Formatted detailed print of a structure;
 str_cubds_fprintf - Formatted summary print of an array of structures;

*******************************************************************************/

 #include <stdio.h>		// FILE
 #include <stdlib.h>		// exit, malloc, free
 #include <string.h>		// strcat, strncpy
 #include <math.h>		// cos, sin

 #include "config.h"		// STRMXLEN, TOLERANCE
 #include "aux.h"		// ABS, SIGN
 #include "constants.h"		// RAD, DEG, K_2PI
 #include "memory.h"		// calloc_2d, free_2d
 #include "geometry.h"		// sph_to_cos_unit, mat_transpose, mat_fprintf
 #include "rotation.h"		// rot_mat_ZYZ
 #include "structures.h"	// str_ellp
 #include "structures_cuboid.h"	// str_cubd

/* Allocate cuboid struct: *****************************************************

 str_cubd_alloc

 Allocate memory for a cuboid structure (and its components).

 INPUT:
 None.

 OUTPUT:
 Pointer to a cuboid structure.

*******************************************************************************/

 struct str_cubd* 
 str_cubd_alloc
 ( void )
 {
   struct str_cubd * cubd = 
     (struct str_cubd *) malloc( sizeof(struct str_cubd) );
   cubd->base = str_rect_alloc();
   cubd->top = str_rect_alloc();
   cubd->M = calloc_2d(3, 3, "cubd->M in str_cubd_alloc (structures_cubd.h)");
   cubd->alp = 0.0;
   cubd->closed_f[0] = 1;
   cubd->closed_f[1] = 1;
   cubd->rotate_f[0] = 0;
   cubd->rotate_f[1] = 0;
   cubd->rotate_f[2] = 0;
   cubd->rotate_f[3] = 0;

   return cubd;
 }

/* Deallocate cuboid struct: ***************************************************

 str_cubd_free

 Release memory of an cuboid structure.

 INPUT:
 cubd - Pointer to an cuboid structure pointer.

 OUTPUT:
 None. Releases memory block and set the cuboid pointer to NULL.

*******************************************************************************/

 void 
 str_cubd_free
 (
   struct str_cubd ** cubd
 )
 {
   str_rect_free( &(*cubd)->base );
   str_rect_free( &(*cubd)->top );
   free_2d( 3, &(*cubd)->M );
   free( *cubd );
   *cubd = NULL;
 }

/* Read cuboid parameters from input file: *************************************

 str_cubd_read

 Reads an cuboid from a formatted input file.

 INPUT:
 fi      - File stream pointer;
 fpos    - File stream position pointer;
 origin  - Pointer to array to receive center of the cuboid;
 axis    - Pointer to array to receive orientation of the cuboid' normal;
 alpha   - Pointer to azimuthal angle around the cuboid' normal;
 lengths - Pointer to array to size in each axis;
 base_s  - Pointer to array to receive the spherical coordinates of the base
           openings relative to cuboid' axis;
 top_s   - Pointer to array to receive the spherical coordinates of the top
           openings relative to cuboid' axis;
 closed  - Pointer to array to receive flag if base and top are open or
           closed;

 OUTPUT:
 None. Updates the values of origin, axis, radius and alpha.

*******************************************************************************/

 void
 str_cubd_read
 (
   FILE * fi, 
   long int * fpos,
   double * origin,
   double * axis,
   double * alpha,
   double * lengths,
   double * s_base,
   double * s_top,
   int * closed
 )
 {
   fseek(fi, *fpos, SEEK_SET);
   fscanf(fi, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d", 
     &origin[0],
     &origin[1],
     &origin[2], 
     &axis[0],
     &axis[1],
     alpha,
     &lengths[0],
     &lengths[1],
     &lengths[2],
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

   if ( (axis[0] < 0.0) || (axis[0] > K_PI) || 
        (axis[1] < 0.0) || (axis[1] > K_2PI) )
   {
     printf("\nERROR: The orientation axis of the cuboid must be between "
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

   // In the current version, open boxes and paralelepipeds are not implemented:
   if ( (s_base[0] != 0.0) ||  (s_base[1] != 0.0) )
   {
      printf("\nERROR: In the current version, only rectangular cuboids are "
        "implemented. The orientation axis of the cuboid' base must be "
       " 0º for the polar and azimuth angles.\n\n");
     exit(-1);
   }
   if ( (s_top[0] != 0.0) ||  (s_top[1] != 0.0) )
   {
      printf("\nERROR: In the current version, only rectangular cuboids are "
        "implemented. The orientation axis of the cuboid' top must be "
       " 0º for the polar and azimuth angles.\n\n");
     exit(-1);
   }
   if ( !closed[0] || !closed[1] )
   {
      printf("\nERROR: In the current version, only closed cuboids are "
        "implemented.\n\n");
     exit(-1);
   }
   /*
   if ( (s_base[0] < 0.0) || (s_base[0] > (60.0 * RAD)) || 
        ( (s_base[1] != 0.0) && (s_base[1] != K_2PI) ) )
   {
      printf("\nERROR: The orientation axis of the cuboid' base must be between"
       " 0º and 60º for the polar angle and 0º or 180º for the azimuthal "
       "angle.\n\n");
     exit(-1);
   }
   if ( (s_top[0] < 0.0) || (s_top[0] > (60.0 * RAD)) || 
        ( (s_top[1] != 0.0) && (s_top[1] != K_2PI) ) )
   {
      printf("\nERROR: The orientation axis of the cuboid' top must be between"
       " 0º and 60º for the polar angle and 0º or 180º for the azimuthal "
       "angle.\n\n");
     exit(-1);
   }
   */
   if ( (lengths[0] < TOLERANCE) || (lengths[1] < TOLERANCE) || 
        (lengths[2] < TOLERANCE) )
   {
     printf("\nERROR: The X, Y and Z lengths must be greater than 0.\n\n");
     exit(-1);
   }
 }

/* Setup cuboid: ***************************************************************

 str_cubd_setup

 Setup the cuboid parameter structure. The details of the cuboid structure
 are presented in the structures.h header file.

 INPUT:
 cubd    - Pointer to the cuboid structure;
 origin  - Pointer to array to receive center of cuboid;
 axis    - Pointer to array to receive orientation of the cuboid normal;
 alpha   - Pointer to azimuthal angle around the cuboid' normal;
 lengths - Pointer to array to receive the sizes in each axis.
 base_s  - Pointer to array to receive orientation of the cuboid' base normal 
           relative to its axis;
 top_s   - Pointer to array to receive orientation of the cuboid' top normal 
           relative to its axis;
 closed  - Pointer to array to receive flag if openings are open or close.

 OUTPUT:
 None. Updates the rectangle structure pointed by cubd.

*******************************************************************************/

 void 
 str_cubd_setup
 (
   struct str_cubd * cubd,
   double * origin,
   double const * axis,
   double const alpha,
   double const * lengths,
   double const * s_base,
   double const * s_top,
   int const    * closed
 )
 {
   for (int i = 0; i < 2; i++)
   {
     cubd->o[i] = origin[i];
     cubd->s[i] = axis[i];
     cubd->hl[i] = 0.5 * lengths[i];
     cubd->closed_f[i] = closed[i];
     cubd->s_base[i] = s_base[i];
     cubd->s_top[i] = s_top[i];
   }
   cubd->o[2] = origin[2];
   cubd->s[2] = 1.0;
   cubd->hl[2] = 0.5 * lengths[2];
   cubd->alp = alpha;
   cubd->x[0] = -lengths[0] * 0.5;
   cubd->x[1] =  lengths[0] * 0.5;
   cubd->y[0] = -lengths[1] * 0.5;
   cubd->y[1] =  lengths[1] * 0.5;
   cubd->z[0] = -lengths[2] * 0.5;
   cubd->z[1] =  lengths[2] * 0.5;

   sph_to_cos_unit(cubd->u, cubd->s);
   if ( (ABS(cubd->u[2]) < TOLERANCE) && (axis[1] < TOLERANCE) &&
        (alpha < TOLERANCE) )
   {
     cubd->u[0] = 0.0;
     cubd->u[1] = 0.0;
     cubd->u[2] = 1.0;
     cubd->M[0][0] = SIGN(cubd->u[2]);
     cubd->M[1][1] = SIGN(cubd->u[2]);
     cubd->M[2][2] = SIGN(cubd->u[2]);
     cubd->rotate_f[0] = 0;
   }
   else
   {
     double **R = calloc_2d(3, 3, "R in str_cubd_setup (structures_cuboid.c)");
     rot_mat_ZYZ (R, cubd->u, cubd->s, alpha);
     mat_transpose(cubd->M, (double const **) R, 3, 3);
     free_2d(3, &R);
     cubd->rotate_f[0] = 1;
   }

   // Set openings:
   str_rect_setup(cubd->base, origin, s_base, 0.0, lengths);
   cubd->rotate_f[1] = (s_base[0] < TOLERANCE)? 0 : 1;
   for (int i = 0; i < 3; i++)
   {
     origin[i] += cubd->u[i] * lengths[2];
   }
   str_rect_setup(cubd->top, origin, s_top, 0.0, lengths);
   cubd->rotate_f[2] = (s_top[0] < TOLERANCE)? 0 : 1;
   cubd->rotate_f[3] = (ABS(s_base[0] - s_top[0]) < TOLERANCE)? 0 : 1;
 }

/* Print summary or detailed information of a cone(s): *************************

 str_cubd_fprintf, str_cubds_fprintf

 Print detailed information of a cuboid or a summary of multiple cuboids.

 INPUT:
 cubd - Pointer to the cuboid structure;

 OUTPUT:
 None. Prints the cuboid parameters.

*******************************************************************************/

 void
 str_cubd_fprintf 
 (
   FILE * odv,
   struct str_cubd const * cubd,
   int const indent 
 )
 {
   char pre_0[STRMXLEN] = "";
   char pre_1[STRMXLEN] = "";
   for (int i = 0; i < indent; i++)
     strcat(pre_0, "  ");
   strncpy(pre_1, pre_0, STRMXLEN);
   strcat(pre_1, "  ");

   fprintf(odv, "%sCuboid parameters:\n", pre_0);
   fprintf(odv, "%sOrigin:   % .2e, % .2e, % .2e (m)\n", pre_1, 
     cubd->o[0], 
     cubd->o[1], 
     cubd->o[2]);
   fprintf(odv, "%sX-range:  % .2e, % .2e  (%.2e m)\n", pre_1,
     cubd->o[0]-cubd->hl[0], 
     cubd->o[0]+cubd->hl[0], 
     2.0 * cubd->hl[0]);
   fprintf(odv, "%sY-range:  % .2e, % .2e  (%.2e m)\n", pre_1, 
     cubd->o[1]-cubd->hl[1],
     cubd->o[1]+cubd->hl[1], 
     2.0 * cubd->hl[1]);
   fprintf(odv, "%sZ-range:  % .2e, % .2e  (%.2e m)\n", pre_1, 
     cubd->o[2]-cubd->hl[2],
     cubd->o[2]+cubd->hl[2], 
     2.0 * cubd->hl[2]);
   fprintf(odv, "%sAxis:      %6.2lfº, %6.2lfº (% .2e, % .2e, % .2e)\n", pre_1, 
     cubd->s[0] * DEG, 
     cubd->s[1] * DEG, 
     cubd->u[0], 
     cubd->u[1],
     cubd->u[2]);
   fprintf(odv, "%sBase axis: %6.2lfº, %6.2lfº (% .2e, % .2e, % .2e)\n", pre_1, 
     cubd->base->s[0] * DEG,
     cubd->base->s[1] * DEG, 
     cubd->base->u[0],
     cubd->base->u[1],
     cubd->base->u[2]);
   fprintf(odv, "%sTop axis:  %6.2lfº, %6.2lfº (% .2e, % .2e, % .2e)\n", pre_1, 
     cubd->top->s[0] * DEG,
     cubd->top->s[1] * DEG, 
     cubd->top->u[0],
     cubd->top->u[1],
     cubd->top->u[2]);
   fprintf(odv, "%sAlpha:     %6.2lfº\n", pre_1, cubd->alp * DEG);
   fprintf(odv, "%sRotate:    %d, %d, %d, %d\n", pre_1, cubd->rotate_f[0],
     cubd->rotate_f[1], cubd->rotate_f[2], cubd->rotate_f[3]);
   fprintf(odv, "%sClosed:    %d, %d\n", pre_1, cubd->closed_f[0], 
     cubd->closed_f[1]);
   mat_fprintf(odv, (double const **) cubd->M, 3, 3, 1.0, "Rotation matrix", 
     indent + 1);
 }

 void
 str_cubds_fprintf
 (
   FILE * odv,
   struct str_cubd const ** cubds,
   int n,
   int const indent 
 )
 {
   char pre_0[STRMXLEN] = "";
   for (int i = 0; i < indent; i++)
     strcat(pre_0, "  ");

   fprintf(odv, "%sNumber of cuboids: %d\n", pre_0, n);
   fprintf(odv, "%s         type                     origin (m)       axis (º) "
     " alpha (º)                     Length (m)  Closed      Rotate\n", pre_0);
   for(int i = 0; i < n; i++)
   {
     fprintf(odv, "%s%02d  %s  % .2e,% .2e,% .2e  %6.2lf,%6.2lf     %6.2lf  "
       "% .2e,% .2e,% .2e    %d, %d  %d, %d, %d, %d\n",
       pre_0,
       i + 1,
       ( ABS(cubds[i]->hl[0] - cubds[i]->hl[1]) < TOLERANCE ) &&
       ( ABS(cubds[i]->hl[0] - cubds[i]->hl[2]) < TOLERANCE ) ? 
         "     Cube" : "   Cuboid",
       cubds[i]->o[0], 
       cubds[i]->o[1], 
       cubds[i]->o[2], 
       cubds[i]->s[0] * DEG, 
       cubds[i]->s[1] * DEG, 
       cubds[i]->alp  * DEG,
       cubds[i]->hl[0] * 2.0,
       cubds[i]->hl[1] * 2.0,
       cubds[i]->hl[2] * 2.0,
       cubds[i]->closed_f[0],
       cubds[i]->closed_f[1],
       cubds[i]->rotate_f[0],
       cubds[i]->rotate_f[1],
       cubds[i]->rotate_f[2],
       cubds[i]->rotate_f[3]);
   }
 }

