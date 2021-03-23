
/*******************************************************************************
 structures_cone.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 Functions to read, write and setup the internal representation of a positive
 finite (truncated) cone. The definition of the cone struct is detailed in the
 structures_cone.h header file.

 The functions are:
 str_cone_alloc   - Allocate structure;
 str_cone_free    - Release memory of structure;
 str_cone_read    - Read parameters from input file;
 str_cone_setup   - Initialize structure with parameter values;
 str_cone_fprintf  - Formatted detailed print of a structure;
 str_cones_fprintf - Formatted summary print of an array of structures;

*******************************************************************************/

 #include <stdio.h>		// FILE
 #include <stddef.h>		// size_t
 #include <stdlib.h>		// exit, malloc, free
 #include <string.h>		// strcat, strncpy
 #include <math.h>		// M_PI, M_PI_2, cos, sin

 #include "config.h"		// STRMXLEN, TOLERANCE
 #include "aux.h"		// ABS, SIGN
 #include "constants.h"		// RAD, DEG, K_2PI
 #include "memory.h"		// calloc_2d, free_2d
 #include "geometry.h"		// sph_to_cos_unit, mat_transpose, mat_fprintf
 #include "rotation.h"		// rot_mat_ZYZ
 #include "structures.h"	// str_ellp
 #include "structures_cone.h"	// str_cone

/* Allocate cone struct: *******************************************************

 str_cone_alloc

 Allocate memory for a cone structure (and its components).

 INPUT:
 None.

 OUTPUT:
 Pointer to a str_cone struct.

*******************************************************************************/

 struct str_cone * 
 str_cone_alloc
 ( void )
 {
   struct str_cone * cone = (struct str_cone *) malloc(sizeof(struct str_cone));
   cone->base = str_ellp_alloc();
   cone->top  = str_ellp_alloc();
   cone->M = calloc_2d(3, 3, "cone->M in str_cone_alloc (structures_cone.h)");
   cone->alp = 0.0;
   cone->closed_f[0] = 1;
   cone->closed_f[1] = 1;
   cone->rotate_f[0] = 0;
   cone->rotate_f[1] = 0;
   cone->rotate_f[2] = 0;
   cone->rotate_f[3] = 0;

   return cone;
 }

/* Deallocate cone struct: *****************************************************

 str_cone_free

 Release memory of a cone structure.

 INPUT:
 cone - Cone parameters
        Pointer to pointer to str_cone struct.

 OUTPUT:
 None. Releases memory block and set the str_cone pointer to NULL.

*******************************************************************************/

 void 
 str_cone_free
 ( struct str_cone ** cone )
 {
   str_ellp_free( &(*cone)->base );
   str_ellp_free( &(*cone)->top );
   free_2d( 3, &(*cone)->M );
   free( *cone );
   *cone = NULL;
 }

/* Read cone parameters from input file: ***************************************

 str_cone_read

 Reads a cone from a formatted standard input file.

 INPUT:
 fi      - File stream pointer;
 fpos    - File stream position pointer;
 origin  - Pointer to array to receive center of the cone's base;
 axis    - Pointer to array to receive orientation of the cone's normal;
 alpha   - Pointer to azimuthal angle around the cone's normal;
 psi
 height  - Pointer to double to receive the height;
 base_s  - Pointer to array to receive the spherical coordinates of the base
           openings relative to cone's axis;
 top_s   - Pointer to array to receive the spherical coordinates of the top
           openings relative to cone's axis;
 closed  - Pointer to array to receive flag if base and top are open or
           closed;

 OUTPUT:
 None. Updates the values of origin, axis, alpha, radius, height, base_s, top_s,
 and closed.

*******************************************************************************/

 void
 str_cone_read
 (
   FILE * fi, 
   long int * fpos,
   double * origin,
   double * axis,
   double * alpha,
   double * psi,
   double * height,
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
     psi,
     &height[0],
     &height[1],
     &s_base[0], 
     &s_base[1],
     &s_top[0],
     &s_top[1],
     &closed[0],
     &closed[1]);
   *fpos = ftell(fi);

   axis[0] *= RAD;
   axis[1] *= RAD;
   *alpha  *= RAD;
   *psi *= RAD;
   s_base[0] *= RAD;
   s_base[1] *= RAD;
   s_top[0] *= RAD;
   s_top[1] *= RAD;

   if ( (axis[0] < 0.0) || (axis[0] > M_PI) || 
        (axis[1] < 0.0) || (axis[1] > K_2PI) )
   {
     printf("\nERROR: The orientation axis of the cone must be between "
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
     printf("\nERROR: The orientation axis of the cone's base must be "
       "between 0º and 60º for the polar angle and 0º or 180º for the "
       "azimuthal angle.\n\n");
     exit(-1);
   }
   if ( (s_top[0] < 0.0) || (s_top[0] > (60.0 * RAD)) || 
        ( (s_top[1] > TOLERANCE) && (ABS(s_top[1] - M_PI) > TOLERANCE) ) )
   {
     printf("\nERROR: The orientation axis of the cone's top must be "
       "between 0º and 60º for the polar angle and 0º or 180º for the "
       "azimuthal angle.\n\n");
     exit(-1);
   }
   if ( (*psi) <= 0.0 || (*psi) >= M_PI_2 )
   {
     printf("\nERROR: The half-angle of the cone must be > 0º and < 90º.\n\n");
     exit(-1);
   }
   if ( height[0] < TOLERANCE || height[1] <= height[0] )
   {
     printf("\nERROR: The cone heights must be positive and height of the top "
       "must be greater than the height of the base.\n\n");
     exit(-1);
   }
 }

/* Setup cone: *****************************************************************

 str_cone_setup

 Setup the cone parameter structure. The details of the cone structure
 are presented in the structures.h header file.

 INPUT:
 cone    - Pointer to the cone structure;
 origin  - Pointer to array containing the center of the cone's base;
 axis    - Pointer to array containing the orientation of the cone's normal;
 alpha   - Pointer to azimuthal angle around the cone's normal;
 radius  - Pointer to array containing the radius of the cone;
 height  - Pointer to array containing the height of the cone along its
           axis;
 base_s  - Pointer to array containing the orientation of the cone's base
           normal relative to its axis;
 top_s   - Pointer to array containing the orientation of the cone's top
           normal relative to its axis;
 closed  - Pointer to array containing the flag if openings are open or close.

 OUTPUT:
 None. Updates the cone structure pointed by cone.

*******************************************************************************/

 void 
 str_cone_setup
 (
   struct str_cone * cone,
   double * origin,
   double const * axis,
   double const alpha,
   double const psi,
   double const * height,
   double const * s_base,
   double const * s_top,
   int const * closed
 )
 {
   for (size_t i = 0; i < 2; i++)
   {
     cone->o[i] = origin[i];
     cone->s[i] = axis[i];
     cone->closed_f[i] = closed[i];
     cone->s_base[i] = s_base[i];
     cone->s_top[i] = s_top[i];
     cone->h[i] = height[i];
   }
   cone->o[2] = origin[2];
   cone->s[2] = 1.0;
   cone->psi = psi;
   cone->mu = cos(psi);
   cone->musq = cone->mu * cone->mu;
   cone->alp = alpha;

   sph_to_cos_unit(cone->u, cone->s);
   if ( (1.0 - ABS(cone->u[2]) < TOLERANCE) && 
        (axis[1] < TOLERANCE) && (alpha < TOLERANCE) )
   {
     cone->M[0][0] = SIGN(cone->u[2]);
     cone->M[1][1] = SIGN(cone->u[2]);
     cone->M[2][2] = SIGN(cone->u[2]);
     cone->rotate_f[0] = 0;
   }
   else
   {
     double ** R = calloc_2d(3, 3, "R in str_cone_setup (structures_cone.c)");
     rot_mat_ZYZ (R, cone->u, cone->s, alpha);
     mat_transpose(cone->M, (double const **) R, 3, 3);
     free_2d(3, &R);
     cone->rotate_f[0] = 1;
   }

   // Set base opening:
   double open_r[2];
   for (size_t i = 0; i < 3; i++)
   {
     origin[i] += cone->u[i] * height[0];
   }
   open_r[0] = sin(psi) * height[0];
   open_r[1] = open_r[0];
   open_r[0] /= cos(s_base[0]);
   cone->dhmx[0] = -sin(s_base[0]);
   str_ellp_setup(cone->base, origin, s_base, 0.0, open_r);
   cone->rotate_f[1] = (s_base[0] < TOLERANCE)? 0 : 1;

   // Set top opening:
   for (size_t i = 0; i < 3; i++)
   {
     origin[i] += cone->u[i] * (height[1] - height[0]);
   }
   open_r[0] = sin(psi) * height[1];
   open_r[1] = open_r[0];
   open_r[0] /= cos(s_top[0]);
   cone->dhmx[1] = -sin(s_top[0]);
   str_ellp_setup(cone->top, origin, s_top, 0.0, open_r);
   cone->rotate_f[2] = (s_top[0] < TOLERANCE)? 0 : 1;
   cone->rotate_f[3] = (ABS(s_base[0] - s_top[0]) < TOLERANCE)? 0 : 1;
 }

/* Print summary or detailed information of a cone(s): *************************

 str_cone_fprintf, str_cones_printf

 Print detailed information of a cone or a summary of multiple cones.

 INPUT:
 cone - Pointer to the cone structure;

 OUTPUT:
 None. Prints the cone parameters.

*******************************************************************************/

 void
 str_cone_fprintf 
 (
   FILE * odv,
   struct str_cone const * cone,
   int const indent
 )
 {
   char pre_0[STRMXLEN] = "";
   char pre_1[STRMXLEN] = "";
   for (size_t i = 0; i < indent; i++)
     strcat(pre_0, "  ");
   strncpy(pre_1, pre_0, STRMXLEN);
   strcat(pre_1, "  ");

   fprintf(odv, "%sCone parameters:\n", pre_0);
   fprintf(odv, "%sOrigin:   % .2e, % .2e, % .2e (m)\n", pre_1, 
     cone->o[0], 
     cone->o[1], 
     cone->o[2]);
   fprintf(odv, "%sPsi:      %6.2lfº\n", pre_1, cone->psi * DEG);
   fprintf(odv, "%sHeight:  % .2e, % .2e (m)\n", pre_1,
     cone->h[0],
     cone->h[1]);
   fprintf(odv, "%sAxis:      %6.2lfº, %6.2lfº (% .2e, % .2e, % .2e)\n", pre_1, 
     cone->s[0] * DEG, 
     cone->s[1] * DEG, 
     cone->u[0], 
     cone->u[1],
     cone->u[2]);
   fprintf(odv, "%sBase axis: %6.2lfº, %6.2lfº (% .2e, % .2e, % .2e)\n", pre_1, 
     cone->base->s[0] * DEG,
     cone->base->s[1] * DEG, 
     cone->base->u[0],
     cone->base->u[1],
     cone->base->u[2]);
   fprintf(odv, "%sTop axis:  %6.2lfº, %6.2lfº (% .2e, % .2e, % .2e)\n", pre_1, 
     cone->top->s[0] * DEG,
     cone->top->s[1] * DEG, 
     cone->top->u[0],
     cone->top->u[1],
     cone->top->u[2]);
   fprintf(odv, "%sAlpha:     %6.2lfº\n", pre_1, cone->alp * DEG);
   fprintf(odv, "%sRotate:    %d, %d, %d, %d\n", pre_1,
     cone->rotate_f[0],
     cone->rotate_f[1], 
     cone->rotate_f[2],
     cone->rotate_f[3]);
   fprintf(odv, "%sClosed:    %d, %d\n", pre_1, 
     cone->closed_f[0],
     cone->closed_f[1]);
   mat_fprintf(odv, (double const **) cone->M, 3, 3, 1.0, "Rotation matrix",
     indent+1);
 }

 void
 str_cones_fprintf
 (
   FILE * odv,
   struct str_cone const ** cones,
   size_t n,
   int const indent
 )
 {
   char pre_0[STRMXLEN] = "";
   for (size_t i = 0; i < indent; i++)
     strcat(pre_0, "  ");

   fprintf(odv, "%sNumber of cones: %d\n", pre_0, n);
   fprintf(odv, "%s                type                     origin (m)       "
     "axis (º)  alpha (º)  Psi (º)           Height (m)  Closed      Rotate\n", 
     pre_0);

   for (size_t i = 0; i < n; i++)
   {
     fprintf(odv, "%s%02d  %s  % .2e,% .2e,% .2e  %6.2lf,%6.2lf     %6.2lf   "
       "%6.2lf  % .2e,% .2e    %d, %d  %d, %d, %d, %d\n",
       pre_0,
       i + 1,
       (cones[i]->rotate_f[1] || cones[i]->rotate_f[2]) ? 
         "Oblique-Circular" : "  Right-Circular",
       cones[i]->o[0], 
       cones[i]->o[1], 
       cones[i]->o[2], 
       cones[i]->s[0] * DEG, 
       cones[i]->s[1] * DEG, 
       cones[i]->alp * DEG,
       cones[i]->psi * DEG,
       cones[i]->h[0],
       cones[i]->h[1],
       cones[i]->closed_f[0],
       cones[i]->closed_f[1],
       cones[i]->rotate_f[0],
       cones[i]->rotate_f[1],
       cones[i]->rotate_f[2],
       cones[i]->rotate_f[3]);
   }
 }

