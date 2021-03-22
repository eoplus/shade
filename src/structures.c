
/*******************************************************************************
 structures.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 General functions to read, write and setup the internal representation of bi- 
 and three-dimensional geometrical structures. The datatype definitions of those
 structures is in the structures.h header file.

 The functions are:
 str_x_alloc   - Allocate structure;
 str_x_free    - Release memory of structure;
 str_x_read    - Read parameters from input file;
 str_x_setup   - Initialize structure with parameter values;
 str_x_printf  - Formatted detailed print of a structure;
 str_xs_printf - Formatted summary print of an array of structures;

 where "x" is:
 rect - rectangles;
 ellp - ellipses;
 trig - triangles; // Not yet implemented;
 cubd - cuboids;
 cyln - cylinders;
 sphd - spheroids; // Not yet implemented;
 prsm - prisms;    // Not yet implemented;

*******************************************************************************/

 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>

 #include "config.h"
 #include "constants.h"
 #include "aux.h"
 #include "memory.h"
 #include "geometry.h"
 #include "rotation.h"
 #include "structures.h"

/* str_rect_alloc

 Allocate memory for a rectangle structure (and its components).

 INPUT:
 None.

 OUTPUT:
 Pointer to a rectangle structure.

 */

 str_rect* 
 str_rect_alloc
 (
  void
 )
 {
   str_rect* rect = (str_rect*) malloc( sizeof(str_rect) );
   rect->M = calloc_2d(3, 3, "rect->M in str_rect_alloc (structures.h)");
   rect->alp = 0.0;
   rect->rotate_f = 0;

   return rect;
 }

/* str_rect_free

 Release memory of a rectangle structure.

 INPUT:
 rect - Pointer to a rectangle structure pointer.

 OUTPUT:
 None. Releases memory block and set the rectangle pointer to NULL.

 */

 void 
 str_rect_free
 (
  str_rect **rect
 )
 {
   free_2d( 3, &(*rect)->M );
   free( (*rect) );
   (*rect) = NULL;
 }

/* str_rect_read

 Reads a rectangle from a formatted input file.

 INPUT:
 fi      - File stream pointer;
 fpos    - File stream position pointer;
 origin  - Pointer to array to receive center of rectangle;
 axis    - Pointer to array to receive orientation of the rectangle normal;
 alpha   - Pointer to azimuthal angle around rectangle's normal;
 lengths - Pointer to array to size in each axis.

 OUTPUT:
 None. Updates the values of origin, axis, xrange, yrange and alpha.

 */

 void 
 str_rect_read
 (
   FILE     *fi, 
   long int *fpos,
   double   *origin,
   double   *axis,
   double   *alpha,
   double   *lengths
 )
 {
   fseek(fi, *fpos, SEEK_SET);
   fscanf(fi, "%lf %lf %lf %lf %lf %lf %lf %lf", 
     &origin[0], &origin[1], &origin[2], 
     &axis[0], &axis[1], alpha, 
     &lengths[0], &lengths[1]);
   *fpos = ftell(fi);

   axis[0] *= RAD;
   axis[1] *= RAD;
   *alpha  *= RAD;
   if ( (axis[0] < 0.0) || (axis[0] > M_PI) || 
        (axis[1] < 0.0) || (axis[1] > K_2PI) )
   {
     printf("\nERROR: The orientation axis of the rectangle must be between "
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
   if ( (lengths[0] < TOLERANCE) || (lengths[1] < TOLERANCE) )
   {
     printf("\nERROR: The X and Y lengths must be greater than 0.\n\n");
     exit(-1);
   }
 }

/* str_rect_setup

 Setup the rectangle parameter structure. The details of the rectangle structure
 are presented in the structures.h header file.

 INPUT:
 rect    - Pointer to the ellipse structure;
 origin  - Pointer to array to receive center of rectangle;
 axis    - Pointer to array to receive orientation of the rectangle normal;
 alpha   - Pointer to azimuthal angle around the rectangle's normal;
 lengths - Pointer to array to receive the sizes in each axis.

 OUTPUT:
 None. Updates the rectangle structure pointed by rect.

 */

 void 
 str_rect_setup
 (
   str_rect     *rect,
   double const *origin,
   double const *axis,
   double const alpha,
   double const *lengths
 )
 {
   for (size_t i = 0; i < 2; i++)
   {
     rect->o[i]  = origin[i];
     rect->s[i]  = axis[i];
     rect->hl[i] = 0.5 * lengths[i];
   }
   rect->o[2] = origin[2];
   rect->s[2] = 1.0;
   rect->alp  = alpha;
   rect->x[0] = -lengths[0] * 0.5;
   rect->x[1] =  lengths[0] * 0.5;
   rect->y[0] = -lengths[1] * 0.5;
   rect->y[1] =  lengths[1] * 0.5;

   if ( (axis[0] < TOLERANCE) && (axis[1] < TOLERANCE) && (alpha < TOLERANCE) )
   {
     rect->u[0] = 0.0;
     rect->u[1] = 0.0;
     rect->u[2] = 1.0;
     rect->M[0][0] = 1.0;
     rect->M[1][1] = 1.0;
     rect->M[2][2] = 1.0;
     rect->rotate_f = 0;
   }
   else
   {
     sph_to_cos_unit(rect->u, rect->s);
     double **R = calloc_2d(3, 3, "R in str_rect_setup (structures.c)");
     rot_mat_ZYZ (R, rect->u, rect->s, alpha);
     mat_transpose(rect->M, (double const **) R, 3, 3);
     free_2d(3, &R);
     rect->rotate_f = 1;
   }
 }

/* str_rect_printf, str_rects_printf

 Print detailed information of a rectangle or a summary of multiple rectangles.

 INPUT:
 rect - Pointer to the rectangle structure;

 OUTPUT:
 None. Prints the rectangle parameters.

 */

 void
 str_rect_fprintf 
 (
   FILE *odv,
   str_rect const *rect,
   int const indent 
 )
 {
   char pre_0[STRMXLEN] = "";
   char pre_1[STRMXLEN] = "";
   for (size_t i = 0; i < indent; i++)
     strcat(pre_0, "  ");
   strncpy(pre_1, pre_0, STRMXLEN);
   strcat(pre_1, "  ");

   fprintf(odv, "%sRectangle parameters:\n", pre_0);
   fprintf(odv, "%sOrigin:   % .2e, % .2e, % .2e (m)\n", pre_1, 
     rect->o[0], 
     rect->o[1], 
     rect->o[2]);
   fprintf(odv, "%sX-range:  % .2e, % .2e  (%.2e m)\n", pre_1,
     rect->x[0], 
     rect->x[1], 
     2.0 * rect->hl[0]);
   fprintf(odv, "%sY-range:  % .2e, % .2e  (%.2e m)\n", pre_1, 
     rect->y[0],
     rect->y[1], 
     2.0 * rect->hl[1]);
   fprintf(odv, "%sAxis:      %6.2lfº, %6.2lfº (% .2e, % .2e, % .2e)\n", pre_1, 
     rect->s[0] * DEG, 
     rect->s[1] * DEG, 
     rect->u[0], 
     rect->u[1],
     rect->u[2]);
   fprintf(odv, "%sAlpha:     %6.2lfº\n", pre_1, rect->alp * DEG);
   fprintf(odv, "%sRotate:    %d\n", pre_1, rect->rotate_f);
   mat_fprintf(odv, (double const **) rect->M, 3, 3, 1.0, "Rotation matrix", indent+1);
 }

 void
 str_rects_printf
 (
  str_rect const **rects,
  size_t n
 )
 {
   printf("  Number of rectangles: %d\n", n);
   printf("           type                     origin (m)       axis (º) "
     " alpha (º)           Length (m)\n");
   for(size_t i = 0; i < n; i++)
   {
     printf("  %02d  %s  % .2e,% .2e,% .2e  %6.2lf,%6.2lf     %6.2lf  "
       "% .2e,% .2e\n",
       i+1,
       ( ABS(rects[i]->hl[0] - rects[i]->hl[1]) < TOLERANCE ) ? 
         "   Square" : "Rectangle",
       rects[i]->o[0], 
       rects[i]->o[1], 
       rects[i]->o[2], 
       rects[i]->s[0] * DEG, 
       rects[i]->s[1] * DEG, 
       rects[i]->alp * DEG,
       2.0 * rects[i]->hl[0],
       2.0 * rects[i]->hl[1]);
   }
 }

/* str_ellp_alloc

 Allocate memory for an ellipse structure (and its components).

 INPUT:
 None.

 OUTPUT:
 Pointer to an ellipse structure.

 */

 str_ellp* 
 str_ellp_alloc
 (
  void
 )
 {
   str_ellp* ellp = (str_ellp*) malloc( sizeof(str_ellp) );
   ellp->M = calloc_2d(3, 3, "ellp->M in str_ellp_alloc (structures.h)");
   ellp->alp = 0.0;
   ellp->rotate_f = 0;

   return ellp;
 }

/* str_ellp_free

 Release memory of an ellipse structure.

 INPUT:
 ellp - Pointer to an ellipse structure pointer.

 OUTPUT:
 None. Releases memory block and set the ellipse pointer to NULL.

 */

 void 
 str_ellp_free
 (
  str_ellp **ellp
 )
 {
   free_2d(3, &(*ellp)->M);
   free(*ellp);
   *ellp = NULL;
 }

/* str_ellp_read

 Reads an ellipse from a formatted input file.

 INPUT:
 fi     - File stream pointer;
 fpos   - File stream position pointer;
 origin - Pointer to array to receive center of the ellipse;
 axis   - Pointer to array to receive orientation of the ellipse's normal;
 alpha  - Pointer to azimuthal angle around the ellipse's normal;
 radius - Pointer to array to receive the X and Y radius of the ellipse.

 OUTPUT:
 None. Updates the values of origin, axis, radius and alpha.

 */

 void
 str_ellp_read
 (
   FILE     *fi, 
   long int *fpos,
   double   *origin,
   double   *axis,
   double   *alpha,
   double   *radius
 )
 {
   fseek(fi, *fpos, SEEK_SET);
   fscanf(fi, "%lf %lf %lf %lf %lf %lf %lf %lf", 
     &origin[0], &origin[1], &origin[2], 
     &axis[0], &axis[1], alpha,
     &radius[0], &radius[1]);
   *fpos = ftell(fi);

   axis[0] *= RAD;
   axis[1] *= RAD;
   *alpha  *= RAD;

   if ( (axis[0] < 0.0) || (axis[0] > M_PI) || 
        (axis[1] < 0.0) || (axis[1] > K_2PI) )
   {
     printf("\nERROR: The orientation axis of the ellipse must be between "
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
   if ( radius[0] <= 0.0 || radius[1] <= 0.0 )
   {
     printf("\nERROR: Ellipse's X and Y axis radius must be greater than "
       "zero.\n\n");
     exit(-1);
   }
 }

/* str_ellp_setup

 Setup the ellipse parameter structure. The details of the ellipse structure
 are presented in the structures.h header file.

 INPUT:
 ellp   - Pointer to the ellipse structure;
 origin - Pointer to the origin (center) of the ellipse;
 axis   - Pointer to the orientation axis of the normal;
 alpha  - Pointer to azimuthal angle around ellipse's normal;
 radius - Pointer to radius of the ellipse for X and Y axis, (on its own 
          reference frame);

 OUTPUT:
 None. Updates the ellipse structure pointed by ellp.

 */

 void
 str_ellp_setup
 (
   str_ellp *ellp,
   double const *origin,
   double const *axis,
   double const alpha,
   double const *radius
 )
 {
   for (size_t i = 0; i < 2; i++)
   {
     ellp->o[i] = origin[i];
     ellp->s[i] = axis[i];
     ellp->r[i] = radius[i];
     ellp->rsq_inv[i] = 1.0 / (radius[i] * radius[i]);
   }
   ellp->o[2] = origin[2];
   ellp->s[2] = 1.0;
   ellp->alp  = alpha;

   if ( (axis[0] < TOLERANCE) && (axis[1] < TOLERANCE) && (alpha < TOLERANCE) )
   {
     ellp->u[0] = 0.0;
     ellp->u[1] = 0.0;
     ellp->u[2] = 1.0;
     ellp->M[0][0] = 1.0;
     ellp->M[1][1] = 1.0;
     ellp->M[2][2] = 1.0;
     ellp->rotate_f = 0;
   }
   else
   {
     sph_to_cos_unit(ellp->u, ellp->s);
     double **R = calloc_2d(3, 3, "R in str_ellp_setup (structures.c)");
     rot_mat_ZYZ (R, ellp->u, ellp->s, alpha);
     mat_transpose(ellp->M, (double const **) R, 3, 3);
     free_2d(3, &R);
     ellp->rotate_f = 1;
   }
 }

/* str_ellp_printf, str_ellp_printf

 Print detailed information of an ellipse or a summary of multiple ellipses.

 INPUT:
 ellp - Pointer to the ellipse structure;

 OUTPUT:
 None. Prints the ellipse parameters.

 */

 void
 str_ellp_fprintf 
 (
   FILE *odv,
   str_ellp const *ellp,
   int const indent 
 )
 {
   char pre_0[STRMXLEN] = "";
   char pre_1[STRMXLEN] = "";
   for (size_t i = 0; i < indent; i++)
     strcat(pre_0, "  ");
   strncpy(pre_1, pre_0, STRMXLEN);
   strcat(pre_1, "  ");

   fprintf(odv, "%sEllipse parameters:\n", pre_0);
   fprintf(odv, "%sOrigin:   % .2e, % .2e, % .2e (m)\n", pre_1, 
     ellp->o[0], 
     ellp->o[1], 
     ellp->o[2]);
   fprintf(odv, "%sRadius:  % .2e, % .2e  (m)\n", pre_1,
     ellp->r[0], 
     ellp->r[1]);
   fprintf(odv, "%sAxis:      %6.2lfº, %6.2lfº (% .2e, % .2e, % .2e)\n", pre_1, 
     ellp->s[0] * DEG, 
     ellp->s[1] * DEG, 
     ellp->u[0], 
     ellp->u[1],
     ellp->u[2]);
   fprintf(odv, "%sAlpha:     %6.2lfº\n", pre_1, ellp->alp * DEG);
   fprintf(odv, "%sRotate:    %d\n", pre_1, ellp->rotate_f);
   mat_fprintf(odv, (double const **) ellp->M, 3, 3, 1.0, "Rotation matrix", indent+1);
 }

 void
 str_ellps_printf
 (
  str_ellp const **ellps,
  size_t n
 )
 {
   printf("  Number of ellipses: %d\n", n);
   printf("           type                     origin (m)       axis (º) "
     " alpha (º)           Radius (m)\n");
   for(size_t i = 0; i < n; i++)
   {
     printf("  %02d  %s  % .2e,% .2e,% .2e  %6.2lf,%6.2lf     %6.2lf  "
       "% .2e,% .2e\n",
       i+1,
       ( ABS(ellps[i]->r[0] -  ellps[i]->r[1]) < TOLERANCE ) ? 
         "   Circle" : "  Ellipse",
       ellps[i]->o[0], 
       ellps[i]->o[1], 
       ellps[i]->o[2], 
       ellps[i]->s[0] * DEG, 
       ellps[i]->s[1] * DEG, 
       ellps[i]->alp  * DEG,
       ellps[i]->r[0],
       ellps[i]->r[1]);
   }
 }

/* str_cubd_alloc

 Allocate memory for a cuboid structure (and its components).

 INPUT:
 None.

 OUTPUT:
 Pointer to a cuboid structure.

 */

 str_cubd* 
 str_cubd_alloc
 (
  void
 )
 {
   str_cubd* cubd = (str_cubd*) malloc( sizeof(str_cubd) );
   cubd->base = str_rect_alloc();
   cubd->top  = str_rect_alloc();
   cubd->M   = calloc_2d(3, 3, "cubd->M in str_cubd_alloc (structures.h)");
   cubd->alp = 0.0;
   cubd->closed_f[0] = 1;
   cubd->closed_f[1] = 1;
   cubd->rotate_f[0] = 0;
   cubd->rotate_f[1] = 0;
   cubd->rotate_f[2] = 0;
   cubd->rotate_f[3] = 0;

   return cubd;
 }

/* str_cubd_free

 Release memory of an cuboid structure.

 INPUT:
 cubd - Pointer to an cuboid structure pointer.

 OUTPUT:
 None. Releases memory block and set the cuboid pointer to NULL.

 */

 void 
 str_cubd_free
 (
  str_cubd **cubd
 )
 {
   str_rect_free(&(*cubd)->base);
   str_rect_free(&(*cubd)->top);
   free_2d(3, &(*cubd)->M);
   free(*cubd);
   *cubd = NULL;
 }

/* str_cubd_read

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

 */

 void
 str_cubd_read
 (
   FILE     *fi, 
   long int *fpos,
   double   *origin,
   double   *axis,
   double   *alpha,
   double   *lengths,
   double   *s_base,
   double   *s_top,
   int      *closed
 )
 {
   fseek(fi, *fpos, SEEK_SET);
   fscanf(fi, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d", 
     &origin[0], &origin[1], &origin[2], 
     &axis[0], &axis[1], alpha,
     &lengths[0], &lengths[1], &lengths[2],
     &s_base[0], &s_base[1],
     &s_top[0], &s_top[1],
     &closed[0], &closed[1]);
   *fpos = ftell(fi);

   axis[0] *= RAD;
   axis[1] *= RAD;
   *alpha  *= RAD;
   s_base[0] *= RAD;
   s_base[1] *= RAD;
   s_top[0] *= RAD;
   s_top[1] *= RAD;

   if ( (axis[0] < 0.0) || (axis[0] > M_PI) || 
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
   if ( (lengths[0] < TOLERANCE) || 
        (lengths[1] < TOLERANCE) || 
        (lengths[2] < TOLERANCE) )
   {
     printf("\nERROR: The X, Y and Z lengths must be greater than 0.\n\n");
     exit(-1);
   }
 }

/* str_cubd_setup

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

 */

 void 
 str_cubd_setup
 (
   str_cubd     *cubd,
   double       *origin,
   double const *axis,
   double const alpha,
   double const *lengths,
   double const *s_base,
   double const *s_top,
   int const    *closed
 )
 {
   for (size_t i = 0; i < 2; i++)
   {
     cubd->o[i]  = origin[i];
     cubd->s[i]  = axis[i];
     cubd->hl[i] = 0.5 * lengths[i];
     cubd->closed_f[i] = closed[i];
     cubd->s_base[i] = s_base[i];
     cubd->s_top[i] = s_top[i];
   }
   cubd->o[2] = origin[2];
   cubd->s[2] = 1.0;
   cubd->hl[2] = 0.5 * lengths[2];
   cubd->alp  = alpha;
   cubd->x[0] = -lengths[0] * 0.5;
   cubd->x[1] =  lengths[0] * 0.5;
   cubd->y[0] = -lengths[1] * 0.5;
   cubd->y[1] =  lengths[1] * 0.5;
   cubd->z[0] = -lengths[2] * 0.5;
   cubd->z[1] =  lengths[2] * 0.5;

   if ( (axis[0] < TOLERANCE) && (axis[1] < TOLERANCE) && (alpha < TOLERANCE) )
   {
     cubd->u[0] = 0.0;
     cubd->u[1] = 0.0;
     cubd->u[2] = 1.0;
     cubd->M[0][0] = 1.0;
     cubd->M[1][1] = 1.0;
     cubd->M[2][2] = 1.0;
     cubd->rotate_f[0] = 0;
   }
   else
   {
     sph_to_cos_unit(cubd->u, cubd->s);
     double **R = calloc_2d(3, 3, "R in str_cubd_setup (structures.c)");
     rot_mat_ZYZ (R, cubd->u, cubd->s, alpha);
     mat_transpose(cubd->M, (double const **) R, 3, 3);
     free_2d(3, &R);
     cubd->rotate_f[0] = 1;
   }

   // Set openings:
   str_rect_setup(cubd->base, origin, s_base, 0.0, lengths);
   cubd->rotate_f[1] = (s_base[0] < TOLERANCE)? 0 : 1;
   for (size_t i = 0; i < 3; i++)
   {
     origin[i] += cubd->u[i] * lengths[2];
   }
   str_rect_setup(cubd->top, origin, s_top, 0.0, lengths);
   cubd->rotate_f[2] = (s_top[0] < TOLERANCE)? 0 : 1;
   cubd->rotate_f[3] = (ABS(s_base[0] - s_top[0]) < TOLERANCE)? 0 : 1;
 }

/* str_cubd_printf, str_cubd_printf

 Print detailed information of a cuboid or a summary of multiple cuboids.

 INPUT:
 cubd - Pointer to the cuboid structure;

 OUTPUT:
 None. Prints the cuboid parameters.

 */

 void
 str_cubd_fprintf 
 (
   FILE *odv,
   str_cubd const *cubd,
   int const indent 
 )
 {
   char pre_0[STRMXLEN] = "";
   char pre_1[STRMXLEN] = "";
   for (size_t i = 0; i < indent; i++)
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
   fprintf(odv, "%sClosed:    %d, %d\n", pre_1, cubd->closed_f[0], cubd->closed_f[1]);
   mat_fprintf(odv, (double const **) cubd->M, 3, 3, 1.0, "Rotation matrix", indent+1);
 }

 void
 str_cubds_printf
 (
  str_cubd const **cubds,
  size_t n
 )
 {
   printf("  Number of cuboids: %d\n", n);
   printf("           type                     origin (m)       axis (º) "
     " alpha (º)                     Length (m)  Closed      Rotate\n");
   for(size_t i = 0; i < n; i++)
   {
     printf("  %02d  %s  % .2e,% .2e,% .2e  %6.2lf,%6.2lf     %6.2lf  "
       "% .2e,% .2e,% .2e    %d, %d  %d, %d, %d, %d\n",
       i+1,
       ( ABS(cubds[i]->hl[0] - cubds[i]->hl[1]) < TOLERANCE ) &&
       ( ABS(cubds[i]->hl[0] - cubds[i]->hl[2]) < TOLERANCE ) ? 
         "     Cube" : "   Cuboid",
       cubds[i]->o[0], 
       cubds[i]->o[1], 
       cubds[i]->o[2], 
       cubds[i]->s[0] * DEG, 
       cubds[i]->s[1] * DEG, 
       cubds[i]->alp  * DEG,
       cubds[i]->hl[0],
       cubds[i]->hl[1],
       cubds[i]->hl[2],
       cubds[i]->closed_f[0],
       cubds[i]->closed_f[1],
       cubds[i]->rotate_f[0],
       cubds[i]->rotate_f[1],
       cubds[i]->rotate_f[2],
       cubds[i]->rotate_f[3]);
   }
 }


/* str_cyln_alloc

 Allocate memory for a cylinder structure (and its components).

 INPUT:
 None.

 OUTPUT:
 Pointer to a cylinder structure.

 */

 str_cyln* 
 str_cyln_alloc
 (
  void
 )
 {
   str_cyln* cyln = (str_cyln*) malloc( sizeof(str_cyln) );
   cyln->base = str_ellp_alloc();
   cyln->top  = str_ellp_alloc();
   cyln->M   = calloc_2d(3, 3, "cyln->M in str_cyln_alloc (structures.h)");
   cyln->alp = 0.0;
   cyln->closed_f[0] = 1;
   cyln->closed_f[1] = 1;
   cyln->rotate_f[0] = 0;
   cyln->rotate_f[1] = 0;
   cyln->rotate_f[2] = 0;
   cyln->rotate_f[3] = 0;

   return cyln;
 }

/* str_cyln_free

 Release memory of an cylinder structure.

 INPUT:
 cyln - Pointer to an cylinder structure pointer.

 OUTPUT:
 None. Releases memory block and set the cylinder pointer to NULL.

 */

 void 
 str_cyln_free
 (
  str_cyln **cyln
 )
 {
   str_ellp_free(&(*cyln)->base);
   str_ellp_free(&(*cyln)->top);
   free_2d(3, &(*cyln)->M);
   free(*cyln);
   *cyln = NULL;
 }

/* str_cyln_read

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

 */

 void
 str_cyln_read
 (
   FILE     *fi, 
   long int *fpos,
   double   *origin,
   double   *axis,
   double   *alpha,
   double   *radius,
   double   *height,
   double   *s_base,
   double   *s_top,
   int      *closed
 )
 {
   fseek(fi, *fpos, SEEK_SET);
   fscanf(fi, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d", 
     &origin[0], &origin[1], &origin[2], 
     &axis[0], &axis[1], alpha,
     radius, height,
     &s_base[0], &s_base[1],
     &s_top[0], &s_top[1],
     &closed[0], &closed[1]);
   *fpos = ftell(fi);

   axis[0] *= RAD;
   axis[1] *= RAD;
   *alpha  *= RAD;
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

/* str_cyln_setup

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

 */

 void 
 str_cyln_setup
 (
   str_cyln     *cyln,
   double       *origin,
   double const *axis,
   double const alpha,
   double const radius,
   double const height,
   double const *s_base,
   double const *s_top,
   int const    *closed
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

   if ( (axis[0] < TOLERANCE) && (axis[1] < TOLERANCE) && (alpha < TOLERANCE) )
   {
     cyln->u[0] = 0.0;
     cyln->u[1] = 0.0;
     cyln->u[2] = 1.0;
     cyln->M[0][0] = 1.0;
     cyln->M[1][1] = 1.0;
     cyln->M[2][2] = 1.0;
     cyln->rotate_f[0] = 0;
   }
   else
   {
     sph_to_cos_unit(cyln->u, cyln->s);
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
   cyln->dhmx[0] = -sin(s_base[0]);
   str_ellp_setup(cyln->base, origin, s_base, 0.0, open_r);
   cyln->rotate_f[1] = (s_base[0] < TOLERANCE)? 0 : 1;
   for (size_t i = 0; i < 3; i++)
   {
     origin[i] += cyln->u[i] * height;
   }
   open_r[0] = radius / cos(s_top[0]);
   cyln->dhmx[1] = -sin(s_top[0]);
//   if ( ABS(s_top[1] - M_PI) < TOLERANCE )
//     cyln->dhmx[1] *= -1.0;
   str_ellp_setup(cyln->top, origin, s_top, 0.0, open_r);
   cyln->rotate_f[2] = (s_top[0] < TOLERANCE)? 0 : 1;
   cyln->rotate_f[3] = (ABS(s_base[0] - s_top[0]) < TOLERANCE)? 0 : 1;
 }

/* str_cyln_printf, str_cyln_printf

 Print detailed information of a cylinder or a summary of multiple cylinders.

 INPUT:
 cyln - Pointer to the cylinder structure;

 OUTPUT:
 None. Prints the cylinder parameters.

 */

 void
 str_cyln_fprintf 
 (
   FILE *odv,
   str_cyln const *cyln,
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
   fprintf(odv, "%sClosed:    %d, %d\n", pre_1, cyln->closed_f[0], cyln->closed_f[1]);
   mat_fprintf(odv, (double const **) cyln->M, 3, 3, 1.0, "Rotation matrix", indent+1);
 }

 void
 str_cylns_printf
 (
  str_cyln const **cylns,
  size_t n
 )
 {
   printf("  Number of cylinders: %d\n", n);
   printf("                  type                     origin (m)       axis (º) "
     " alpha (º)  Radius (m)  Height (m)  Closed      Rotate\n");
   for(size_t i = 0; i < n; i++)
   {
     printf("  %02d  %s  % .2e,% .2e,% .2e  %6.2lf,%6.2lf     %6.2lf   "
       "% .2e   % .2e    %d, %d  %d, %d, %d, %d\n",
       i+1,
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

