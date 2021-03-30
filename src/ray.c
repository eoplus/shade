
/*******************************************************************************
 ray.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 Version: 1.6
 Date: 2021-03-25
 License: GPL-3.0

 Functions to allocate, free, set initial conditions, update position, update 
 direction and print a light_ray datatype.

*******************************************************************************/

 #include <stdio.h>
 #include <string.h>
 #include <stdlib.h>
 #include <math.h>

 #include "config.h"
 #include "constants.h"
 #include "aux.h"
 #include "memory.h"
 #include "geometry.h"
 #include "rotation.h"
 #include "ray.h"

/* Allocate ray: ***************************************************************

 ray_alloc

 Allocates internal components of a light_ray struct. Details of the light_ray 
 data type are provided in the ray.h header file.

 INPUT:
 real_f  - Flag to indicate if it is a real ray (1) or virtual (0).
 iop_nw0 - Number of single scatteirng albedos of the simulation.
 btt_nbr - Number of bottom bi-hemispherical reflectances of the simulation.

 OUTPUT:
 None. Memmory is allocated to members of a light_ray.

*******************************************************************************/

 void
 ray_alloc
 (
   int const real_f,
   int const iop_nw0,
   int const btt_nbr,
   struct light_ray * ray
 )
 {
   ray->nw0 = 0;
   ray->nbr = 0;
   ray->nsdw_f = 1;
   ray->real_f = 0;

   if ( real_f )
   {
     ray->stks = calloc_3d(iop_nw0, btt_nbr, 4, 
       "ray->stks in ray_alloc (ray.c)");
     ray->real_f = 1;
     ray->nw0  = iop_nw0;
     ray->nbr  = btt_nbr;
   }
 }

/* Deallocate ray: *************************************************************

 ray_free

 Free the memory block pointed members of a light_ray struct and set the 
 pointer(s) value to NULL.

 INPUT:
 ray - Ray properties
       Pointer to light_ray struct; 

 OUTPUT:
 None. Free memory and update light_ray member pointer(s) value to NULL.

*******************************************************************************/

 void
 ray_free
 (
   struct light_ray * ray
 )
 {
   if ( ray->real_f )
     free_3d(ray->nw0, ray->nbr, &(ray)->stks);
 }


/* Print ray: ******************************************************************

 ray_fprintf

 Print the properties of the light_ray with formatted text.

 INPUT:
 odv    - Output device. Should be a file stream (e.g. 'stdout' to print to
          default);
 ray    - Pointer to a light_ray.
 indent - Number of indentations at the header level (each indent = 2 spaces)
          Constant int;

 OUTPUT:
 Print components.

*******************************************************************************/

 void
 ray_fprintf
 (
   FILE *odv,
   struct light_ray *ray,
   int const indent 
 )
 {
   // Set indentation:
   char pre_0[STRMXLEN] = "";
   char pre_1[STRMXLEN] = "";
   for (int i = 0; i < indent; i++)
     strcat(pre_0, "  ");
   strncpy(pre_1, pre_0, STRMXLEN);
   strcat(pre_1, "  ");

   double l = sqrt(
     pow(ray->b[0] - ray->a[0], 2.0) + 
     pow(ray->b[1] - ray->a[1], 2.0) + 
     pow(ray->b[2] - ray->a[2], 2.0) );

   fprintf(odv, "%sRay parameters:\n", pre_0);
   fprintf(odv, "%sOrigin:     % .2e, % .2e, % .2e (m)\n",
     pre_1,
     ray->a[0], 
     ray->a[1], 
     ray->a[2]);
   fprintf(odv, "%sPosition:   % .2e, % .2e, % .2e (m)\n", 
     pre_1,
     ray->b[0], 
     ray->b[1], 
     ray->b[2]);
   fprintf(odv, "%sPathlength: % .2e (m)\n", pre_1, l);
   fprintf(odv, "%sDirection:   %6.2fº, %6.2fº (% .2e, % .2e, % .2e)\n", 
     pre_1,
     ray->s[0] * DEG, 
     ray->s[1] * DEG, 
     ray->u[0], 
     ray->u[1], 
     ray->u[2]);
   fprintf(odv, "%sActive:  %d\n", pre_1, ray->nsdw_f);
   fprintf(odv, "%sVirtual: %d\n", pre_1, !ray->real_f);

   if ( ray->real_f )
   {
     fprintf(odv, "%sStokes (only first vector is shown):\n", pre_1);
     fprintf(odv, "%s| I | % .3e\n", pre_1, ray->stks[0][0][0]);
     fprintf(odv, "%s| Q | % .3e\n", pre_1, ray->stks[0][0][1]);
     fprintf(odv, "%s| U | % .3e\n", pre_1, ray->stks[0][0][2]);
     fprintf(odv, "%s| V | % .3e\n", pre_1, ray->stks[0][0][3]);
   }
 }

