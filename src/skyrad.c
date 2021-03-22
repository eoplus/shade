
/*******************************************************************************
 skyrad.c

 Version: 2021-01-11

 In the current version, it reads in external sky radiance distributions 
 (i.e., from the disk). Data is allocated in three dimensional arrays with the 
 following ordering:
 
 [sim_ns][NPl][NAz]

 where sim_ns is the number of sun zenith angles, NPl is the number of polar 
 angles and NAz the number of azimuthal angles in the sky radiance data. NPl 
 and NAz are defined by the resolution of the sky radiance data. The 
 simulations generating those distributions are external. Included with this 
 version of the code are the simulations of Castagna et al. (2019).

 Data should present the relative contribution to the diffuse irradiance. 
 Therefore, the per solid angle radiance is saled by cos(theta) * sin(theta), 
 and normalized to integrate to 1.

 References:
 Castagna, A.; Johnson, C. B.; Voss, K.; Dierssen, H. M.; Patrick, H.; Germer, 
 T. A.; Sabbe, K.; Vyverman, W. 2019. Uncertainty in global downwelling plane 
 irradiance estimates from sintered polytetrafluoroethylene plaque radiance 
 measurements. Applied Optics 58, 16, 4497. DOI: 10.1364/AO.58.004497

*******************************************************************************/

 #include <stdio.h>
 #include <stdlib.h>
 #include <math.h>

 #include "config.h"
 #include "memory.h"
 #include "skyrad.h"

/* Allocate and deallocate a sky radiance struct: ******************************

 skr_alloc, skr_free

 Allocate and deallocate memory for a skyradiance struct. When allocating, 
 values are initialized to zero and when deallocating, the pointer is set to
 NULL.

 INPUT:
 skr - Sky radiance parameters
       Pointer to pointer to sky radiance struct;

 OUTPUT:
 The adress of a sky radiance struct memory block or none, when deallocating.
 
*******************************************************************************/

 struct skyradiance*
 skr_alloc
 ( void )
 { 
   struct skyradiance *skr = 
     (struct skyradiance*) calloc(1, sizeof(struct skyradiance));
   #ifdef CHCK_MEM_ALLOC
   if ( !skr )
   {
     printf("\nERROR: Could not allocate memory for skr in skr_alloc"
       " (skyrad.c)\n");
     exit(-1);
   }
   #endif // CHCK_MEM_ALLOC

   return skr;
 }

 void
 skr_free
 ( struct skyradiance **skr )
 {
   free_1d( &(*skr)->ybrks );
   free_1d( &(*skr)->xbrks );
   free_1d( &(*skr)->vector );
   free_3d( (*skr)->ns, (*skr)->ny, &(*skr)->grid );
   free( *skr );
   *skr = NULL;
 }

/* Setup a sky radiance distribution: ******************************************

 skr_setup

 Setup a skyradiance struct.

 INPUT:
 skr      - Sky radiance parameters
            Pointer to skyradiance struct;
 sim_ns   - Number of Sun zenith angles
            Range: (0,Inf), unitless
            Constant int;
 skr_ny   - Number of rows in Sky radiance data
            Range: (0,Inf), unitless
            Constant int;
 skr_nx   - Number of columns in Sky radiance data
            Range: (0,Inf), unitless
            Constant int;
 skr_resy - Angular resolution of Sky radiance rows
            Range: (0,Inf), degree
            Constant double;
 skr_resx - Angular resolution of Sky radiance columns
            Range: (0,Inf), degree
            Constant double;
 skr_fls  - Path to sky radiance files
            Pointer to pointer to constant char;
 sim_sza  - Sun zenith angles
            Range: (0,Inf), meters
            Constant double;

 OUTPUT:
 None. Updates the values pointed by skr.
 
*******************************************************************************/

 void
 skr_setup
 (
   struct skyradiance *skr,
   int const    sim_ns,
   int const    skr_ny,
   int const    skr_nx,
   double const skr_resy,
   double const skr_resx,
   char const   **skr_fls,
   double const *sim_sza
 )
 {
   skr->ns = sim_ns + 2;
   skr->ny = skr_ny;
   skr->nx = skr_nx;
   skr->ybrks = calloc_1d(skr_ny + 1, "skr->ybrks in skr_setup (skyrad.c)");
   skr->xbrks = calloc_1d(skr_nx + 1, "skr->xbrks in skr_setup (skyrad.c)");
   skr->grid  = calloc_3d (skr->ns, skr->ny, skr->nx, 
     "skr->grid in skr_setup (skyrad.c)");
   skr->minx = 0.0;
   skr->miny = 0.0;
   skr->kx_inv = 1.0 / skr_resx;
   skr->ky_inv = 1.0 / skr_resy;

   skr->ybrks[0] = 0.0;
   for (size_t i = 1; i < (skr_ny + 1); i++)
     skr->ybrks[i] = skr->ybrks[i - 1] + skr_resy;

   skr->xbrks[0] = 0.0;
   for(size_t i = 1; i < (skr_nx + 1); i++)
     skr->xbrks[i] = skr->xbrks[i - 1] + skr_resx;

   skr_read(skr, skr_fls);
 }

/* Read sky radiance data from disk: *******************************************

 skr_read

 Read values opf sky radiance from disk.

 INPUT:
 skr     - Sky radiance parameters
           Pointer to sky radiance struct;
 skr_fls - File paths to sky radiance files
           Pointer to pointer to constant char;

 OUTPUT:
 None. Updates the values pointed by skr.
 
*******************************************************************************/

 void
 skr_read
 (
   struct skyradiance *skr,
   char const **skr_fls
 )
 {
   int NPl, NAz, skrid;
   int cs, cy, cx;
   FILE* fpi;

   for (size_t cs = 0; cs < skr->ns; cs++) 
   {
     fpi = fopen(skr_fls[cs], "r");
     if ( !fpi )
     {
       printf("\nERROR: Failed to open input file %s\n", skr_fls[cs]);
       exit(-1);
     }

     for (size_t cy = 0; cy < skr->ny; cy++)
     {
       for (size_t cx = 0; cx < skr->nx; cx++)
       {
         fscanf(fpi, "%lf", &skr->grid[cs][cy][cx]);
       }
     }
   }
 }

