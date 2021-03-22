
/*******************************************************************************
 accumulators.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 Version: 2021-01-12

 Grid and sectorial geometries are available to provide spatially resolved 
 integration. If no spatially resolved integration is requested, any accumulator
 will provide the total spatial integral. Here are provided functions to setup 
 the accumulators and to perform the integration.

 The accumulator structure, is declared and defined in the accumulators.h header
 file. The spatially resolved integral is tracked in the "grid" member of the
 accumulators. In essence each accumulator provides "bins" to accumulate 
 contributions (or weights) of rays exiting the system in a particular range of
 coordinates. The "grid" members are five-dimensional arrays, with the first
 dimension being the Sun zenith angle, the second for the single scattering
 albedo, the third for the bottom reflectance and the last two dimensions being
 the rows and columns. More information is provided in the accumulators.h header
 file.

 The accumulator geometry can be: (1) 'sectorial', which divide the annuli in 
 sectors for azimuthal asymmetry; or (2) 'grid', which will keep track of the 
 {x,y} position. An 'annular' accumulator is in general not appropriate since
 azimuthal effects are expected except when the Sun in is the zenith and under
 a sky radiance distribution for clear skies.

 The spatial extent of the accumulator is defined by the variable "acc_ext", 
 defining a radius (that is, extent is from -"acc_ext" to +"acc_ext"). Note that 
 regardless of extent of the accumulator, a final bin will be added between 
 "acc_ext" and INFINITY, such as to capture all contributions beyond the 
 accumulator spatially resolved extent. The resolutions of the accumulator are
 defined by the variables "acc_resy" and "acc_resx", and have different meaning 
 depending on geometry. For sectorial geometries the "acc_resy" identify the 
 metric resolution of the annuli, while "acc_resx" identify the angular 
 resolution (in degrees) of the sectors. Sectors will always cover the angular 
 extent from 0 to 360 degrees. For the grid simulations, normally "acc_resy" and 
 "acc_resx" have the same value and define the length of the grid cell in y and 
 x dimensions.

 Note that the resolution will have an impact on the performance, as a too high 
 resolution, specially in grid geometry will require more photons to be traced 
 to reduce statistical sampling noise.

 For the case to study the immediate impact of shadowing structures, the extent
 might be small, somewhat larger than the maximum horizontal dimension of the 
 shadowing structures (e.g., length of boat + distance to deployment of sensor).
 It should be somewhat larger since the system also considering shadow from the 
 in air structures, such that the projected shadow from direct and diffuse 
 components should be accounted for.

 When simulation has an axis of symmetry (e.g., symmetric structures centered
 on an axis to the Sun under clear skies) the user my opt for a post-processing
 that will reduce the Monte Carlo sample noise: that is, "fold" the grid on its
 axis of symmetry, and calculate the average. Since this processing is specific
 and can easily be done afterwards, its is not implemented in this code to keep
 it generic without too much complexity.

*******************************************************************************/

 #include <math.h>
 #include <stdlib.h>
 #include <string.h>
 #include <stdio.h>

 #include "config.h"
 #include "constants.h"
 #include "aux.h"
 #include "memory.h"
 #include "skyrad.h"
 #include "accumulators.h"

/* Allocate backward MC accumulator: *******************************************
 
 accm_b_alloc



*******************************************************************************/

 struct accumulator_bmc *
 accm_b_alloc
 ( void )
 {
   struct accumulator_bmc * accm = 
     (struct accumulator_bmc*) calloc(1, sizeof(struct accumulator_bmc));
   #ifdef CHCK_MEM_ALLOC
   if ( !accm )
   {
     printf("\nERROR: Could not allocate memory for accm in accm_alloc"
       " (accumulators.c)\n");
     exit(-1);
   }
   #endif // CHCK_MEM_ALLOC

   return accm;
 }

/* Deallocate backward MC accumulator: *****************************************
 
 accm_b_free



*******************************************************************************/

 void
 accm_b_free
 ( struct accumulator_bmc ** accm )
 {
   // Free spatial integral array:
   free_3d((*accm)->ns, (*accm)->nw0, &(*accm)->vector);

   #ifdef SPRES
   // Free spatially resolved integral array:
   free_1d(&(*accm)->xbrks);
   free_1d(&(*accm)->ybrks);
   free_5d((*accm)->ns, (*accm)->nw0, (*accm)->nbr, (*accm)->ny, &(*accm)->grid);
   #endif // SPRES

   free( *accm );
   *accm = NULL;
 }

/* Setup backward MC accumulator: *****************************************
 
 accm_setup



*******************************************************************************/

 void
 accm_setup
 (
   struct accumulator_bmc * accm,
   char const * accm_tp,
   int const sim_ns,
   int const iop_nw0,
   int const btt_nbr,
   double const acc_ext,
   double const acc_resy,
   double const acc_resx
 )
 {
   strncpy(accm->type, accm_tp, STRMXLEN);
   if ( strcmp (accm_tp, "sectorial") == 0 )
   {
/*     // Sectorial:
     accm_b_setup_sect(accm, sim_ns, iop_nw0, btt_nbr, acc_ext, acc_resy,
       acc_resx);
     return;
*/   }
   else if ( strcmp(accm_tp, "grid") == 0 )
   {
     // Grid:
     accm_b_setup_grid(accm, sim_ns, iop_nw0, btt_nbr, acc_ext, acc_resy,
       acc_resx);
     return;
   }
 }

/* Setup sectorial accumulator: ************************************************

 accm_b_setup_sect

 Setup the sectorial geometry accumulator.

 INPUT:
 sim_ns   - Number of illuminations (sun zeniths or sky radiance distributions);
 iop_nw0  - Number of single scattering albedos;
 btt_nbr  - Number of bottom reflectances;
 acc_ext  - Spatial extent for spatially resolved accumulator, meters, 
            (0,INFINITY).
 acc_resy - Spatial resolution for the annuli radius in meters, (0, INFINITY);
 acc_resx - Angular resolution of the sectors in radiance, (0, 2PI);

 OUTPUT:
 Returns a pointer to a initialized accm_bmc. 

*******************************************************************************/

 void
 accm_b_setup_sect
 (
   struct accumulator_bmc *accm,
   int const    sim_ns,
   int const    iop_nw0,
   int const    btt_nbr,
   double const acc_ext,
   double const acc_resy,
   double const acc_resx
 )
 {
   accm->ns  = sim_ns;
   accm->nw0 = iop_nw0;
   accm->nbr = btt_nbr;

   accm->fun.add = &accm_b_add_sect;

   // Allocate spatial integral array:
   accm->vector = calloc_3d(accm->ns, accm->nw0, accm->nbr, 
     "accm->vector in accm_b_setup_sect (accumulators.c)");

   #ifdef SPRES
   accm->ny    = floor((acc_ext - acc_resy) / acc_resy) + 2;
   accm->nx    = floor(K_2PI / acc_resx);
   accm->ybrks = calloc_1d(accm->ny + 1,
     "accm->ybrks in accm_b_setup_sect (accumulators.c)"); 
   accm->xbrks = calloc_1d(accm->nx + 1, 
     "accm->xbrks in accm_b_setup_sect (accumulators.c)");

   // Define break points in rows (y):
   accm->ybrks[0] = 0.0;
   accm->ybrks[1] = acc_resy / 2.0;
   accm->ybrks[accm->ny] = INFINITY;
   for (size_t ci = 2; ci < accm->ny; ci++)
   {
     accm->ybrks[ci] = accm->ybrks[ci - 1] + acc_resy;
   }

   // Define break points in columns (x):
   accm->xbrks[0] = 0.0;
   for (size_t ci = 1; ci <= accm->nx; ci++)
   {
     accm->xbrks[ci] = accm->xbrks[ci - 1] + acc_resx;
   }
    
   // Allocate spatially resolved integral array:
   accm->grid = calloc_5d(accm->ns, accm->nw0, accm->nbr, accm->ny, 
     accm->nx, "accm->grid in accm_b_setup_sect (accumulators.c)");

   // Calculate parameters for findBin function:
   accm->miny = 0.0;
   accm->ky_inv = accm->ny / accm->ybrks[accm->ny - 1];
   accm->minx = 0.0;
   accm->kx_inv = accm->nx / K_2PI;
   #endif // SPRES
 }

/* Add to sectorial accumulator: ***********************************************

 accm_b_add_sect

 Add contribution of rays to a sectorial accumulator.

 For the diffuse component, this function is called only once and loop is made 
 over all sky radiance distributions. For the direct component the situation 
 is more challenging as the ray has to be propagated to the surface for each
 sun zenith angle (with the local estimate method), so this function is called
 several times with different arguments.

 The dirf flag is available to simplify the condition when the accumulator is 
 for the direct component of irradiance, since propagation to sun is 
 calculated exactly.

 INPUT:
 acm   - Accumulator parameters
         Pointer to accumulator_bmc struct;
 skr   - Sky radiance distribution parameters
         Pointer to constant skyradiance struct;
 p     - Terminal point of the ray
         Range: (-Inf,Inf), meters
         [0]X, [1]Y, [2]Z
         Pointer to constant double;
 s     - Refracted spherical directions of the ray
         Range: [0,PI], [0,2PI], radians
         [0]Theta, [1]Phi, [2]Radius (always 1.0 and never used)
         Pointer to constant double;
 stks  - Diffuse Stokes vector
         [N w0][N bhr][4]: [][][0]I, [][][1]Q, [][][2]U, [][][3]V
         Pointer to pointer to pointer to constant double;
 scale - Arbitrary multiplier (e.g., Fresnel transmittance or probability to 
         Sun)
         Range: (0,Inf), unitless
         Constant double;
 dirf  - Flag to indicate if the contribution is for the direct component
         Values: (0) No, (1) Yes
         Constant int
 cs    - Index of the accumulator layer for the SZA in case of direct component
         Values: [0, N SZA]
         Constant int.

 OUTPUT:
 None. Updates the accumulator structure.

*******************************************************************************/

 void
 accm_b_add_sect
 (
   struct accumulator_bmc *accm,
   struct skyradiance const *skr,
   double const *p,
   double const *s,
   double const ***stks,
   double const scale,
   int const    dirf,
   int const    cs
 )
 {
   #ifdef SPRES
   // Find radius and azimuth:
   double azmt = 0.0;
   double radius = sqrt(p[0] * p[0] + p[1] * p[1]);
   if ( radius < TOLERANCE )
   {
     azmt = K_2PI;
   }
   else
   {
     azmt = acos(p[1] / radius);
     if (p[0] > 0.0)
     {
       azmt = K_2PI - azmt;
     }
   }

   int rid = findBin(radius, accm->miny, accm->ky_inv, accm->ny);
   int cid = findBin(azmt, accm->minx, accm->kx_inv, accm->nx);
   #endif // SPRES

   if ( dirf )
   {
     for (size_t cw = 0; cw < accm->nw0; cw++)
     {
       for (size_t cb = 0; cb < accm->nbr; cb++)
       {
         accm->vector[cs][cw][cb] += (stks[cw][cb][0] * scale);

         #ifdef SPRES
         accm->grid[cs][cw][cb][rid][cid] += (stks[cw][cb][0] * scale);
         #endif
       }
     }
   }
   else
   {
     int rids = findBin(s[0], skr->miny, skr->ky_inv, skr->ny);
     int cids = findBin(s[1], skr->minx, skr->kx_inv, skr->nx);

     for (size_t cz = 0; cz < skr->ns; cz++)
     {
       double skw = skr->grid[cz][rids][cids];
       for (size_t cw = 0; cw < accm->nw0; cw++)
       {
         for (size_t cb = 0; cb < accm->nbr; cb++)
         {
           accm->vector[cz][cw][cb] += (stks[cw][cb][0] * scale * skw);

           #ifdef SPRES
           accm->grid[cz][cw][cb][rid][cid] += (stks[cw][cb][0] * scale * skw);
           #endif
         } // Loop: accm->nbr
       } // Loop: accm->nw0
     } // Loop: accm->ns
   } // Logical: Direct flag on?
 }


/* Setup a grid accumulator: ***************************************************

 accm_b_setup_grid

 Setup the grid geometry accumulator.

 INPUT:
 sim_ns   - Number of illuminations (sun zeniths or sky radiance distributions);
 iop_nw0  - Number of single scattering albedos;
 btt_nbr  - Number of bottom reflectances;
 acc_ext  - Extent in meters, (0, INFINITY);
 acc_resy - Spatial resolution in the Y axis in meters, (0, INFINITY);
 acc_resx - Spatial resolution in the X axis in meters, (0, INFINITY);

 OUTPUT:
 Returns a pointer to a initialized accm_bmc. 

*******************************************************************************/

 void
 accm_b_setup_grid
 (
   struct accumulator_bmc * accm,
   int const sim_ns,
   int const iop_nw0,
   int const btt_nbr,
   double const acc_ext,
   double const acc_resy,
   double const acc_resx
 )
 {
   accm->ns  = sim_ns;
   accm->nw0 = iop_nw0;
   accm->nbr = btt_nbr;

   accm->fun.add = &accm_b_add_grid;

   // Allocate spatial integral array:
   accm->vector = calloc_3d (accm->ns, accm->nw0, accm->nbr, 
     "accm->vector in accm_b_setup_sect (accumulators.c)");

   #ifdef SPRES
   accm->ny = 2 * floor(((acc_ext - acc_resy / 2) / acc_resy + 2)) - 1;
   accm->nx = 2 * floor(((acc_ext - acc_resx / 2) / acc_resx + 2)) - 1;
   accm->ybrks = calloc_1d(accm->ny + 1, 
     "accm->ybrks in accm_b_setup_sect (accumulators.c)"); 
   accm->xbrks = calloc_1d(accm->nx + 1, 
     "accm->xbrks in accm_b_setup_sect (accumulators.c)");

   // Define break points in rows (y):
   accm->ybrks[0] = -INFINITY;
   accm->ybrks[(accm->ny + 1) / 2] = acc_resy / 2;
   accm->ybrks[accm->ny] = INFINITY;
   for (size_t ci = ((accm->ny + 1) / 2) + 1; ci < accm->ny; ci++)
   {
     accm->ybrks[ci] = accm->ybrks[ci - 1] + acc_resy;
   }
   for (size_t ci = ((accm->ny + 1) / 2) - 1; ci > 0; ci--)
   {
     accm->ybrks[ci] = accm->ybrks[ci + 1] - acc_resy;
   }

   // Define break points in columns (x):
   accm->xbrks[0] = -INFINITY;
   accm->xbrks[(accm->nx + 1) / 2] = acc_resx / 2;
   accm->xbrks[accm->nx] = INFINITY;
   for (size_t ci = ((accm->nx + 1) / 2) + 1; ci < accm->nx; ci++)
   {
     accm->xbrks[ci] = accm->xbrks[ci - 1] + acc_resx;
   }
   for (size_t ci = ((accm->nx + 1) / 2) - 1; ci > 0; ci--)
   {
     accm->xbrks[ci] = accm->xbrks[ci + 1] - acc_resx;
   }

   // Allocate spatially resolved integral array:
   accm->grid = calloc_5d (accm->ns, accm->nw0, accm->nbr, accm->ny, 
     accm->nx, "accm->grid in accm_b_setup_sect (acm.geom.c)");

   // Calculate parameters for findBin function:
   accm->miny = accm->ybrks[1] - acc_resy;
   accm->ky_inv = accm->ny / (accm->ybrks[accm->ny - 1] + 
     acc_resy - accm->miny);
   accm->minx = accm->xbrks[1] - acc_resx;
   accm->kx_inv = accm->nx / (accm->xbrks[accm->nx - 1] + 
     acc_resx - accm->minx);
   #endif // SPRES
 }

/* Add to regular grid accumulator: ********************************************

 accm_b_add_grid

 Add contribution of rays to a regular grid accumulator.

 For the diffuse component, this function is called only once and loop is made 
 over all sky radiance distributions. For the direct component the situation 
 is more challenging as the ray has to be propagated to the surface for each
 sun zenith angle (with the local estimate method), so this function is called
 several times with different arguments.

 The dirf flag is available to simplify the condition when the accumulator is 
 for the direct component of irradiance, since propagation to sun is 
 calculated exactly.

 INPUT:
 acm   - Accumulator parameters
         Pointer to accumulator_bmc struct;
 skr   - Sky radiance distribution parameters
         Pointer to constant skyradiance struct;
 p     - Terminal point of the ray
         Range: (-Inf,Inf), meters
         [0]X, [1]Y, [2]Z
         Pointer to constant double;
 s     - Refracted spherical directions of the ray
         Range: [0,PI], [0,2PI], radians
         [0]Theta, [1]Phi, [2]Radius (always 1.0 and never used)
         Pointer to constant double;
 stks  - Diffuse Stokes vector
         [N w0][N bhr][4]: [][][0]I, [][][1]Q, [][][2]U, [][][3]V
         Pointer to pointer to pointer to constant double;
 scale - Arbitrary multiplier (e.g., Fresnel transmittance or probability to 
         Sun)
         Range: (0,Inf), unitless
         Constant double;
 dirf  - Flag to indicate if the contribution is for the direct component
         Values: (0) No, (1) Yes
         Constant int
 cs    - Index of the accumulator layer for the SZA in case of direct component
         Values: [0, N SZA]
         Constant int.

 OUTPUT:
 None. Updates the accumulator structure.

*******************************************************************************/


 void 
 accm_b_add_grid
 (
   struct accumulator_bmc * accm,
   struct skyradiance const * skr,
   double const * p,
   double const * s,
   double const *** stks,
   double const scale,
   int const dirf,
   int const cs
 )
 {
   #ifdef SPRES
   // int rid = findBin(p[1], accm->miny, accm->ky_inv, accm->ny);
   // int cid = findBin(p[0], accm->minx, accm->kx_inv, accm->nx);
   int rid, cid;
   FINDBIN(p[1], accm->miny, accm->ky_inv, accm->ny, rid);
   FINDBIN(p[0], accm->minx, accm->kx_inv, accm->nx, cid);
   #endif

   if ( dirf )
   {
     for (size_t cw = 0; cw < accm->nw0; cw++)
     {
       for (size_t cb = 0; cb < accm->nbr; cb++)
       {
         accm->vector[cs][cw][cb] += (stks[cw][cb][0] * scale);

         #ifdef SPRES
         accm->grid[cs][cw][cb][rid][cid] += (stks[cw][cb][0] * scale);
         #endif
       }
     }
   } 
   else
   {
     // int rids = findBin(s[0], skr->miny, skr->ky_inv, skr->ny);
     // int cids = findBin(s[1], skr->minx, skr->kx_inv, skr->nx);
     int rids, cids;
     FINDBIN(s[0], skr->miny, skr->ky_inv, skr->ny, rids);
     FINDBIN(s[1], skr->minx, skr->kx_inv, skr->nx, cids);

     for(size_t cz = 0; cz < skr->ns; cz++)
     {
       double skw = skr->grid[cz][rids][cids];
       for (size_t cw = 0; cw < accm->nw0; cw++)
       {
         for (size_t cb = 0; cb < accm->nbr; cb++)
         {
           accm->vector[cz][cw][cb] += (stks[cw][cb][0] * scale * skw);

           #ifdef SPRES
           accm->grid[cz][cw][cb][rid][cid] += (stks[cw][cb][0] * scale * skw);
           #endif
         } // Loop: accm->nbr
       } // Loop: accm->nw0
     } // Loop: accm->ns
   } // Logical: Direct flag on?
 }

/* Normalize accumulator: ******************************************************
 
 accm_norm

 Normalize the accumulator. Normalization factors are dependent on simulation 
 conditions and are calculated inside the MC solver.

 Note: the normalization factor depends on each sensor type (L, E, E0) and at 
 the curent version are specified directly in the MC code. It will be in the
 source code for sources (src.c) in future versions.

 INPUT:
 accm  - Pointer to accumulator structure.
 f0    - Incident irradiance (normally set to 1 to get reflectance);
 normf - Sensor dependent normalization factor.

 OUTPUT:
 None. Updates the values pointed by accm.
 
*******************************************************************************/

 void
 accm_b_norm
 (
   struct accumulator_bmc *accm,
   double const f0,
   double const *normf
 )
 {
   for (size_t cz = 0; cz < accm->ns; cz++)
   {
     for (size_t cw = 0; cw < accm->nw0; cw++)
     {
       for (size_t cb = 0; cb < accm->nbr; cb++)
       {
         accm->vector[cz][cw][cb] *= (f0 * normf[cz]);

         #ifdef SPRES
         for (size_t cr = 0; cr < accm->ny; cr++)
         {
           for (size_t cc = 0; cc < accm->nx; cc++)
           {
             accm->grid[cz][cw][cb][cr][cc] *= (f0 * normf[cz]);
           }
         }
         #endif // SPRES
       }
     }
   }
 }

/* Add value to a single cell: *************************************************

 accm_add_k

 Add a fixed value to the accumulator. This function is used to add the direct
 transmission contribution when the Sun is in the FOV of the sensor.

 INPUT:
 accm - Pointer to accumulator structure;
 val  - Value to be added;
 cs   - Index of the sun zenith angle layer;
 cr   - Index of the row for the spatially resolved accumulator;
 cc   - Index of the column for the spatially resolved accumulator.

 OUTPUT:
 None. Updates the values pointed by accm.

*******************************************************************************/

 void 
 accm_b_add_k
 (
   struct accumulator_bmc *accm,
   double const val,
   int const cs,
   int const cr,
   int const cc
 )
 {
   for (size_t cw = 0; cw < accm->nw0; cw++)
   {
     for (size_t cb = 0; cb < accm->nbr; cb++)
     {
       accm->vector[cs][cw][cb] += val;

       #ifdef SPRES
       accm->grid[cs][cw][cb][cr][cc] += val;
       #endif
     }
   }
 }

/* Sum accumulators: ***********************************************************

 accm_sum

 Sum accumulators. This function is used to calculate the average by summing 
 over the results of the subsets.

 INPUT:
 accm_sum - Pointer to accumulator receiving values;
 accm_in  - Pointer to the accumulator from simulation;
 scale    - Multiplier (1/number of subdivisions).

 OUTPUT:
 None. Updates the values pointed by accm_sum.

*******************************************************************************/

 void
 accm_b_sum
 (
   struct accumulator_bmc *accm_sum,
   struct accumulator_bmc const *accm_in,
   double const scale
 )
 {
   for (size_t cz = 0; cz < accm_in->ns; cz++)
   {
     for (size_t cw = 0; cw < accm_in->nw0; cw++)
     {
       for (size_t cb = 0; cb < accm_in->nbr; cb++)
       {
         accm_sum->vector[cz][cw][cb] += accm_in->vector[cz][cw][cb] * scale;

         #ifdef SPRES
         for (size_t cr = 0; cr < accm_in->ny; cr++)
         {
           for (size_t cc = 0; cc < accm_in->nx; cc++)
           {
             accm_sum->grid[cz][cw][cb][cr][cc] += 
               accm_in->grid[cz][cw][cb][cr][cc] * scale;
           }
         }
         #endif // SPRES
       }
     }
   }
 }

/* Reset accumulator: **********************************************************

 accm_reset

 Reset accumulator to 0.0.
 
 INPUT:
 accm - Pointer to accumulator.

 OUTPUT:
 None. Updates values of the accumulator.

*******************************************************************************/

 void
 accm_b_reset
 ( struct accumulator_bmc *accm )
 {
   for (size_t cz = 0; cz < accm->ns; cz++)
   {
     for (size_t cw = 0; cw < accm->nw0; cw++)
     {
       for (size_t cb = 0; cb < accm->nbr; cb++)
       {
         accm->vector[cz][cw][cb] = 0.0;

         #ifdef SPRES
         for (size_t cr = 0; cr < accm->ny; cr++)
         {
           for (size_t cc = 0; cc < accm->nx; cc++)
           {
             accm->grid[cz][cw][cb][cr][cc] = 0.0;
           }
         }
         #endif // SPRES
       }
     }
   }
 }

/* Write the grid or vector member of an accumulator: **************************

 accm_b_write_grid,  accm_b_write_vec

 INPUT:
 accm    - Pointer to accumulator
 sim_ns  - Number of Sun zenith angles
 ofbn    - Base file name
 sufx    - Suffix
 iop_w0  - Pointer to single scattering albedos
 btt_bhr - Pointer to bottom reflectance
 sim_sza - Pointer to Sun zenith angles

*******************************************************************************/

 void
 accm_b_write_grid
 (
   struct accumulator_bmc *accm,
   int const    sim_ns,
   char const   *ofbn,
   char const   *sufx,
   double const *iop_w0,
   double const *btt_bh,
   double const *sim_sza
 )
 {
   char wf[12], bf[12], sz[12];
   char ofn[200];
   FILE *fpo;

   for (size_t cs = 0; cs < accm->ns; cs++)
   {
     if (cs < sim_ns)
     {
       sprintf(sz, "_S%03d", (int) ceil(sim_sza[cs] * DEG));
     }
     else if (cs == sim_ns)
     {
       sprintf(sz, "_Isotropic");
     } 
     else
     {
       sprintf(sz, "_Cardioidal");
     }
     for (size_t cw = 0; cw < accm->nw0; cw++)
     {
       sprintf(wf, "_W%05.2f", iop_w0[cw]);
       for (size_t cb = 0; cb < accm->nbr; cb++)
       {
         sprintf(bf, "_B%05.3f", btt_bh[cb]);
         strncpy(ofn, ofbn, STRMXLEN);
         strcat(ofn, sufx);
         strcat(ofn, sz);
         strcat(ofn, wf);
         strcat(ofn, bf);
         strcat(ofn, ".txt");

         fpo = fopen (ofn, "w");
         if (fpo == NULL)
         {
           printf("\nERROR: Failed to create output file %s\n", ofn);
           exit(-3);
         }
         for (size_t cr = 0; cr < accm->ny; cr++)
         {
           for (size_t cc = 0; cc < (accm->nx - 1); cc++)
           {
             fprintf(fpo, "%.6E ", accm->grid[cs][cw][cb][cr][cc]);
           }
           fprintf(fpo, "%.6E\n", accm->grid[cs][cw][cb][cr][accm->nx - 1]);
         }
         fclose (fpo);
       }
     }
   }
 }


 void
 accm_b_write_vec
 (
   struct accumulator_bmc *accm,	// Pointer to accumulator 
   int const    sim_ns,		// Number of Sun zenith angles
   char const   *ofbn,		// Base file name
   char const   *sufx,		// Suffix
   double const *iop_w0,	// Pointer to single scattering albedos
   double const *btt_bhr,	// Pointer to bottom reflectance
   double const *sim_sza	// Pointer to Sun zenith angles
 )
 {
   char wf[10], bf[10];
   char ofn[200];
   FILE *fpo;

   for (size_t cb = 0; cb < accm->nbr; cb++)
   {
     sprintf(bf, "_B%05.3f", btt_bhr[cb]);
     strncpy(ofn, ofbn, STRMXLEN);
     strcat(ofn, sufx);
     strcat(ofn, bf);
     strcat(ofn, ".txt");

     fpo = fopen (ofn, "w");
     if (fpo == NULL)
     {
       printf("\nERROR: Failed to create output file %s\n", ofn);
       exit(-3);
     }

     fprintf(fpo, "\t");
     for (size_t cw = 0; cw < (accm->nw0 - 1); cw++)
     {
       fprintf(fpo, "W_%6.4f ", iop_w0[cw]);
     }
     fprintf(fpo, "W_%6.4f\n", iop_w0[(accm->nw0 - 1)]);
     for (size_t cz = 0; cz < accm->ns; cz++)
     {
       if (cz < sim_ns)
       {
         fprintf(fpo, "I_%04.1f\t", sim_sza[cz] * DEG);
       } else if (cz == sim_ns) {
         fprintf(fpo, "UNIF\t");
       } else {
         fprintf(fpo, "CARD\t");
       }
       for (size_t cw = 0; cw < (accm->nw0 - 1); cw++)
       {
         fprintf(fpo, "%.6E ", accm->vector[cz][cw][cb]);
       }
       fprintf(fpo, "%.6E\n", accm->vector[cz][(accm->nw0 - 1)][cb]);
     }
     fclose (fpo);
   }
 }

