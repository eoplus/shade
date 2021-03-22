
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

 Allocate memory and check allocation for an accumulator_bmc struct.

 INPUT:
 None.

 OUTPUT:
 Pointer to accumulator_bmc struct.

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

   // Those will be allocated and setup as necessary by the setup function:
   accm->xbrks = NULL;
   accm->ybrks = NULL;
   accm->vector = NULL;
   accm->grid = NULL;

   return accm;
 }

/* Deallocate backward MC accumulator: *****************************************
 
 accm_b_free

 Release memory of an accumulator_bmc struct.

 INPUT:
 accm - Accumulator parameters
        Pointer to pointer to accumulator_bmc struct.

 OUTPUT:
 None. Releases memory and updates pointer to NULL.

*******************************************************************************/

 void
 accm_b_free
 ( struct accumulator_bmc ** accm )
 {
   // Free spatial integral array:
   free_3d( (*accm)->ns, (*accm)->nw0, &(*accm)->vector );

   #ifdef SPATIALLY_RESOLVED
   // Free spatially resolved integral array:
   free_1d( &(*accm)->xbrks );
   free_1d( &(*accm)->ybrks );
   free_5d( (*accm)->ns, (*accm)->nw0, (*accm)->nbr, (*accm)->ny, 
     &(*accm)->grid );
   #endif // SPATIALLY_RESOLVED

   free( *accm );
   *accm = NULL;
 }

/* Setup MC accumulator: *******************************************************
 
 accm_setup

 This is a wrapper function that will call the appropriate setup function for
 the accumulator type.

 INPUT:
 acm      - Accumulator parameters
            Pointer to accumulator_bmc struct;
 accm_tp  - Accumulator type
            Values: "sectorial", "grid"
            Pointer to constant char;
 sim_ns   - Number of illuminations
            Integer range: [1,Inf), unitless
            Constant int;
 iop_nw0  - Number of single scattering albedos
            Integer range: [1,Inf), unitless
            Constant int;
 btt_nbr  - Number of bottom bi-hemispherical reflectances
            Integer range: [1,Inf), unitless
            Constant int;
 acc_ext  - Spatial extent for the spatially resolved accumulator
            Range: (0,Inf), meters
            Constant double;
 acc_resy - Y-axis resolution
            Range: (0,Inf), meters
            Constant double;
 acc_resx - X-axis resolution
            Range: (0,Inf), meters
            Constant double;

 OUTPUT:
 None. Updates the accumulator pointed by accm.

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
     // Sectorial:
     accm_b_setup_sect(accm, sim_ns, iop_nw0, btt_nbr, acc_ext, acc_resy,
       acc_resx);
     return;
   }
   else if ( strcmp(accm_tp, "grid") == 0 )
   {
     // Regular grid:
     accm_b_setup_grid(accm, sim_ns, iop_nw0, btt_nbr, acc_ext, acc_resy,
       acc_resx);
     return;
   }
 }

/* Setup sectorial accumulator: ************************************************

 accm_b_setup_sect

 Setup the sectorial geometry accumulator. For details on INPUT and OUTPUT, see
 the function 'accm_setup'.

*******************************************************************************/

 void
 accm_b_setup_sect
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

   accm->fun.add = &accm_b_add_sect;

   // Allocate spatial integral array:
   accm->vector = calloc_3d(accm->ns, accm->nw0, accm->nbr, 
     "accm->vector in accm_b_setup_sect (accumulators.c)");

   #ifdef SPATIALLY_RESOLVED
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
   #endif //  SPATIALLY_RESOLVED
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

 For details of INPUT and OUTPUT, see the generic 'accm_b_add' in the
 accumulators.h header file.

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
   #ifdef SPATIALLY_RESOLVED
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

   int rid, cid;
   FINDBIN(radius, accm->miny, accm->ky_inv, accm->ny, rid);
   FINDBIN(azmt, accm->minx, accm->kx_inv, accm->nx, cid);
   #endif // SPATIALLY_RESOLVED

   if ( dirf )
   {
     for (size_t cw = 0; cw < accm->nw0; cw++)
     {
       for (size_t cb = 0; cb < accm->nbr; cb++)
       {
         accm->vector[cs][cw][cb] += (stks[cw][cb][0] * scale);

         #ifdef SPATIALLY_RESOLVED
         accm->grid[cs][cw][cb][rid][cid] += (stks[cw][cb][0] * scale);
         #endif // SPATIALLY_RESOLVED
       }
     }
   }
   else
   {
     int rids, cids;
     FINDBIN(s[0], skr->miny, skr->ky_inv, skr->ny, rids);
     FINDBIN(s[1], skr->minx, skr->kx_inv, skr->nx, cids);

     for (size_t cz = 0; cz < skr->ns; cz++)
     {
       double skw = skr->grid[cz][rids][cids];
       for (size_t cw = 0; cw < accm->nw0; cw++)
       {
         for (size_t cb = 0; cb < accm->nbr; cb++)
         {
           accm->vector[cz][cw][cb] += (stks[cw][cb][0] * scale * skw);

           #ifdef SPATIALLY_RESOLVED
           accm->grid[cz][cw][cb][rid][cid] += (stks[cw][cb][0] * scale * skw);
           #endif // SPATIALLY_RESOLVED
         } // Loop: accm->nbr
       } // Loop: accm->nw0
     } // Loop: accm->ns
   } // Logical: Direct flag on?
 }


/* Setup a grid accumulator: ***************************************************

 accm_b_setup_grid

 Setup the regular grid geometry accumulator. For details on INPUT and OUTPUT,
 see the function 'accm_setup'.

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

   #ifdef SPATIALLY_RESOLVED
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
   #endif // SPATIALLY_RESOLVED
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

 For details of INPUT and OUTPUT, see the generic 'accm_b_add' in the
 accumulators.h header file.

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
   #ifdef SPATIALLY_RESOLVED
   int rid, cid;
   FINDBIN(p[1], accm->miny, accm->ky_inv, accm->ny, rid);
   FINDBIN(p[0], accm->minx, accm->kx_inv, accm->nx, cid);
   #endif // SPATIALLY_RESOLVED

   if ( dirf )
   {
     for (size_t cw = 0; cw < accm->nw0; cw++)
     {
       for (size_t cb = 0; cb < accm->nbr; cb++)
       {
         accm->vector[cs][cw][cb] += (stks[cw][cb][0] * scale);

         #ifdef SPATIALLY_RESOLVED
         accm->grid[cs][cw][cb][rid][cid] += (stks[cw][cb][0] * scale);
         #endif // SPATIALLY_RESOLVED
       }
     }
   } 
   else
   {
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

           #ifdef SPATIALLY_RESOLVED
           accm->grid[cz][cw][cb][rid][cid] += (stks[cw][cb][0] * scale * skw);
           #endif // SPATIALLY_RESOLVED
         } // Loop: accm->nbr
       } // Loop: accm->nw0
     } // Loop: accm->ns
   } // Logical: Direct flag on?
 }

/* Normalize backward accumulator: *********************************************
 
 accm_b_norm

 Normalize the backward accumulator. Normalization factors are dependent on
 simulation conditions and are calculated inside the MC solver (function bmc in
 the mc_backward.c source file).

 INPUT:
 acm   - Accumulator parameters
         Pointer to accumulator_bmc struct;
 f0    - Incident irradiance (normally set to 1.0 to get reflectance)
         Range: (0,Inf), W/m^2
         Constant double;
 normf - Sensor dependent normalization factor
         Pointer to constant double.

 OUTPUT:
 None. Updates the values pointed by accm.
 
*******************************************************************************/

 void
 accm_b_norm
 (
   struct accumulator_bmc * accm,
   double const f0,
   double const * normf
 )
 {
   for (size_t cz = 0; cz < accm->ns; cz++)
   {
     for (size_t cw = 0; cw < accm->nw0; cw++)
     {
       for (size_t cb = 0; cb < accm->nbr; cb++)
       {
         accm->vector[cz][cw][cb] *= (f0 * normf[cz]);

         #ifdef SPATIALLY_RESOLVED
         for (size_t cr = 0; cr < accm->ny; cr++)
         {
           for (size_t cc = 0; cc < accm->nx; cc++)
           {
             accm->grid[cz][cw][cb][cr][cc] *= (f0 * normf[cz]);
           }
         }
         #endif // SPATIALLY_RESOLVED
       }
     }
   }
 }

/* Add value to a single cell of a backward accumulator: ***********************

 accm_b_add_k

 Add a fixed value to a single cell of the accumulator. This function is used to
 add the direct transmission contribution when the Sun is in the field-of-view
 of the sensor.

 INPUT:
 acm  - Accumulator parameters
        Pointer to accumulator_bmc struct;
 val  - Value to be added
        Constant double;
 cs   - Index of the Sun zenith angle layer
        Integer range: [0, N SZA], unitless
        Constant int;
 cr   - Index of the row for the spatially resolved accumulator
        Integer range: [0, N Y-axis bin], unitless
        Constant int;
 cc   - Index of the column for the spatially resolved accumulator
        Integer range: [0, N X-axis bin], unitless
        Constant int;

 OUTPUT:
 None. Updates the values pointed by accm.

*******************************************************************************/

 void 
 accm_b_add_k
 (
   struct accumulator_bmc * accm,
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

       #ifdef SPATIALLY_RESOLVED
       accm->grid[cs][cw][cb][cr][cc] += val;
       #endif // SPATIALLY_RESOLVED
     }
   }
 }

/* Sum backward accumulators: **************************************************

 accm_b_sum

 Sum accumulators. This function is used to calculate the average by summing 
 over the results of the subsets.

 INPUT:
 accm_sum - Accumulator receiving values;
            Pointer to accumulator_bmc struct;
 accm_in  - Accumulator providing values to be added;
            Pointer to accumulator_bmc struct;
 scale    - Multiplier (1/number of subdivisions)
            Constant double.

 OUTPUT:
 None. Updates the values pointed by accm_sum.

*******************************************************************************/

 void
 accm_b_sum
 (
   struct accumulator_bmc * accm_sum,
   struct accumulator_bmc const * accm_in,
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

         #ifdef SPATIALLY_RESOLVED
         for (size_t cr = 0; cr < accm_in->ny; cr++)
         {
           for (size_t cc = 0; cc < accm_in->nx; cc++)
           {
             accm_sum->grid[cz][cw][cb][cr][cc] += 
               accm_in->grid[cz][cw][cb][cr][cc] * scale;
           }
         }
         #endif // SPATIALLY_RESOLVED
       }
     }
   }
 }

/* Reset backward accumulator: *************************************************

 accm_b_reset

 Reset accumulator to 0.0.
 
 INPUT:
 acm - Accumulator parameters
       Pointer to accumulator_bmc struct;

 OUTPUT:
 None. Updates values of the accumulator.

*******************************************************************************/

 void
 accm_b_reset
 ( struct accumulator_bmc * accm )
 {
   for (size_t cz = 0; cz < accm->ns; cz++)
   {
     for (size_t cw = 0; cw < accm->nw0; cw++)
     {
       for (size_t cb = 0; cb < accm->nbr; cb++)
       {
         accm->vector[cz][cw][cb] = 0.0;

         #ifdef SPATIALLY_RESOLVED
         for (size_t cr = 0; cr < accm->ny; cr++)
         {
           for (size_t cc = 0; cc < accm->nx; cc++)
           {
             accm->grid[cz][cw][cb][cr][cc] = 0.0;
           }
         }
         #endif // SPATIALLY_RESOLVED
       }
     }
   }
 }

/* Write the grid or vector member of an accumulator: **************************

 accm_b_write_grid,  accm_b_write_vec

 Although the accumulators carry the number of illumination layers, the
 accumulators for the diffuse component has two additional layers by default 
 (uniform and cardioidal distributions). The number of Sun zenith angle is
 provided as an argument to make simpler find out the extra layers and set their
 name.

 INPUT:
 acm     - Accumulator parameters
           Pointer to accumulator_bmc struct;
 sim_ns  - Number of Sun zenith angles
 ofbn    - Base file name
           Pointer to constant char;
 sufx    - File name suffix
           Pointer to constant char;
 iop_w0  - Single scattering albedos
           Dimension: [N w0]
           Range: [0.0,1.0], unitless
           Pointer to constant double;
 btt_bhr - Bottom reflectance
           Dimension: [N BHR]
           Range: [0.0,1.0], unitless
           Pointer to constant double;
 sim_sza - Sun zenith angles
           Dimension: [N SZA]
           Range: [0.0,PI/2.0], radians
           Pointer to constant double;

 OUTPUT:
 None. Write values to disk.

*******************************************************************************/

 void
 accm_b_write_grid
 (
   struct accumulator_bmc * accm,
   int const sim_ns,
   char const * ofbn,
   char const * sufx,
   double const * iop_w0,
   double const * btt_bh,
   double const * sim_sza
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
   struct accumulator_bmc * accm,
   int const sim_ns,
   char const * ofbn,
   char const * sufx,
   double const * iop_w0,
   double const * btt_bhr,
   double const * sim_sza
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

